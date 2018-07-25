/*
 * Host part for calling GPUGaussTransform from the CPU 
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <cuda_runtime.h>

#include "gausstransform.h"

#ifndef block_size_x
    #define block_size_x 128    //best for GTX 690 and K40
#endif

#include "kernels.cu"



GPUGaussTransform::GPUGaussTransform(int n) {
    //allocate GPU memory for size max_n
    max_n = n;
    dim = 2;
    int elems = max_n * dim;

    cudaError_t err;

    //pseudo load balancing across available GPUs
    int count;
    cudaGetDeviceCount(&count);
    int id = rand() % count;
    cudaSetDevice(id);

    err = cudaMalloc((void **)&d_A, elems*sizeof(double));
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaMalloc: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMalloc((void **)&d_B, elems*sizeof(double));
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaMalloc: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMalloc((void **)&d_grad, elems*sizeof(double));
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaMalloc: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMalloc((void **)&d_cross_term, max_n*sizeof(double));
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaMalloc: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    err = cudaStreamCreate(&stream);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaStreamCreate: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    cudaDeviceSynchronize();
} 


GPUGaussTransform::~GPUGaussTransform() {
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_grad);
    cudaFree(d_cross_term);
    cudaStreamDestroy(stream);
} 

double GPUGaussTransform::compute(const double *A, const double *B,
    int m, int n, double scale, double *grad) {

    double energy;
    cudaError_t err;

    //move data to the GPU
    err = cudaMemcpyAsync(d_A, A, m*dim*sizeof(double), cudaMemcpyHostToDevice, stream);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaMemcpyAsync: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemcpyAsync(d_B, B, n*dim*sizeof(double), cudaMemcpyHostToDevice, stream);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaMemcpyAsync: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    //setup kernel execution parameters
    dim3 threads(block_size_x, 1, 1);
    dim3 grid(m, 1, 1);
    
    //call the first kernel
    double scale_sq = scale * scale;
    GaussTransform<<<grid, threads, 0, stream>>>(d_A, d_B, m, n, scale_sq, d_grad, d_cross_term); 

    //call the second kernel
    dim3 grid2(1, 1, 1);
    reduce_cross_term<<<grid2, threads, 0, stream>>>(d_cross_term, d_cross_term, m, n, m);

    //copy result from GPU memory to host memory
    err = cudaMemcpyAsync(grad, d_grad, m*dim*sizeof(double), cudaMemcpyDeviceToHost, stream);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaMemcpyDeviceToHost: %s\n", cudaGetErrorString (err));
        exit(1);
    }

    err = cudaMemcpyAsync(&energy, d_cross_term, 1*sizeof(double), cudaMemcpyDeviceToHost, stream);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaMemcpyDeviceToHost: %s\n", cudaGetErrorString (err));
        exit(1);
    }

    //wait for the GPU stuff to have finished
    cudaStreamSynchronize(stream);

    return energy;
}



extern "C"
float test_GaussTransformHost(double *cost, const double* A, const double* B,
            int m, int n, int dim, double scale, double* grad) {

    GPUGaussTransform gpu_gt(1000000);

    *cost = gpu_gt.compute(A, B, m, n, scale, grad);

    return 0.0;
}


