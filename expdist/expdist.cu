/*
 * Host part for calling GPUExpDist from the CPU 
 */

#include <stdint.h>
#include <cuda_runtime.h>
#include <math.h>

#include "expdist.h"

//tuned for Nvidia K40
#ifndef block_size_x //if not using kernel tuner
#define block_size_x 32
#define block_size_y 4
#define tile_size_x 2
#define tile_size_y 4
#define use_shared_mem 1

#endif
#define reduce_block_size 256


#include "kernels.cu"


GPUExpDist::GPUExpDist(int n) {
    //allocate GPU memory for size max_n
    max_n = n;
    dim = 2;
    int elems = max_n * dim;

    cudaError_t err;

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
    err = cudaMalloc((void **)&d_scale_A, max_n*sizeof(double));
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaMalloc: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMalloc((void **)&d_scale_B, max_n*sizeof(double));
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

GPUExpDist::~GPUExpDist() {
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_scale_A);
    cudaFree(d_scale_B);
    cudaFree(d_cross_term);
    cudaStreamDestroy(stream);
} 

double GPUExpDist::compute(const double *A, const double *B, int m, int n, const double *scale_A, const double *scale_B) {

    double cost;
    cudaError_t err;

    //move data to the GPU
    err = cudaMemcpyAsync(d_A, A, m*dim*sizeof(double), cudaMemcpyHostToDevice, stream);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaMemcpyAsync d_A: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemcpyAsync(d_B, B, n*dim*sizeof(double), cudaMemcpyHostToDevice, stream);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaMemcpyAsync d_B: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemcpyAsync(d_scale_A, scale_A, m*sizeof(double), cudaMemcpyHostToDevice, stream);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaMemcpyAsync d_scale_A: %s\n", cudaGetErrorString(err));
        exit(1);
    }
    err = cudaMemcpyAsync(d_scale_B, scale_B, n*sizeof(double), cudaMemcpyHostToDevice, stream);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaMemcpyAsync d_scale_B: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    //compute number of thread blocks that would be used by the ExpDist kernel for this m and n
    int nblocks = ((int) ceil(m / (block_size_x*tile_size_x)) * (int) ceil(n / (block_size_y*tile_size_y)));

    //setup kernel execution parameters
    dim3 threads(block_size_x, block_size_y, 1);
    dim3 grid(1, 1, 1); //to be overwritten

    //check if the number of thread blocks does not exceed the allocated space
    //if it does, run the ExpDist_column kernel that uses fewer thread blocks
    if (nblocks < max_n) {
        //setup kernel execution parameters
        grid.x = (int) ceilf(m / (float)(block_size_x * tile_size_x));
        grid.y = (int) ceilf(n / (float)(block_size_y * tile_size_y));
    
        //call the first kernel
        ExpDist<<<grid, threads, 0, stream>>>(d_A, d_B, m, n, d_scale_A, d_scale_B, d_cross_term); 

    } else {
        //setup kernel execution parameters
        grid.x = (int) ceilf(m / (float)(block_size_x * tile_size_x));
    
        //call the first kernel
        ExpDist_column<<<grid, threads, 0, stream>>>(d_A, d_B, m, n, d_scale_A, d_scale_B, d_cross_term); 
    }

    //call the second kernel
    dim3 threads2(reduce_block_size, 1, 1);
    dim3 grid2(1, 1, 1);
    reduce_cross_term<<<grid2, threads2, 0, stream>>>(d_cross_term, d_cross_term, m, n, grid.x*grid.y);

    err = cudaMemcpyAsync(&cost, d_cross_term, 1*sizeof(double), cudaMemcpyDeviceToHost, stream);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in cudaMemcpyDeviceToHost cross_term: %s\n", cudaGetErrorString (err));
        exit(1);
    }

    return cost;
}


extern "C"
float test_GPUExpDistHost(double *cost, const double* A, const double* B,
            int m, int n, int dim, const double *scale_A, const double *scale_B, int max_n) {

    GPUExpDist gpu_expdist(max_n);

    *cost = gpu_expdist.compute(A, B, m, n, scale_A, scale_B);

    cudaDeviceSynchronize();
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "Error in test_GPUExpDistHost: %s\n", cudaGetErrorString (err));
        exit(1);
    }

    return 0.0;
}


