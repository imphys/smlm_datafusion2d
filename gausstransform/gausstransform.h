#pragma once

#include <stdint.h>

#include <cuda_runtime.h>

class GPUGaussTransform {
    public:
        GPUGaussTransform(int max_n);
        ~GPUGaussTransform();
        double compute(const double *A, const double *B, int m, int n, double scale, double *grad);
    private:
        int dim;
        int max_n;
        double *d_A;
        double *d_B;
        double *d_grad;
        double *d_cross_term;

        cudaStream_t stream;
};


