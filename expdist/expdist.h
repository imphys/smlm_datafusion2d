#pragma once

#include <stdint.h>
#include <cuda_runtime.h>

#ifdef WIN32
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT 
#endif

class DLL_EXPORT GPUExpDist {
    public:
        GPUExpDist(int max_n);
        ~GPUExpDist();
        double compute(const double *A, const double *B, int m, int n, const double *scale_A, const double *scale_B);
    private:
        int dim;
        int max_n;
        double *d_A;
        double *d_B;
        double *d_scale_A;
        double *d_scale_B;
        double *d_cross_term;

        cudaStream_t stream;
};


