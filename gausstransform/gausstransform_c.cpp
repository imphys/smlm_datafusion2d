#include <math.h>

//definition
#include "gausstransform_ref.cpp"



//use openmp just for measuring time
#include <omp.h>


//extern C interface
extern "C" {

float call_GaussTransform(double* cost, const double* A, const double* B, int m, int n, int dim, double 
        scale, double* grad) {
    *cost = GaussTransform(A, B, m, n, dim, scale, grad);
    return 0.0;
}


float time_GaussTransform(double* cost, const double* A, const double* B, int m, int n, int dim, double 
        scale, double* grad) {
    double start = omp_get_wtime();
    *cost = GaussTransform(A, B, m, n, dim, scale, grad);
    return (float)((omp_get_wtime() - start)*1e3);
}

}

