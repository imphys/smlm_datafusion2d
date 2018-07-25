#include <math.h>

//definition
#include "expdist_ref.c"



//use openmp just for measuring time
#include <omp.h>


//extern C interface
extern "C" {

float call_expdist(double *cost, const double *A, const double *B, int m, int n, int dim, const double *scale_A, const double *scale_B) {
    *cost = expdist(A, B, m, n, dim, scale_A, scale_B);
    return 0.0;
}


float time_expdist(double *cost, const double *A, const double *B, int m, int n, int dim, const double *scale_A, const double *scale_B) {
    double start = omp_get_wtime();
    *cost = expdist(A, B, m, n, dim, scale_A, scale_B);
    return (float)((omp_get_wtime() - start)*1e3);
}

}

