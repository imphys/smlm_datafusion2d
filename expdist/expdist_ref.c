/*
    expdist   Computes the bhattacharya cost function for two given point set

    SYNOPSIS:
    D = expdist(A, B, m, n, dim, scale_A, scale_B);

    INPUT
        A
            The first particle containing the list of coordinates
        B
            The second particle containing the list of coordinates
        m
            Size of the particle A
        n
            Size of the particle B
        dim
            particles dimension (2D or 3D) 
        scale_A
            uncertainties of particle A 
        scale_B
            uncertainties of particle B 

    OUTPUT
        result
            The distance between particle A and B

    (C) Copyright 2017              Quantitative Imaging Group
        All rights reserved         Faculty of Applied Physics
                                    Delft University of Technology
                                    Lorentzweg 1
                                    2628 CJ Delft
                                    The Netherlands
    Hamidreza Heydarian, Feb 2017
*/

#define SQR(X)  ((X)*(X))

#include <math.h>
/* #include <stdio.h> */

#ifdef WIN32
__declspec( dllexport )
#endif
double expdist(const double* A, const double* B,  int m, int n, int dim, const double* scale_A, const double* scale_B)
{
    int i,j,d;
    int id, jd;
    double dist_ij, cross_term = 0;

    for (i=0;i<m;++i)
    {
        for (j=0;j<n;++j)
        {
            dist_ij = 0;
            for (d=0;d<dim;++d)
            {
                id = i + d * m;
                jd = j + d * n;
                dist_ij = dist_ij + SQR( A[id] - B[jd]);
            }
            cross_term += exp(-dist_ij/(scale_A[i] + scale_B[j]));
        }
    }

    return cross_term;
}

