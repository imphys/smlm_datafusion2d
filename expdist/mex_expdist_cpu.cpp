#include <stdio.h>
#include <cmath>

/*%%=====================================================================
 * %% Project:   Pointset Registration using Gaussian Mixture Model
 * %% Module:    $RCSfile: mex_GaussTransform.c,v $
 * %% Language:  C
 * %% Author:    $Author: bing.jian $
 * %% Date:      $Date: 2008-11-13 16:34:29 -0500 (Thu, 13 Nov 2008) $
 * %% Version:   $Revision: 109 $
 * %%=====================================================================*/

#include "mex.h"

#include "expdist_ref.c"

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
    /* Declare variables */ 
    int m, n, dim; 
    double *A, *B, *result, *scale_A, *scale_B;

    //prevent clearing from memory
    void mexLock(void);
    
    /* Check for proper number of input and output arguments */    
    if (nrhs != 4) {
	mexErrMsgTxt("Three input arguments required.");
    } 
    if (nlhs > 2){
	mexErrMsgTxt("Too many output arguments.");
    }
    
    /* Check data type of input argument */
    if (!(mxIsDouble(prhs[0]))) {
      mexErrMsgTxt("Input array must be of type double.");
    }
    if (!(mxIsDouble(prhs[1]))) {
      mexErrMsgTxt("Input array must be of type double.");
    }
    if (!(mxIsDouble(prhs[2]))) {
      mexErrMsgTxt("Input array must be of type double.");
    }
    if (!(mxIsDouble(prhs[3]))) {
      mexErrMsgTxt("Input array must be of type double.");
    }     
     
    /* Get the number of elements in the input argument */
    /* elements=mxGetNumberOfElements(prhs[0]); */
    /* Get the data */
    A = (double *)mxGetPr(prhs[0]);
    B = (double *)mxGetPr(prhs[1]);
    scale_A = (double *)mxGetPr(prhs[2]);
    scale_B = (double *)mxGetPr(prhs[3]);
  	/* Get the dimensions of the matrix input A&B. */
  	m = mxGetM(prhs[0]);
  	n = mxGetM(prhs[1]);

  	dim = mxGetN(prhs[0]);
  	if (mxGetN(prhs[1])!=dim)
  	{
  		mexErrMsgTxt("The two input point sets should have same dimension.");
  	}
    
  	if (mxGetM(prhs[0])!=mxGetM(prhs[2]))
  	{
  		mexErrMsgTxt("localizations and uncertainties (A) should have same dimension.");
  	}    
    
  	if (mxGetM(prhs[1])!=mxGetM(prhs[3]))
  	{
  		mexErrMsgTxt("localizations and uncertainties (B) should have same dimension.");
  	}        
    
  	if (mxGetN(prhs[2])!=1)
  	{
  		mexErrMsgTxt("uncertainties should be a m*1 matrix.");
  	}        
    
  	if (mxGetN(prhs[3])!=1)
  	{
  		mexErrMsgTxt("uncertainties should be a n*1 matrix.");
  	}            
    /* Allocate the space for the return argument */

    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    result = mxGetPr(plhs[0]);

    *result = expdist(A, B, m, n, dim, scale_A, scale_B);

}

