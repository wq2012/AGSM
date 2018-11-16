/**
 * Copyright (C) 2012 Quan Wang <wangq10@rpi.edu>,
 * Signal Analysis and Machine Perception Laboratory,
 * Department of Electrical, Computer, and Systems Engineering,
 * Rensselaer Polytechnic Institute, Troy, NY 12180, USA
 *
 * Quan Wang, Kim L. Boyer, 
 * The active geometric shape model: A new robust deformable shape model and its applications, 
 * Computer Vision and Image Understanding, Volume 116, Issue 12, December 2012, Pages 1178-1194, 
 * ISSN 1077-3142, 10.1016/j.cviu.2012.08.004. 
 */

/**
 * This function ensures mirror boundary condition. 
 * For example, 
 * 
 * A = [
 *     X  X  X  X  X   X
 *     X  1  2  3  11  X
 *     X  4  5  6  12  X 
 *     X  7  8  9  13  X 
 *     X  X  X  X  X   X
 *     ]
 *
 * B = BoundMirrorEnsure(A) will yield
 *
 *     5  4  5  6  12  6
 *     2  1  2  3  11  3
 *     5  4  5  6  12  6 
 *     8  7  8  9  13  9 
 *     5  4  5  6  12  6
 */
        
#include "mex.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

void boundMirrorEnsure(double *A, int m, int n)
{
    A[0+0*m]=A[2+2*m];
    A[m-1+0*m]=A[m-3+2*m];
    A[0+(n-1)*m]=A[2+(n-3)*m];
    A[m-1+(n-1)*m]=A[(m-3)+(n-3)*m];
    for(int i=1;i<n-1;i++)
    {
        A[0+i*m]=A[2+i*m];
        A[m-1+i*m]=A[m-3+i*m];
    }
    for(int i=1;i<m-1;i++)
    {
        A[i+0*m]=A[i+2*m];
        A[i+(n-1)*m]=A[i+(n-3)*m];
    }
}

void copyMatrix(double *A, double *B, int m, int n)
{
    for(int i=0;i<m*n;i++)
    {
        B[i]=A[i];
    }
}

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    if(nrhs!=1)
    {
        mexErrMsgIdAndTxt( "MATLAB:BoundMirrorEnsure:invalidNumInputs",
                "One input required.");
    }
    if(nlhs!=1)
    {
        mexErrMsgIdAndTxt( "MATLAB:BoundMirrorEnsure:invalidNumOutputs",
                "One output required.");
    }
    
    double *A=mxGetPr(prhs[0]);
    int m=mxGetM(prhs[0]);
    int n=mxGetN(prhs[0]);
    
    if(m<4 || n<4)
    {
        mexErrMsgIdAndTxt( "MATLAB:BoundMirrorEnsure:matrixTooSmall",
                "Input matrix is too small.");
    }
    
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    
    double *B=mxGetPr(plhs[0]);
    
    boundMirrorEnsure(A,m,n);
    
    copyMatrix(A, B, m, n);

}