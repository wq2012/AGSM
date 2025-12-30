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
 * B = BoundMirrorExpand(A)
 * This function expands the matrix using mirror boundary conditions. 
 * For example, if:
 * A = [ 1  2  3  11;
 *       4  5  6  12;
 *       7  8  9  13 ]
 * B = BoundMirrorExpand(A) will yield:
 *     [ 5  4  5  6  12  6;
 *       2  1  2  3  11  3;
 *       5  4  5  6  12  6;
 *       8  7  8  9  13  9;
 *       5  4  5  6  12  6 ]
 */

#include "mex.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

using namespace std;

/**
 * Expand the input matrix by one pixel on all sides using mirror reflection.
 */
void boundMirrorExpand(double *A, double *B, int m, int n)
{
    B[0 + 0 * (m + 2)] = A[1 + 1 * m];
    B[(m + 1) + 0 * (m + 2)] = A[(m - 2) + 1 * m];
    B[0 + (n + 1) * (m + 2)] = A[1 + (n - 2) * m];
    B[(m + 1) + (n + 1) * (m + 2)] = A[(m - 2) + (n - 2) * m];
    
    for (int i = 1; i < n + 1; i++)
    {
        B[0 + i * (m + 2)] = A[1 + (i - 1) * m];
        B[m + 1 + i * (m + 2)] = A[m - 2 + (i - 1) * m];
    }
    for (int i = 1; i < m + 1; i++)
    {
        B[i + 0 * (m + 2)] = A[(i - 1) + 1 * m];
        B[i + (n + 1) * (m + 2)] = A[i - 1 + (n - 2) * m];
    }
    for (int i = 1; i < m + 1; i++)
    {
        for (int j = 1; j < n + 1; j++)
        {
            B[i + j * (m + 2)] = A[i - 1 + (j - 1) * m];
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Check for proper number of arguments */
    if (nrhs != 1)
    {
        mexErrMsgIdAndTxt("AGSM:BoundMirrorExpand:invalidNumInputs", "One input required.");
    }
    if (nlhs > 1)
    {
        mexErrMsgIdAndTxt("AGSM:BoundMirrorExpand:invalidNumOutputs", "At most one output.");
    }

    double *A = mxGetPr(prhs[0]);
    int m = (int)mxGetM(prhs[0]);
    int n = (int)mxGetN(prhs[0]);

    if (m < 3 || n < 3)
    {
        mexErrMsgIdAndTxt("AGSM:BoundMirrorExpand:matrixTooSmall", "Input matrix must be at least 3x3.");
    }

    plhs[0] = mxCreateDoubleMatrix(m + 2, n + 2, mxREAL);
    double *B = mxGetPr(plhs[0]);

    boundMirrorExpand(A, B, m, n);
}