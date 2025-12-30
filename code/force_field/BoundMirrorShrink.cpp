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
 * B = BoundMirrorShrink(A)
 * This function shrinks the matrix to remove the padded mirror boundaries. 
 * It removes one pixel from all sides of the input matrix.
 */

#include "mex.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

using namespace std;

/**
 * Remove the padding from all sides of the matrix.
 */
void boundMirrorShrink(double *A, double *B, int m, int n)
{
    for (int i = 0; i < m - 2; i++)
    {
        for (int j = 0; j < n - 2; j++)
        {
            B[i + j * (m - 2)] = A[i + 1 + (j + 1) * m];
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Check for proper number of arguments */
    if (nrhs != 1)
    {
        mexErrMsgIdAndTxt("AGSM:BoundMirrorShrink:invalidNumInputs", "One input required.");
    }
    if (nlhs > 1)
    {
        mexErrMsgIdAndTxt("AGSM:BoundMirrorShrink:invalidNumOutputs", "At most one output.");
    }

    double *A = mxGetPr(prhs[0]);
    int m = (int)mxGetM(prhs[0]);
    int n = (int)mxGetN(prhs[0]);

    if (m < 4 || n < 4)
    {
        mexErrMsgIdAndTxt("AGSM:BoundMirrorShrink:matrixTooSmall", "Input matrix must be at least 4x4.");
    }

    plhs[0] = mxCreateDoubleMatrix(m - 2, n - 2, mxREAL);
    double *B = mxGetPr(plhs[0]);

    boundMirrorShrink(A, B, m, n);
}