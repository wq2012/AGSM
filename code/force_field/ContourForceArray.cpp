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
 * F = ContourForceArray(u, v, x, y)
 * This function returns the force vectors at specified contour points.
 */

#include "mex.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

using namespace std;

/**
 * Extract force field values at specified contour coordinates.
 */
void getForceArray(double *u, double *v, double *x, double *y, int m, int n, int N, double *F)
{
    int r, c;
    for (int i = 0; i < N; i++)
    {
        r = (int)y[i] - 1;
        c = (int)x[i] - 1;
        
        if (r >= 0 && r < m && c >= 0 && c < n)
        {
            F[i] = u[r + c * m];
            F[i + N] = v[r + c * m];
        }
        else
        {
            F[i] = 0;
            F[i + N] = 0;
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Check for proper number of arguments */
    if (nrhs != 4)
    {
        mexErrMsgIdAndTxt("AGSM:ContourForceArray:invalidNumInputs", "Four inputs required: u, v, x, y.");
    }
    if (nlhs > 1)
    {
        mexErrMsgIdAndTxt("AGSM:ContourForceArray:invalidNumOutputs", "At most one output.");
    }

    double *u = mxGetPr(prhs[0]);
    double *v = mxGetPr(prhs[1]);
    int m = (int)mxGetM(prhs[0]);
    int n = (int)mxGetN(prhs[0]);

    if (m != mxGetM(prhs[1]) || n != mxGetN(prhs[1]))
    {
        mexErrMsgIdAndTxt("AGSM:ContourForceArray:invalidInputDimension", "Input u and v must have same dimensions.");
    }

    double *x = mxGetPr(prhs[2]);
    double *y = mxGetPr(prhs[3]);
    int N = (int)mxGetN(prhs[2]);

    if (mxGetM(prhs[2]) != 1 || mxGetM(prhs[3]) != 1 || mxGetN(prhs[3]) != N)
    {
        mexErrMsgIdAndTxt("AGSM:ContourForceArray:invalidInputDimension", "Input x and y must have same dimensions 1xN.");
    }

    plhs[0] = mxCreateDoubleMatrix(N, 2, mxREAL);
    double *F = mxGetPr(plhs[0]);

    getForceArray(u, v, x, y, m, n, N, F);
}