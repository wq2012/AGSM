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
 * force = CircleNormalForce(u, v, x, y, theta)
 * This function calculates the normal force on a circle.
 */

#include "mex.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

using namespace std;

/**
 * Calculate the average normal force along the given points of a circle.
 */
double getNormalForce(double *u, double *v, double *x, double *y, double *theta, int m, int n, int N)
{
    double force = 0;
    int r, c;
    double fx, fy;
    
    for (int i = 0; i < N; i++)
    {
        r = (int)y[i] - 1;
        c = (int)x[i] - 1;
        
        if (r >= 0 && r < m && c >= 0 && c < n)
        {
            fx = u[r + c * m];
            fy = v[r + c * m];
            force += (fx * cos(theta[i]) + fy * sin(theta[i]));
        }
    }
    return N > 0 ? force / N : 0;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Check for proper number of arguments */
    if (nrhs != 5)
    {
        mexErrMsgIdAndTxt("AGSM:CircleNormalForce:invalidNumInputs", "Five inputs required: u, v, x, y, theta.");
    }
    if (nlhs > 1)
    {
        mexErrMsgIdAndTxt("AGSM:CircleNormalForce:invalidNumOutputs", "At most one output.");
    }

    double *u = mxGetPr(prhs[0]);
    double *v = mxGetPr(prhs[1]);
    int m = (int)mxGetM(prhs[0]);
    int n = (int)mxGetN(prhs[0]);

    if (m != mxGetM(prhs[1]) || n != mxGetN(prhs[1]))
    {
        mexErrMsgIdAndTxt("AGSM:CircleNormalForce:invalidInputDimension", "Input u and v must have same dimensions.");
    }

    double *x = mxGetPr(prhs[2]);
    double *y = mxGetPr(prhs[3]);
    double *theta = mxGetPr(prhs[4]);
    int N = (int)mxGetN(prhs[2]);

    if (mxGetM(prhs[2]) != 1 || mxGetM(prhs[3]) != 1 || mxGetM(prhs[4]) != 1 || 
        mxGetN(prhs[3]) != N || mxGetN(prhs[4]) != N)
    {
        mexErrMsgIdAndTxt("AGSM:CircleNormalForce:invalidInputDimension", 
            "Input x, y, and theta should have same dimensions 1xN.");
    }

    double force = getNormalForce(u, v, x, y, theta, m, n, N);

    plhs[0] = mxCreateDoubleScalar(force);
}