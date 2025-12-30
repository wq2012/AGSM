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
 * [torque, k] = LineTorque(u, v, x, y, theta)
 * This function calculates the maximum torque and the pivot point index.
 */

#include "mex.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

using namespace std;

/**
 * Calculate the maximum torque along the line and find the best pivot index k.
 */
double getTorque(double *u, double *v, double *x, double *y, int m, int n, int N, double theta, int &k)
{
    double maxTorque = 0;
    double torque;
    int bestK = 1;
    double cos_theta = cos(theta);
    double sin_theta = sin(theta);
    
    for (int curK = 1; curK <= N; curK++)
    {
        torque = 0;
        for (int i = 0; i < N; i++)
        {
            int r = (int)y[i] - 1;
            int c = (int)x[i] - 1;
            
            if (r >= 0 && r < m && c >= 0 && c < n)
            {
                double fx = u[r + c * m];
                double fy = v[r + c * m];
                double d = sqrt((x[i] - x[curK - 1]) * (x[i] - x[curK - 1]) + 
                                (y[i] - y[curK - 1]) * (y[i] - y[curK - 1]));
                
                int sgn = 0;
                if (curK - 1 == i) sgn = 0;
                else if (x[curK - 1] < x[(i > curK - 1 ? N - 1 : 0)]) sgn = -1;
                else sgn = 1;

                torque += ((fx * cos_theta + fy * sin_theta) * d * sgn);
            }
        }
        torque /= (N * (double)N);
        if (fabs(torque) > fabs(maxTorque))
        {
            maxTorque = torque;
            bestK = curK;
        }
    }
    
    k = bestK;
    return maxTorque;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Check for proper number of arguments */
    if (nrhs != 5)
    {
        mexErrMsgIdAndTxt("AGSM:LineTorque:invalidNumInputs", "Five inputs required: u, v, x, y, theta.");
    }
    if (nlhs > 2)
    {
        mexErrMsgIdAndTxt("AGSM:LineTorque:invalidNumOutputs", "At most two outputs: torque, k.");
    }

    double *u = mxGetPr(prhs[0]);
    double *v = mxGetPr(prhs[1]);
    int m = (int)mxGetM(prhs[0]);
    int n = (int)mxGetN(prhs[0]);

    if (m != mxGetM(prhs[1]) || n != mxGetN(prhs[1]))
    {
        mexErrMsgIdAndTxt("AGSM:LineTorque:invalidInputDimension", "Input u and v must have same dimensions.");
    }

    double *x = mxGetPr(prhs[2]);
    double *y = mxGetPr(prhs[3]);
    int N = (int)mxGetN(prhs[2]);

    if (mxGetM(prhs[2]) != 1 || mxGetM(prhs[3]) != 1 || mxGetN(prhs[3]) != N)
    {
        mexErrMsgIdAndTxt("AGSM:LineTorque:invalidInputDimension", "Input x and y must be 1xN vectors of same size.");
    }

    if (!mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) || mxGetN(prhs[4]) * mxGetM(prhs[4]) != 1)
    {
        mexErrMsgIdAndTxt("AGSM:LineTorque:thetaNotScalar", "Input theta must be a scalar.");
    }
    double theta = mxGetScalar(prhs[4]);

    int k;
    double torque = getTorque(u, v, x, y, m, n, N, theta, k);

    plhs[0] = mxCreateDoubleScalar(torque);
    if (nlhs > 1)
    {
        plhs[1] = mxCreateDoubleScalar((double)k);
    }
}