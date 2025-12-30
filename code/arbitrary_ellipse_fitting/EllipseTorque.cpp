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
 * torque = EllipseTorque(u, v, x, y, theta, xc, yc, phi)
 * This function calculates the torque on an ellipse about its center.
 */

#include "mex.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

using namespace std;

/**
 * Calculate the average torque along the ellipse contour about (xc, yc).
 */
double getTorque(double *u, double *v, double *x, double *y, double *theta, int m, int n, int N, double xc, double yc, double phi)
{
    double torque = 0;
    int r, c;
    double fx, fy;
    double d;
    
    for (int i = 0; i < N; i++)
    {
        r = (int)y[i] - 1;
        c = (int)x[i] - 1;
        
        if (r >= 0 && r < m && c >= 0 && c < n)
        {
            fx = u[r + c * m];
            fy = v[r + c * m];
            d = sqrt((x[i] - xc) * (x[i] - xc) + (y[i] - yc) * (y[i] - yc));
            
            torque += ((-fx * sin(theta[i] + phi) + fy * cos(theta[i] + phi)) * d);
        }
    }
    return N > 0 ? torque / (N * (double)N) : 0;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Check for proper number of arguments */
    if (nrhs != 8)
    {
        mexErrMsgIdAndTxt("AGSM:EllipseTorque:invalidNumInputs", "Eight inputs required: u, v, x, y, theta, xc, yc, phi.");
    }
    if (nlhs > 1)
    {
        mexErrMsgIdAndTxt("AGSM:EllipseTorque:invalidNumOutputs", "At most one output.");
    }

    double *u = mxGetPr(prhs[0]);
    double *v = mxGetPr(prhs[1]);
    int m = (int)mxGetM(prhs[0]);
    int n = (int)mxGetN(prhs[0]);

    if (m != mxGetM(prhs[1]) || n != mxGetN(prhs[1]))
    {
        mexErrMsgIdAndTxt("AGSM:EllipseTorque:invalidInputDimension", "Input u and v must have same dimensions.");
    }

    double *x = mxGetPr(prhs[2]);
    double *y = mxGetPr(prhs[3]);
    double *theta = mxGetPr(prhs[4]);
    int N = (int)mxGetN(prhs[2]);

    if (mxGetM(prhs[2]) != 1 || mxGetM(prhs[3]) != 1 || mxGetM(prhs[4]) != 1 || 
        mxGetN(prhs[3]) != N || mxGetN(prhs[4]) != N)
    {
        mexErrMsgIdAndTxt("AGSM:EllipseTorque:invalidInputDimension", 
            "Input x, y, and theta should have same dimensions 1xN.");
    }

    double xc = mxGetScalar(prhs[5]);
    double yc = mxGetScalar(prhs[6]);
    double phi = mxGetScalar(prhs[7]);

    double torque = getTorque(u, v, x, y, theta, m, n, N, xc, yc, phi);

    plhs[0] = mxCreateDoubleScalar(torque);
}