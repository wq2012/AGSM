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
 * [xc, yc, r] = InitialCircle(f)
 * This function computes initial guesses for circle parameters (center and radius) 
 * given a 2D intensity or energy field.
 */

#include "mex.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

using namespace std;

/**
 * Core logic for initial guess based on center of gravity and average distance.
 */
void guess(double *I, int m, int n, double *xc_, double *yc_, double *r_)
{
    double xc = 0;
    double yc = 0;
    double r = 0;
    double sumI = 0;
    
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            xc += (j * I[i + j * m]);
            yc += (i * I[i + j * m]);
            sumI += I[i + j * m];
        }
    }
    
    if (sumI != 0)
    {
        xc /= sumI;
        yc /= sumI;
    }
    else
    {
        xc = n / 2.0;
        yc = m / 2.0;
    }
    
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            r += (sqrt((i - yc) * (i - yc) + (j - xc) * (j - xc)) * I[i + j * m]);
        }
    }
    
    if (sumI != 0) r /= sumI;
    else r = min(m, n) / 4.0;
    
    *xc_ = xc;
    *yc_ = yc;
    *r_ = r;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Check for proper number of arguments */
    if (nrhs != 1)
    {
        mexErrMsgIdAndTxt("AGSM:InitialCircle:invalidNumInputs", "One input required.");
    }
    if (nlhs != 3)
    {
        mexErrMsgIdAndTxt("AGSM:InitialCircle:invalidNumOutputs", "Three outputs required.");
    }

    double *f = mxGetPr(prhs[0]);
    int m = (int)mxGetM(prhs[0]);
    int n = (int)mxGetN(prhs[0]);

    if (m < 5 || n < 5)
    {
        mexErrMsgIdAndTxt("AGSM:InitialCircle:matrixTooSmall", "Input matrix is too small.");
    }

    plhs[0] = mxCreateDoubleScalar(0);
    plhs[1] = mxCreateDoubleScalar(0);
    plhs[2] = mxCreateDoubleScalar(0);
    
    double *xc = mxGetPr(plhs[0]);
    double *yc = mxGetPr(plhs[1]);
    double *r = mxGetPr(plhs[2]);

    guess(f, m, n, xc, yc, r);
}
