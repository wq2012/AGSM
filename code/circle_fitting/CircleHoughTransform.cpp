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
 * H = CircleHoughTransform(I, rmin, rmax, P)
 * This function performs a 3D Circle Hough Transform on an image.
 */

#include "mex.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

using namespace std;

/**
 * Custom rounding function.
 */
double custom_round(double a)
{
    if (a - floor(a) > 0.5) return floor(a) + 1;
    else return floor(a);
}

/**
 * Calculate the flat index for the 3D Hough accumulator.
 */
long H_index(long y, long x, long r, const mwSize *dim)
{
    long i = (y - 1) + (x - 1) * (long)dim[0] + (r - 1) * (long)dim[0] * (long)dim[1];
    if (i >= (long)dim[0] * (long)dim[1] * (long)dim[2])
    {
        return -1; // Out of bounds
    }
    return i;
}

/**
 * Core Hough Transform logic.
 */
void hough(double *I, int m, int n, double rmin, double rmax, const mwSize* dim, double P, double *H)
{
    mexPrintf("Running Hough transform ...\n");
    
    for (int i = 1; i <= m; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            if (I[(i - 1) + (j - 1) * m] < 0.5) continue;
            
            for (int x_index = 1; x_index <= (int)dim[1]; x_index++)
            {
                double xx = (x_index - 1) / P + 1;
                for (int y_index = 1; y_index <= (int)dim[0]; y_index++)
                {
                    double yy = (y_index - 1) / P + 1;
                    double r = sqrt((xx - j) * (xx - j) + (yy - i) * (yy - i));
                    r = custom_round(r * P) / P;
                    
                    if (r < rmin || r > rmax) continue;
                    
                    int r_index = (int)((r - rmin) * P + 1.5); // Adjusted for stability
                    if (r_index > (int)dim[2]) continue;
                    
                    long idx = H_index((long)y_index, (long)x_index, (long)r_index, dim);
                    if (idx != -1) H[idx]++;
                }
            }
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Check input */
    if (nrhs != 4)
    {
        mexErrMsgIdAndTxt("AGSM:CircleHoughTransform:invalidNumInputs", "Four inputs required: I, rmin, rmax, P.");
    }
    if (nlhs != 1)
    {
        mexErrMsgIdAndTxt("AGSM:CircleHoughTransform:invalidNumOutputs", "One output required: H.");
    }
    
    double *I = mxGetPr(prhs[0]);
    int m = (int)mxGetM(prhs[0]);
    int n = (int)mxGetN(prhs[0]);
    
    if (m < 5 || n < 5)
    {
        mexErrMsgIdAndTxt("AGSM:CircleHoughTransform:matrixTooSmall", "Input matrix is too small.");
    }
    
    double rmin = mxGetScalar(prhs[1]);
    double rmax = mxGetScalar(prhs[2]);
    double P = mxGetScalar(prhs[3]);
    
    /* Set output dimensions */
    mwSize *dim = new mwSize[3];
    dim[0] = (mwSize)((m - 1) * P + 1);
    dim[1] = (mwSize)((n - 1) * P + 1);
    dim[2] = (mwSize)((rmax - rmin) * P + 1);
    
    mexPrintf("Circle Hough transform, size of volume: %d * %d * %d\n", dim[0], dim[1], dim[2]);
    plhs[0] = mxCreateNumericArray(3, dim, mxDOUBLE_CLASS, mxREAL);
    double *H = mxGetPr(plhs[0]);
    
    /* Run Hough algorithm */
    hough(I, m, n, rmin, rmax, dim, P, H);
    
    delete[] dim;
}