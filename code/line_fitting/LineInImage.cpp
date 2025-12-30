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
 * [x, y] = LineInImage(m, n, theta, s)
 * This function generates points on a line within image boundaries.
 */

#include "mex.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

#define PI 3.14159265

using namespace std;

/**
 * Custom rounding function for compatibility.
 */
double custom_round(double a)
{
    int b = (int)a;
    a -= b;
    if (a >= 0.5) b++;
    return (double)b;
}

/**
 * Copy a coordinate array.
 */
void copyCoordinates(double *src, double *dst, int N)
{
    for (int i = 0; i < N; i++)
    {
        dst[i] = src[i];
    }
}

/**
 * Calculate coordinates of a line within (m, n) boundaries.
 */
int getLineCoordinates(int m, int n, double theta, double s, double *&x, double *&y)
{
    double theta_mod = theta;
    while (theta_mod > PI) theta_mod -= PI;
    while (theta_mod < 0) theta_mod += PI;
    int N = 0;

    if (theta_mod > PI / 4 && theta_mod < PI / 4 * 3) // More horizontal line
    {
        int x1, x2;
        double temp;
        for (x1 = 1; x1 <= n; x1++)
        {
            temp = custom_round((s - x1 * cos(theta)) / sin(theta));
            if (temp >= 1 && temp <= m) break;
        }
        if (x1 > n) return 0;

        for (x2 = n; x2 >= 1; x2--)
        {
            temp = custom_round((s - x2 * cos(theta)) / sin(theta));
            if (temp >= 1 && temp <= m) break;
        }

        N = x2 - x1 + 1;
        x = new double[N];
        y = new double[N];
        
        for (int i = x1; i <= x2; i++)
        {
            x[i - x1] = i;
            y[i - x1] = custom_round((s - i * cos(theta)) / sin(theta));
        }
    }
    else // More vertical line
    {
        int y1, y2;
        double temp;
        for (y1 = 1; y1 <= m; y1++)
        {
            temp = custom_round((s - y1 * sin(theta)) / cos(theta));
            if (temp >= 1 && temp <= n) break;
        }
        if (y1 > m) return 0;

        for (y2 = m; y2 >= 1; y2--)
        {
            temp = custom_round((s - y2 * sin(theta)) / cos(theta));
            if (temp >= 1 && temp <= n) break;
        }

        N = y2 - y1 + 1;
        x = new double[N];
        y = new double[N];
        
        for (int i = y1; i <= y2; i++)
        {
            y[i - y1] = i;
            x[i - y1] = custom_round((s - i * sin(theta)) / cos(theta));
        }
    }
    return N;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 4)
    {
        mexErrMsgIdAndTxt("AGSM:LineInImage:invalidNumInputs", "Four inputs required: m, n, theta, s.");
    }
    if (nlhs != 2)
    {
        mexErrMsgIdAndTxt("AGSM:LineInImage:invalidNumOutputs", "Two outputs required: x, y.");
    }

    int m = (int)mxGetScalar(prhs[0]);
    int n = (int)mxGetScalar(prhs[1]);
    double theta = mxGetScalar(prhs[2]);
    double s = mxGetScalar(prhs[3]);

    double *xx, *yy;
    int N = getLineCoordinates(m, n, theta, s, xx, yy);

    plhs[0] = mxCreateDoubleMatrix(1, N, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, N, mxREAL);

    double *outX = mxGetPr(plhs[0]);
    double *outY = mxGetPr(plhs[1]);

    if (N > 0)
    {
        copyCoordinates(xx, outX, N);
        copyCoordinates(yy, outY, N);
        delete[] xx;
        delete[] yy;
    }
}