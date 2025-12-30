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
 * [u, v] = GVF(f, alpha, mu, iter)
 * This function computes the Gradient Vector Flow (GVF) field.
 *
 * Inputs:
 *    f - 2D energy field (or gradient field).
 *    alpha - Step size of update, between 0 and 1.
 *    mu - GVF regularization coefficient.
 *    iter - Number of iterations.
 *
 * Outputs:
 *    u - Resulting GVF in x direction.
 *    v - Resulting GVF in y direction.
 */

#include "mex.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

using namespace std;

/**
 * Expand matrix with mirror boundary conditions. 
 */
double *boundMirrorExpand(double *A, int m, int n)
{
    double *B = new double[(m + 2) * (n + 2)];
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
    return B;
}

/**
 * Shrink matrix to remove padding.
 */
double *boundMirrorShrink(double *A, int m, int n)
{
    double *B = new double[(m - 2) * (n - 2)];
    for (int i = 0; i < m - 2; i++)
    {
        for (int j = 0; j < n - 2; j++)
        {
            B[i + j * (m - 2)] = A[i + 1 + (j + 1) * m];
        }
    }
    return B;
}

/**
 * Ensure mirror reflection at boundaries during iteration. 
 */
void boundMirrorEnsure(double *A, int m, int n)
{
    A[0 + 0 * m] = A[2 + 2 * m];
    A[m - 1 + 0 * m] = A[m - 3 + 2 * m];
    A[0 + (n - 1) * m] = A[2 + (n - 3) * m];
    A[m - 1 + (n - 1) * m] = A[(m - 3) + (n - 3) * m];
    
    for (int i = 1; i < n - 1; i++)
    {
        A[0 + i * m] = A[2 + i * m];
        A[m - 1 + i * m] = A[m - 3 + i * m];
    }
    for (int i = 1; i < m - 1; i++)
    {
        A[i + 0 * m] = A[i + 2 * m];
        A[i + (n - 1) * m] = A[i + (n - 3) * m];
    }
}

/**
 * General matrix copy.
 */
void copyMatrix(double *A, double *B, int m, int n)
{
    for (int i = 0; i < m * n; i++)
    {
        B[i] = A[i];
    }
}

/**
 * Compute central difference gradient.
 */
void gradient(double *f, double *fx, double *fy, int m, int n)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (j == 0) fx[i + j * m] = f[i + (j + 1) * m] - f[i + j * m];
            else if (j == n - 1) fx[i + j * m] = f[i + j * m] - f[i + (j - 1) * m];
            else fx[i + j * m] = (f[i + (j + 1) * m] - f[i + (j - 1) * m]) / 2.0;
            
            if (i == 0) fy[i + j * m] = f[i + 1 + j * m] - f[i + j * m];
            else if (i == m - 1) fy[i + j * m] = f[i + j * m] - f[i - 1 + j * m];
            else fy[i + j * m] = (f[i + 1 + j * m] - f[i - 1 + j * m]) / 2.0;
        }
    }
}

/**
 * Compute Laplacian operator.
 */
void del2(double *f, double *del2f, int m, int n)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == 0 || i == m - 1 || j == 0 || j == n - 1) del2f[i + j * m] = 0;
            else del2f[i + j * m] = f[i + 1 + j * m] + f[i - 1 + j * m] + f[i + (j + 1) * m] + f[i + (j - 1) * m] - 4 * f[i + j * m];
        }
    }
}

/**
 * Main GVF iterative computation.
 */
void gvf(double *f, int m, int n, double alpha, double mu, int iter, double *u, double *v)
{
    double fmin = f[0], fmax = f[0];
    for (int i = 1; i < m * n; i++)
    {
        if (f[i] > fmax) fmax = f[i];
        if (f[i] < fmin) fmin = f[i];
    }
    
    if (fmax <= fmin)
    {
        mexErrMsgIdAndTxt("AGSM:GVF:invalidField", "Input field f is constant or degenerate.");
    }
    
    for (int i = 0; i < m * n; i++) f[i] = (f[i] - fmin) / (fmax - fmin);

    double *f_exp = boundMirrorExpand(f, m, n);
    int m2 = m + 2, n2 = n + 2;
    double *fx = new double[m2 * n2];
    double *fy = new double[m2 * n2];
    gradient(f_exp, fx, fy, m2, n2);
    
    double *u_int = new double[m2 * n2];
    double *v_int = new double[m2 * n2];
    copyMatrix(fx, u_int, m2, n2);
    copyMatrix(fy, v_int, m2, n2);
    
    double *sqr_mag = new double[m2 * n2];
    for (int i = 0; i < m2 * n2; i++) sqr_mag[i] = fx[i] * fx[i] + fy[i] * fy[i];
    
    double *del_u = new double[m2 * n2];
    double *del_v = new double[m2 * n2];
    
    for (int it = 0; it < iter; it++)
    {
        boundMirrorEnsure(u_int, m2, n2);
        boundMirrorEnsure(v_int, m2, n2);
        del2(u_int, del_u, m2, n2);
        del2(v_int, del_v, m2, n2);
        
        for (int i = 0; i < m2 * n2; i++)
        {
            u_int[i] += alpha * (mu * del_u[i] - sqr_mag[i] * (u_int[i] - fx[i]));
            v_int[i] += alpha * (mu * del_v[i] - sqr_mag[i] * (v_int[i] - fy[i]));
        }
    }
    
    double *u_shrk = boundMirrorShrink(u_int, m2, n2);
    double *v_shrk = boundMirrorShrink(v_int, m2, n2);
    copyMatrix(u_shrk, u, m, n);
    copyMatrix(v_shrk, v, m, n);
   
    delete[] f_exp; delete[] fx; delete[] fy;
    delete[] u_int; delete[] v_int;
    delete[] del_u; delete[] del_v;
    delete[] u_shrk; delete[] v_shrk;
    delete[] sqr_mag;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Check input arguments */
    if (nrhs != 4)
    {
        mexErrMsgIdAndTxt("AGSM:GVF:invalidNumInputs", "Four inputs required: f, alpha, mu, iter.");
    }
    if (nlhs != 2)
    {
        mexErrMsgIdAndTxt("AGSM:GVF:invalidNumOutputs", "Two outputs required: [u, v].");
    }
    
    double *f_in = mxGetPr(prhs[0]);
    int m = (int)mxGetM(prhs[0]);
    int n = (int)mxGetN(prhs[0]);
    
    if (m < 10 || n < 10)
    {
        mexErrMsgIdAndTxt("AGSM:GVF:fTooSmall", "Input matrix f is too small.");
    }
    
    double *f = new double[m * n];
    for (int i = 0; i < m * n; i++) f[i] = f_in[i];
    
    double alpha = mxGetScalar(prhs[1]);
    if (alpha <= 0 || alpha >= 1.0)
    {
        mexWarnMsgIdAndTxt("AGSM:GVF:invalidAlpha", "Alpha should be in range (0, 1). Using default if necessary.");
    }
    
    double mu = mxGetScalar(prhs[2]);
    int iter = (int)mxGetScalar(prhs[3]);
    
    /* Create results */
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(m, n, mxREAL);
    double *u = mxGetPr(plhs[0]);
    double *v = mxGetPr(plhs[1]);
    
    /* Compute GVF */
    gvf(f, m, n, alpha, mu, iter, u, v);
    
    delete[] f;
}
