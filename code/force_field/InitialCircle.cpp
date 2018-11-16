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
 * This is the C++/MEX code for computing initial guesses of circle parameters
 *
 * compile:
 *     mex InitialCircle.cpp
 *
 * usage:
 *     [xc,yc,r] = InitialCircle(f)
 *       f: 2D energy field (or gradient field)
 *       xc: initial guess of x-coordinate of center
 *       yc: initial guess of y-coordinate of center
 *       r: initial guess of radius
 */

#include "mex.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

void guess(double *I, int m, int n, double *xc_, double *yc_, double *r_)
{
    double xc=0;
    double yc=0;
    double r=0;
    double sumI=0;
    
    for(int i=0;i<m;i++)
    {
        for(int j=0;j<n;j++)
        {
            xc+=(j*I[i+j*m]);
            yc+=(i*I[i+j*m]);
            sumI+=I[i+j*m];
        }
    }
    xc/=sumI;
    yc/=sumI;
    
    for(int i=0;i<m;i++)
    {
        for(int j=0;j<n;j++)
        {
            r+=(sqrt((i-yc)*(i-yc)+(j-xc)*(j-xc))*I[i+j*m]);
        }
    }
    r/=sumI;
    
    *xc_=xc;
    *yc_=yc;
    *r_=r;
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    double *f; // input matrix
    int m;
    int n;
    
    double *xc;
    double *yc;
    double *r;

    /*  check for proper number of arguments */
    if(nrhs!=1)
    {
        mexErrMsgIdAndTxt( "MATLAB:InitialCircle:invalidNumInputs",
                "One input required.");
    }
    if(nlhs!=3)
    {
        mexErrMsgIdAndTxt( "MATLAB:InitialCircle:invalidNumOutputs",
                "Three outputs required.");
    }
    
    /*  get f */
    f=mxGetPr(prhs[0]);
    m=(int)mxGetM(prhs[0]);
    n=(int)mxGetN(prhs[0]);
    if(m<5 || n<5)
    {
        mexErrMsgIdAndTxt( "MATLAB:InitialCircle:fWrongSize",
                "Input f is too small.");
    }

    /*  set the output pointers to the output matrices */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    
    /*  create C++ pointers to copies of the output matrices */
    xc=mxGetPr(plhs[0]);
    yc=mxGetPr(plhs[1]);
    r=mxGetPr(plhs[2]);
    
    /*  call the C++ subroutine */
    guess(f,m,n,xc,yc,r);
    
    return;
}

