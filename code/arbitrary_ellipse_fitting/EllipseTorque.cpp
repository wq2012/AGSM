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
 * torque=EllipseTorque(u,v,x,y,theta,xc,yc,phi)
 */

#include "mex.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

using namespace std;

/* the main part */
double getTorque(double *u, double *v, double *x, double *y, double *theta, int m, int n, int N, double xc, double yc, double phi)
{
    double torque;
    double d;
    int r,c;
    double fx,fy;

    torque=0;
    for(int i=0;i<N;i++)
    {
        r=(int)y[i]-1;
        c=(int)x[i]-1;
        fx=u[r+c*m];
        fy=v[r+c*m];
        d=sqrt( (x[i]-xc)*(x[i]-xc) + (y[i]-yc)*(y[i]-yc) );
        
        torque += ( (-fx*sin(theta[i]+phi)+fy*cos(theta[i]+phi))*d );
    }
    torque/=(N*(double)N);
    
    return torque;
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    double *u;
    double *v;
    double *x;
    double *y;
    double *theta;
    double xc, yc, phi;
    int m; // number of rows
    int n; // number of columns
    int N; // number of points on line
    double torque;
    
    /*  check for proper number of arguments */
    if(nrhs!=8)
    {
        mexErrMsgIdAndTxt( "MATLAB:EllipseTorque:invalidNumInputs",
                "Eight inputs required.");
    }
    if(nlhs>1)
    {
        mexErrMsgIdAndTxt( "MATLAB:EllipseTorque:invalidNumOutputs",
                "At most one output.");
    }
    
    /*  get u and v  */
    u=mxGetPr(prhs[0]);
    v=mxGetPr(prhs[1]);
    m=(int)mxGetM(prhs[0]);
    n=(int)mxGetN(prhs[0]);
    if(m!=mxGetM(prhs[1]) || n!=mxGetN(prhs[1]))
    {
        mexErrMsgIdAndTxt( "MATLAB:EllipseTorque:invalidInputDimension",
                "Input u and v should have same dimensions.");
    }
    
    /*  get x, y and theta  */
    x=mxGetPr(prhs[2]);
    y=mxGetPr(prhs[3]);
    theta=mxGetPr(prhs[4]);
    N=mxGetN(prhs[2]);
    if(mxGetM(prhs[2])!=1 || mxGetM(prhs[3])!=1 || mxGetM(prhs[4])!=1 || mxGetN(prhs[3])!=N || mxGetN(prhs[4])!=N)
    {
        mexErrMsgIdAndTxt( "MATLAB:EllipseTorque:invalidInputDimension",
                "Input x, y and theta should have same dimensions 1*N.");
    }
    
    /*  get xc  */
    if( !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) ||
            mxGetN(prhs[5])*mxGetM(prhs[5])!=1 )
    {
        mexErrMsgIdAndTxt( "MATLAB:EllipseTorque:xcNotScalar",
                "Input xc must be a scalar.");
    }
    xc=mxGetScalar(prhs[5]);
    
    /*  get phi  */
    if( !mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) ||
            mxGetN(prhs[6])*mxGetM(prhs[6])!=1 )
    {
        mexErrMsgIdAndTxt( "MATLAB:EllipseTorque:ycNotScalar",
                "Input yc must be a scalar.");
    }
    yc=mxGetScalar(prhs[6]);
    
    /*  get phi  */
    if( !mxIsDouble(prhs[7]) || mxIsComplex(prhs[7]) ||
            mxGetN(prhs[7])*mxGetM(prhs[7])!=1 )
    {
        mexErrMsgIdAndTxt( "MATLAB:EllipseTorque:phiNotScalar",
                "Input phi must be a scalar.");
    }
    phi=mxGetScalar(prhs[7]);
    
    /*  call the C++ subroutine  */
    torque=getTorque(u,v,x,y,theta,m,n,N,xc,yc,phi);
    
    /*  output  */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *result=mxGetPr(plhs[0]);
    result[0]=torque;
}