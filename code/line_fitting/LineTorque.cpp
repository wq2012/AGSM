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
 * [torque,k]=LineTorque(u,v,x,y,theta)
 */

#include "mex.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

using namespace std;

/* the main part */
double getTorque(double *u, double *v, double *x, double *y, int m, int n, int N, double theta, int &k)
{
    double maxTorque=0;
    double torque;
    int bestK;
    double d;
    int sgn;
    int r,c;
    double fx,fy;
    
    for(k=1;k<=N;k++)
    {
        torque=0;
        for(int i=0;i<N;i++)
        {
            r=(int)y[i]-1;
            c=(int)x[i]-1;
            fx=u[r+c*m];
            fy=v[r+c*m];
            d=sqrt( (x[i]-x[k-1])*(x[i]-x[k-1]) + (y[i]-y[k-1])*(y[i]-y[k-1]) );
    
            if(k-1==i)
            {
                sgn=0;
            }
            else if(x[k-1]<x[(i>k-1?N-1:0)])
            {
                sgn=-1;
            }
            else
            {
                sgn=1;
            }

            torque += ( (fx*cos(theta)+fy*sin(theta))*d*sgn );
        }
        torque/=(N*(double)N);
        if(fabs(torque)>fabs(maxTorque))
        {
            maxTorque=torque;
            bestK=k;
        }
    }
    
    k=bestK;
    return maxTorque;
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    double *u;
    double *v;
    double *x;
    double *y;
    double theta;
    int m; // number of rows
    int n; // number of columns
    int N; // number of points on line
    double torque;
    int k;
    
    /*  check for proper number of arguments */
    if(nrhs!=5)
    {
        mexErrMsgIdAndTxt( "MATLAB:LineNormalForce:invalidNumInputs",
                "Five inputs required.");
    }
    if(nlhs>2)
    {
        mexErrMsgIdAndTxt( "MATLAB:LineNormalForce:invalidNumOutputs",
                "At most two outputs.");
    }
    
    /*  get u and v  */
    u=mxGetPr(prhs[0]);
    v=mxGetPr(prhs[1]);
    m=(int)mxGetM(prhs[0]);
    n=(int)mxGetN(prhs[0]);
    if(m!=mxGetM(prhs[1]) || n!=mxGetN(prhs[1]))
    {
        mexErrMsgIdAndTxt( "MATLAB:LineNormalForce:invalidInputDimension",
                "Input u and v should have same dimensions.");
    }
    
    /*  get x and y  */
    x=mxGetPr(prhs[2]);
    y=mxGetPr(prhs[3]);
    N=mxGetN(prhs[2]);
    if(mxGetM(prhs[2])!=1 || mxGetM(prhs[3])!=1 || mxGetN(prhs[3])!=N)
    {
        mexErrMsgIdAndTxt( "MATLAB:LineNormalForce:invalidInputDimension",
                "Input x and y should have same dimensions 1*N.");
    }
    
    /*  get theta  */
    if( !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) ||
            mxGetN(prhs[4])*mxGetM(prhs[4])!=1 )
    {
        mexErrMsgIdAndTxt( "MATLAB:LineInImage:thetaNotScalar",
                "Input theta must be a scalar.");
    }
    theta=mxGetScalar(prhs[4]);
    
    /*  call the C++ subroutine  */
    torque=getTorque(u,v,x,y,m,n,N,theta,k);
    
    /*  output  */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *result1=mxGetPr(plhs[0]);
    double *result2=mxGetPr(plhs[1]);
    result1[0]=torque;
    result2[0]=(double)k;
}