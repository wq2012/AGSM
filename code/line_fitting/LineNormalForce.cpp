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
 * force=LineNormalForce(u,v,x,y,theta)
 */

#include "mex.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

using namespace std;

/* the main part */
double getNormalForce(double *u, double *v, double *x, double *y, int m, int n, int N, double theta)
{
    double force=0;
    int r,c;
    double fx,fy;
    for(int i=0;i<N;i++)
    {
        r=(int)y[i]-1;
        c=(int)x[i]-1;
        fx=u[r+c*m];
        fy=v[r+c*m];
        force+=(fx*cos(theta)+fy*sin(theta));
    }
    return force/N;
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
    double force;
    
    /*  check for proper number of arguments */
    if(nrhs!=5)
    {
        mexErrMsgIdAndTxt( "MATLAB:LineNormalForce:invalidNumInputs",
                "Five inputs required.");
    }
    if(nlhs>1)
    {
        mexErrMsgIdAndTxt( "MATLAB:LineNormalForce:invalidNumOutputs",
                "At most one output.");
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
    force=getNormalForce(u,v,x,y,m,n,N,theta);
    
    /*  output  */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *result=mxGetPr(plhs[0]);
    result[0]=force;
}