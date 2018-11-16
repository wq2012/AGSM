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
 * [x,y]=LineInImage(m,n,theta,s)
 */

#include "mex.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

#define PI 3.14159265

using namespace std;

/* rounding is only supported in C++11 */
double round(double a) 
{
    int b=(int)a;
    a-=b;
    if(a>=0.5) b++;
    return (double)b;
}

/* cpoy a matrix */
void copyMatrix(double *a, double *b, int N)
{
    for(int i=0;i<N;i++)
    {
        b[i]=a[i];
    }
}

/* the main part */
int getLineCoordinates(int m, int n, double theta, double s, double *&x, double *&y)
{
    double theta_mod=theta;
    while(theta_mod>PI) theta_mod-=PI;
    while(theta_mod<0) theta_mod+=PI;
    int N=0; // number of points

    if(theta_mod>PI/4 && theta_mod<PI/4*3) // flat line
    {
        int x1;
        int x2;
        double temp;
        for(x1=1;x1<=n;x1++)
        {
            temp=round((s-x1*cos(theta))/sin(theta));
            if(temp>=1 && temp<=m) break;
        }
        if(x1>n)
        {
            return N;
        }
        for(x2=n;x2>=1;x2--)
        {
            temp=round((s-x2*cos(theta))/sin(theta));
            if(temp>=1 && temp<=m) break;
        }
        N=x2-x1+1;
        x=new double[N];
        y=new double[N];
        
        for(int i=x1;i<=x2;i++)
        {
            x[i-x1]=i;
            y[i-x1]=round((s-i*cos(theta))/sin(theta));
        }
    }
    else // steep line
    {
        int y1;
        int y2;
        double temp;
        for(y1=1;y1<=m;y1++)
        {
            temp=round((s-y1*sin(theta))/cos(theta));
            if(temp>=1 && temp<=n) break;
        }
        if(y1>m)
        {
            return N;
        }
        for(y2=m;y2>=1;y2--)
        {
            temp=round((s-y2*sin(theta))/cos(theta));
            if(temp>=1 && temp<=n) break;
        }
        N=y2-y1+1;
        x=new double[N];
        y=new double[N];
        for(int i=y1;i<=y2;i++)
        {
            y[i-y1]=i;
            x[i-y1]=round((s-i*sin(theta))/cos(theta));
        }
    }
    return N;
}


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    int m; // number of rows
    int n; // number of columns
    double theta; // line parameter
    double s; // line parameter
    double *x; // x-coordinates of line
    double *y; // y-coordinates of line
    
    /*  check for proper number of arguments */
    if(nrhs!=4)
    {
        mexErrMsgIdAndTxt( "MATLAB:LineInImage:invalidNumInputs",
                "Four inputs required.");
    }
    if(nlhs!=2)
    {
        mexErrMsgIdAndTxt( "MATLAB:LineInImage:invalidNumOutputs",
                "Two outputs required.");
    }
    
    /*  get m */
    if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
            mxGetN(prhs[0])*mxGetM(prhs[0])!=1 )
    {
        mexErrMsgIdAndTxt( "MATLAB:LineInImage:mNotScalar",
                "Input m must be a scalar.");
    }
    m=(int)mxGetScalar(prhs[0]);
    if(m<0)
    {
        mexErrMsgIdAndTxt( "MATLAB:LineInImage:mWrongRange",
                "Input m must be non-negative.");
    }
    
    /*  get n */
    if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
            mxGetN(prhs[1])*mxGetM(prhs[1])!=1 )
    {
        mexErrMsgIdAndTxt( "MATLAB:LineInImage:nNotScalar",
                "Input n must be a scalar.");
    }
    n=(int)mxGetScalar(prhs[1]);
    if(n<0)
    {
        mexErrMsgIdAndTxt( "MATLAB:LineInImage:nWrongRange",
                "Input n must be non-negative.");
    }
    
    /*  get theta */
    if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
            mxGetN(prhs[2])*mxGetM(prhs[2])!=1 )
    {
        mexErrMsgIdAndTxt( "MATLAB:LineInImage:thetaNotScalar",
                "Input theta must be a scalar.");
    }
    theta=mxGetScalar(prhs[2]);
    
    /*  get s  */
    if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
            mxGetN(prhs[3])*mxGetM(prhs[3])!=1 )
    {
        mexErrMsgIdAndTxt( "MATLAB:LineInImage:sNotScalar",
                "Input s must be a scalar.");
    }
    s=mxGetScalar(prhs[3]);
    
    /*  call the C++ subroutine  */
    double *xx;
    double *yy;
    int N=getLineCoordinates(m,n,theta,s,xx,yy);
    
    /*  output  */
    plhs[0] = mxCreateDoubleMatrix(1, N, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, N, mxREAL);
    
    /*  create C++ pointers to a copies of the output matrix  */
    x=mxGetPr(plhs[0]);
    y=mxGetPr(plhs[1]);
    
    copyMatrix(xx,x,N);
    copyMatrix(yy,y,N);

    if(N>0)
    {
        delete[] xx;
        delete[] yy;
    }

}