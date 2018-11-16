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
 * This is the C++/MEX code for computing Gradient Vector Flow (GVF)
 *
 * compile:
 *     mex GVF.cpp
 *
 * usage:
 *     [u,v] = GVF(f, alpha, mu, iter)
 *       f: 2D energy field (or gradient field)
 *       alpha: step size of update, between 0 and 1
 *       mu: GVF regularization coefficient
 *       iter: number of iterations
 *       u: resulting GVF in x direction
 *       v: resulting GVF in y direction
 */

#include "mex.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

double *boundMirrorExpand(double *A, int m, int n)
{
    double *B=new double[(m+2)*(n+2)];
    B[0+0*(m+2)]=A[1+1*m];
    B[(m+1)+0*(m+2)]=A[(m-2)+1*m];
    B[0+(n+1)*(m+2)]=A[1+(n-2)*m];
    B[(m+1)+(n+1)*(m+2)]=A[(m-2)+(n-2)*m];
    for(int i=1;i<n+1;i++)
    {
        B[0+i*(m+2)]=A[1+(i-1)*m];
        B[m+1+i*(m+2)]=A[m-2+(i-1)*m];
    }
    for(int i=1;i<m+1;i++)
    {
        B[i+0*(m+2)]=A[(i-1)+1*m];
        B[i+(n+1)*(m+2)]=A[i-1+(n-2)*m];
    }
    for(int i=1;i<m+1;i++)
    {
        for(int j=1;j<n+1;j++)
        {
            B[i+j*(m+2)]=A[i-1+(j-1)*m];
        }
    }
    return B;
}

double *boundMirrorShrink(double *A, int m, int n)
{
    double *B=new double[(m-2)*(n-2)];
    for(int i=0;i<m-2;i++)
    {
        for(int j=0;j<n-2;j++)
        {
            B[i+j*(m-2)]=A[i+1+(j+1)*m];
        }
    }
    return B;
}

void boundMirrorEnsure(double *A, int m, int n)
{
    A[0+0*m]=A[2+2*m];
    A[m-1+0*m]=A[m-3+2*m];
    A[0+(n-1)*m]=A[2+(n-3)*m];
    A[m-1+(n-1)*(m-1)]=A[(m-3)+(n-3)*m];
    for(int i=1;i<n-1;i++)
    {
        A[0+i*m]=A[2+i*m];
        A[m-1+i*m]=A[m-3+i*m];
    }
    for(int i=1;i<m-1;i++)
    {
        A[i+0*m]=A[i+2*m];
        A[i+(n-1)*m]=A[i+(n-3)*m];
    }
}

void copyMatrix(double *A, double *B, int m, int n)
{
    for(int i=0;i<m*n;i++)
    {
        B[i]=A[i];
    }
}

void gradient(double *f, double *fx, double *fy, int m, int n)
{
    for(int i=0;i<m;i++)
    {
        for(int j=0;j<n;j++)
        {
            // fx
            if(j==0)
            {
                fx[i+j*m]=f[i+(j+1)*m]-f[i+j*m];
            }
            else if(j==n-1)
            {
                fx[i+j*m]=f[i+j*m]-f[i+(j-1)*m];
            }
            else
            {
                fx[i+j*m]=(f[i+(j+1)*m]-f[i+(j-1)*m])/2;
            }
            
            // fy
            if(i==0)
            {
                fy[i+j*m]=f[i+1+j*m]-f[i+j*m];
            }
            else if(i==m-1)
            {
                fy[i+j*m]=f[i+j*m]-f[i-1+j*m];
            }
            else
            {
                fy[i+j*m]=(f[i+1+j*m]-f[i-1+j*m])/2;
            }
        }
    }
}

void del2(double *f, double *del2f, int m, int n)
{
    for(int i=0;i<m;i++)
    {
        for(int j=0;j<n;j++)
        {
            if(i==0 || i==m-1 || j==0 || j==n-1)
            {
                del2f[i+j*m]=0;
            }
            else
            {
                del2f[i+j*m]=f[i+1+j*m]+f[i-1+j*m]+f[i+(j+1)*m]+f[i+(j-1)*m]-4*f[i+j*m];
            }
        }
    }
}

void gvf(double *f, int m, int n, double alpha, double mu, int iter, double *u, double *v)
{
    // normalization
    double fmin,fmax;
    fmin=f[0];
    fmax=f[0];
    for(int i=1;i<m*n;i++)
    {
        if(f[i]>fmax) fmax=f[i];
        if(f[i]<fmin) fmin=f[i];
    }
    if(fmax<=fmin)
    {
        mexErrMsgIdAndTxt( "MATLAB:GVF:fWrongRange",
                "Input f is constant.");
    }
    for(int i=0;i<m*n;i++)
    {
        f[i]=(f[i]-fmin)/(fmax-fmin);
    }

    // gradient
    double *f2=boundMirrorExpand(f,m,n);
    double *fx=new double[(m+2)*(n+2)];
    double *fy=new double[(m+2)*(n+2)];
    gradient(f2,fx,fy,m+2,n+2);
    
    // initialization
    double *u1=new double[(m+2)*(n+2)];
    double *v1=new double[(m+2)*(n+2)];
    
    copyMatrix(fx,u1,m+2,n+2);
    copyMatrix(fy,v1,m+2,n+2);
    
    
    
    // square of magnitude
    double *SqrMagf=new double[(m+2)*(n+2)];
    for(int i=0;i<(m+2)*(n+2);i++)
    {
        SqrMagf[i]=fx[i]*fx[i]+fy[i]*fy[i];
    }
    
    double *del2u1=new double[(m+2)*(n+2)];
    double *del2v1=new double[(m+2)*(n+2)];
    
    for(int it=0;it<iter;it++)
    {
        boundMirrorEnsure(u1,m+2,n+2);
        boundMirrorEnsure(v1,m+2,n+2);
        del2(u1,del2u1,m+2,n+2);
        del2(v1,del2v1,m+2,n+2);
        
        for(int i=0;i<(m+2)*(n+2);i++)
        {
            u1[i] = u1[i] + alpha * ( mu*del2u1[i] - SqrMagf[i]*(u1[i]-fx[i]) );
            v1[i] = v1[i] + alpha * ( mu*del2v1[i] - SqrMagf[i]*(v1[i]-fy[i]) );
        }
    }
    
    double *u2=boundMirrorShrink(u1, m+2, n+2);
    double *v2=boundMirrorShrink(v1, m+2, n+2);
    
    copyMatrix(u2,u,m,n);
    copyMatrix(v2,v,m,n);
   
    delete[] f2;  
    delete[] fx;
    delete[] fy;
    delete[] u1;
    delete[] v1;
    delete[] del2u1;
    delete[] del2v1;
    delete[] u2;
    delete[] v2;
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    double *f0; // input matrix
    double *f; // a copy of input matrix
    int m;
    int n;
    
    double alpha;
    double mu;
    
    int iter;
    
    double *u;
    double *v;
    
    
    /*  check for proper number of arguments */
    if(nrhs!=4)
    {
        mexErrMsgIdAndTxt( "MATLAB:GVF:invalidNumInputs",
                "Four inputs required.");
    }
    if(nlhs!=2)
    {
        mexErrMsgIdAndTxt( "MATLAB:GVF:invalidNumOutputs",
                "Two outputs required.");
    }
    
    /*  get f */
    f0=mxGetPr(prhs[0]);
    m=mxGetM(prhs[0]);
    n=mxGetN(prhs[0]);
    if(m<10 || n<10)
    {
        mexErrMsgIdAndTxt( "MATLAB:GVF:fWrongSize",
                "Input f is too small.");
    }
    f=new double[m*n];
    for(int i=0;i<m*n;i++)
    {
        f[i]=f0[i];
    }
    
    /*  get alpha */
    if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
            mxGetN(prhs[1])*mxGetM(prhs[1])!=1 )
    {
        mexErrMsgIdAndTxt( "MATLAB:GVF:alphaNotScalar",
                "Input alpha must be a scalar.");
    }
    alpha=mxGetScalar(prhs[1]);
    if(alpha<0)
    {
        mexErrMsgIdAndTxt( "MATLAB:GVF:alphaWrongRange",
                "Input alpha must be non-negative.");
    }
    
    /*  get mu */
    if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
            mxGetN(prhs[2])*mxGetM(prhs[2])!=1 )
    {
        mexErrMsgIdAndTxt( "MATLAB:GVF:muNotScalar",
                "Input mu must be a scalar.");
    }
    mu=mxGetScalar(prhs[2]);
    if(mu<0)
    {
        mexErrMsgIdAndTxt( "MATLAB:GVF:muWrongRange",
                "Input mu must be non-negative.");
    }
    
    /*  get iter */
    if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
            mxGetN(prhs[3])*mxGetM(prhs[3])!=1 )
    {
        mexErrMsgIdAndTxt( "MATLAB:GVF:iterNotScalar",
                "Input iter must be a scalar.");
    }
    iter=(int)mxGetScalar(prhs[3]);
    if(iter<1)
    {
        mexErrMsgIdAndTxt( "MATLAB:GVF:iterWrongRange",
                "Input iter must be positive.");
    }
    
    
    /*  set the output pointers to the output matrices */
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(m, n, mxREAL);
    
    /*  create C++ pointers to copies of the output matrices */
    u=mxGetPr(plhs[0]);
    v=mxGetPr(plhs[1]);
    
    /*  call the C++ subroutine */
    gvf(f,m,n,alpha,mu,iter,u,v);
    
    delete[] f;
    return;
}

