% This is a simple example showing how to fit a circle 
% to data using active geometric shape model.

% Copyright (C) 2012 Quan Wang <wangq10@rpi.edu>,
% Signal Analysis and Machine Perception Laboratory,
% Department of Electrical, Computer, and Systems Engineering,
% Rensselaer Polytechnic Institute, Troy, NY 12180, USA
%
% You are free to use this software for academic purposes if you cite our paper:
% Quan Wang, Kim L. Boyer,
% The active geometric shape model: A new robust deformable shape model and its applications,
% Computer Vision and Image Understanding, Volume 116, Issue 12, December 2012, Pages 1178-1194,
% ISSN 1077-3142, 10.1016/j.cviu.2012.08.004.
%
% For commercial use, please contact the authors.

clear;clc;close all;
addpath('..');
addpath('../force_field');
addpath('../math');
if is_octave()
    rand('state',sum(100*clock));
else
    rng('shuffle');
end

if is_octave()
    % prevent opening figures in headless mode
    set(0, 'defaultfigurevisible', 'off');
end

%% 1. experiment set up
% image size = 400*500
rows=400;
cols=500;

% number of data points and outliers
num_data=50;
num_outlier=10;
noise=5;

% ground truth parameters
x0=250+(rand(1)-0.5)*50;
y0=200+(rand(1)-0.5)*50;
r0=50+rand(1)*50;

% image blur parameters
sigma=20;

%% 2. generate data
% generate data points
theta=rand(num_data,1)*2*pi;
x=round(x0+r0*cos(theta)+noise*randn(size(theta)));
y=round(y0+r0*sin(theta)+noise*randn(size(theta)));

% generate outliers
x=[x;rand(num_outlier,1)*cols];
y=[y;rand(num_outlier,1)*rows];
x=round(x);
y=round(y);
x(x<=1)=2;
x(x>=cols)=cols-1;
y(y<=1)=2;
y(y>=rows)=rows-1;


% display data points 
plot(x,y,'.');
axis equal;
axis([1 cols 1 rows]);
title('data points and outliers');
xlabel('x');
ylabel('y');
drawnow;

%% generate image
I0=zeros(rows,cols);
for i=1:length(y)
    I0(y(i),x(i))=100;
end
I0=double(I0);

%% 3. GVF field
I=gaussianBlur(I0,sigma);
[u,v] = GVF(I, 1 , 0.1, 50);
dx=u;dy=v;

%% 4. circle fitting
[xc, yc, r] = InitialCircle(I); % initial guess
init=[xc,yc,r]; % initial parameters
increment=[0.2,0.2,0.2]; % increment in each iteration
threshold=[0.000001,0.000001,0.000001]; % threshold for forces/torques

disp('Fitting the circle...');
tic;
[xc, yc, r]=fit_circle_force(init,increment,threshold,dx,dy,500);
t=toc;
fprintf('Running AGSM takes %f seconds \n',t);

% correction of curvature
r_before_correction=r;
r=correctCurve(r,sigma,100);

% display fitting results
figure;
plot(x,y,'.');
axis equal;
axis([1 cols 1 rows]);
hold on;
theta=0:0.01:2*pi;
plot(xc+r*cos(theta),yc+r*sin(theta),'-.r','LineWidth',2);
legend('noisy data points and outliers','active geometric shape model fit');
title('fitting a circle using active geometric shape model');
xlabel('x');
ylabel('y');
drawnow;

%% 5. error analysis
% circle Hough transform
rmin=50;
rmax=100;
P=1;
FS=5;
[xc_Hough,yc_Hough,r_Hough]=circleHough(I0,rmin,rmax,P,FS);

% display errors
fprintf('num_data=%d    num_outlier=%d \n',num_data,num_outlier);
fprintf('------------------------------------------------------\n');
fprintf('                     xc          yc          r \n');
fprintf('True:                %.4f    %.4f    %.4f  \n',x0,y0,r0);
fprintf('AGSM fit:            %.4f    %.4f    %.4f  \n',xc,yc,r);
fprintf('Hough transform:     %.4f    %.4f    %.4f  \n',xc_Hough,yc_Hough,r_Hough);
