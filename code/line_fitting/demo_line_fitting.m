% This is a simple example showing how to fit a straight line
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
addpath('../force_field');
rng('shuffle');

%% 1. experiment set up
% image size = 400*500
rows=400;
cols=500;

% number of data points and outliers
num_data=50;
num_outlier=5;
noise=5;

% ground truth parameters
theta0=rand(1)*2*pi;
s0=cols/2*cos(theta0)+rows/2*sin(theta0)+20*(rand(1)-0.5);

% image blur parameters
sigma=50;

%% 2. generate data
% ranges
x_min=max(1, min(  ceil((s0-sin(theta0))/cos(theta0)), ceil((s0-rows*sin(theta0))/cos(theta0) ) ) );
x_max=min(cols, max(  floor((s0-sin(theta0))/cos(theta0)), floor((s0-rows*sin(theta0))/cos(theta0) ) ) );
y_min=max(1, min(  ceil((s0-cos(theta0))/sin(theta0)), ceil((s0-cols*cos(theta0))/sin(theta0) ) ) );
y_max=min(rows, max(  floor((s0-cos(theta0))/sin(theta0)), floor((s0-cols*cos(theta0))/sin(theta0) ) ) );

% generate data points 
if mod(theta0,pi)>pi/4 && mod(theta0,pi)<pi*3/4 % flat line
    x=rand(num_data,1);
    x=round(x*(x_max-x_min))+x_min;
    y=s0-x*cos(theta0);y=y/sin(theta0);
else % steep line
    y=rand(num_data,1);
    y=round(y*(y_max-y_min))+y_min;
    x=s0-y*sin(theta0);x=x/cos(theta0);
end
x=x+noise*randn(size(x));
y=y+noise*randn(size(y));

% generate outliers
x=[x;rand(num_outlier,1)*cols];
y=[y;rand(num_outlier,1)*rows];
x=round(x);
y=round(y);
x(x<=1)=2;
x(x>=cols)=cols-1;
y(y<=1)=2;
y(y>=rows)=rows-1;

% display data points and outliers
plot(x,y,'.');
axis equal;
axis([1 cols 1 rows]);
title('data points and outliers');
xlabel('x');
ylabel('y');
drawnow;

% generate image
I=zeros(rows,cols);
for i=1:length(y)
    I(y(i),x(i))=100;
end
I=double(I);

%% 3. GVF field
I2=gaussianBlur(I,sigma);
[u,v] = GVF(I2, 1 , 0.1, 50);
dx=u;dy=v;

%% 4. line fitting
guess=rand(1)*2*pi; % random guess of theta
init=[guess,cols/2*cos(guess)+rows/2*sin(guess)]; % initial parameters
increment=[pi/180*0.2,0.2]; % increment in each iteration
threshold=[0.000001,0.000001]; % threshold for forces/torques

disp('Fitting the line...');
tic;
[theta, s]=fit_line_force(init,increment,threshold,dx,dy,500);
t=toc;
fprintf('Running AGSM takes %f seconds \n',t);

if s*s0<0
    s=-s;
    theta=theta+pi;
end
theta=mod(theta,2*pi);

% display fitting results
figure;
plot(x,y,'.');
axis equal;
axis([1 cols 1 rows]);
hold on;
[xx,yy]=LineInImage(rows,cols,theta,s);
xx=[xx(1) xx(end)];
yy=[yy(1) yy(end)];
plot(xx,yy,'-.r','LineWidth',2);
plot(x,y,'.');
legend('noisy data points and outliers','active geometric shape model fit');
title('fitting a line using active geometric shape model');
xlabel('x');
ylabel('y');
drawnow;

%% 5. error analysis
% total least squares analysis
[theta_TLS, s_TLS, error_TLS]=totalLeastSquares(x,y);
if s_TLS*s0<0
    s_TLS=-s_TLS;
    theta_TLS=theta_TLS+pi;
end
theta_TLS=mod(theta_TLS,2*pi);

% ransac total least squares
[theta_RTLS,s_RTLS]=ransac_TLS(x,y);
if s_RTLS*s0<0
    s_RTLS=-s_RTLS;
    theta_RTLS=theta_RTLS+pi;
end
theta_RTLS=mod(theta_RTLS,2*pi);

% Hough transform
[theta_Hough,s_Hough]=hough_line(I);

if s_Hough*s0<0
    s_Hough=-s_Hough;
    theta_Hough=theta_Hough+pi;
end
theta_Hough=mod(theta_Hough,2*pi);

% get errors
error_theta_AGSM=abs(theta-theta0);
error_s_AGSM=abs(s-s0);
error_theta_TLS=abs(theta_TLS-theta0);
error_s_TLS=abs(s_TLS-s0);
error_theta_RTLS=abs(theta_RTLS-theta0);
error_s_RTLS=abs(s_RTLS-s0);
error_theta_Hough=abs(theta_Hough-theta0);
error_s_Hough=abs(s_Hough-s0);

% display errors
fprintf('num_data=%d    num_outlier=%d \n',num_data,num_outlier);
fprintf('------------------------------------------------------\n');
fprintf('                   theta       s \n');
fprintf('True:              %.4f      %.2f  \n',theta0,s0);
fprintf('TLS fit:           %.4f      %.2f  \n',theta_TLS,s_TLS);
fprintf('RANSAC TLS fit:    %.4f      %.2f  \n',theta_RTLS,s_RTLS);
fprintf('HT fit:            %.4f      %.2f  \n',theta_Hough,s_Hough);
fprintf('AGSM fit:          %.4f      %.2f  \n',theta,s);
fprintf('TLS error:         %.4f      %.2f  \n',error_theta_TLS,error_s_TLS);
fprintf('RANSAC TLS error:  %.4f      %.2f  \n',error_theta_RTLS,error_s_RTLS);
fprintf('HT error:          %.4f      %.2f  \n',error_theta_Hough,error_s_Hough);
fprintf('AGSM error:        %.4f      %.2f  \n',error_theta_AGSM,error_s_AGSM);


