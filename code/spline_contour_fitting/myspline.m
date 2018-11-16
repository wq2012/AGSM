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

function DD=myspline(D,theta)
%%  Get the distances to center of a cubic spline contour
%   D: distances from landmarks to center
%   theta: angles of all points on contour
%   DD: resulting distances
NN=length(D);
lm_theta=linspace(0,2*pi,NN+1);
lm_theta=[lm_theta(end-3:end-1)-2*pi, lm_theta, lm_theta(2:3)+2*pi];
DD=spline(lm_theta,[D(end-2:end), D, D(1:3)],theta);