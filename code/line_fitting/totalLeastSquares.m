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

function [theta, s, mse] = totalLeastSquares(x, y)
% TOTALLEASTSQUARES Fit a line to points using Total Least Squares.
%    [theta, s, mse] = totalLeastSquares(x, y) calculates the line
%    parameters that minimize the perpendicular distances.
%
%    Inputs:
%        x, y - Coordinates of points.
%
%    Outputs:
%        theta - Angle of the fitted line.
%        s - Offset of the fitted line.
%        mse - Mean squared error.

    x_bar = mean(x);
    y_bar = mean(y);
    x2_bar = mean(x.^2);
    y2_bar = mean(y.^2);
    xy_bar = mean(x.*y);
    
    X = [x2_bar - x_bar^2, xy_bar - x_bar*y_bar; ...
         xy_bar - x_bar*y_bar, y2_bar - y_bar^2];
    
    [V, D] = eig(X);
    
    % The eigenvector corresponding to the smallest eigenvalue gives the normal direction
    [~, idx] = min(diag(D));
    a = V(1, idx);
    b = V(2, idx);
    
    theta = atan2(b, a);
    s = x_bar * cos(theta) + y_bar * sin(theta);
    
    % Mean squared error
    mse = mean((x * cos(theta) + y * sin(theta) - s).^2);
end
