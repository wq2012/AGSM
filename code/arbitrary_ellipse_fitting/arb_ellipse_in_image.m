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

function [x, y, theta] = arb_ellipse_in_image(m, n, xc, yc, a, b, phi)
% ARB_ELLIPSE_IN_IMAGE Return discrete coordinates of an arbitrary ellipse.
%    [x, y, theta] = arb_ellipse_in_image(m, n, xc, yc, a, b, phi) 
%    generates points on an ellipse that fall within the image boundaries.
%
%    The parametric equations of the ellipse:
%       x = xc + a * cos(theta) * cos(phi) - b * sin(theta) * sin(phi)
%       y = yc + a * cos(theta) * sin(phi) + b * sin(theta) * cos(phi)
%
%    Inputs:
%        m - Number of image rows.
%        n - Number of image columns.
%        xc - X-coordinate of center.
%        yc - Y-coordinate of center.
%        a - Semi-major axis.
%        b - Semi-minor axis.
%        phi - Orientation angle in radians.
%
%    Outputs:
%        x - X-coordinates of points on ellipse.
%        y - Y-coordinates of points on ellipse.
%        theta - Angles corresponding to points.

    theta = linspace(0, 2 * pi, round(2 * pi * max([a, b])));
    data = [xc + a * cos(theta) * cos(phi) - b * sin(theta) * sin(phi); ...
            yc + a * cos(theta) * sin(phi) + b * sin(theta) * cos(phi)];
    data = round(data);
    [data, index] = unique(data', 'rows'); % Remove repeated points
    data = data';
    x = data(1, :);
    y = data(2, :);
    theta = theta(index);
    [theta, index] = sort(theta);
    x = x(index);
    y = y(index);
    
    if min(x) < 1 || min(y) < 1 || max(x) > n || max(y) > m
        warning('AGSM:ArbEllipseInImage:outOfBounds', 'Contour out of image boundaries!');
    end
end