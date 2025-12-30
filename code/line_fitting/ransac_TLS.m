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

function [theta, s, inliers] = ransac_TLS(x, y, n, t, d)
% RANSAC_TLS Robust line fitting using RANSAC and Total Least Squares.
%    [theta, s, inliers] = ransac_TLS(x, y, n, t, d) fits a line to
%    points (x, y) while being robust to outliers.
%
%    Inputs:
%        x, y - Coordinates of points.
%        n - Number of iterations (default: 100).
%        t - Distance threshold for inliers (default: 2.0).
%        d - Minimum number of inliers required (default: round(0.5*length(x))).
%
%    Outputs:
%        theta - Angle of the fitted line.
%        s - Offset of the fitted line.
%        inliers - Logical array of inlier indices.

    if nargin < 3, n = 100; end
    if nargin < 4, t = 2.0; end
    if nargin < 5, d = round(0.5 * length(x)); end

    inliers = [];
    theta = 0;
    s = 0;
    max_inliers = 0;
    N = length(x);
    if N < 2, return; end

    for i = 1:n
        % Randomly select 2 points
        idx = randperm(N, 2);
        x_sample = x(idx);
        y_sample = y(idx);
        
        % Fit line to these 2 points
        [theta_tmp, s_tmp] = totalLeastSquares(x_sample, y_sample);
        
        % Check inliers
        dist = abs(x * cos(theta_tmp) + y * sin(theta_tmp) - s_tmp);
        inliers_tmp = dist < t;
        num_inliers = sum(inliers_tmp);
        
        if num_inliers > max_inliers && num_inliers >= d
            max_inliers = num_inliers;
            inliers = inliers_tmp;
            theta = theta_tmp;
            s = s_tmp;
        end
    end
    
    % Re-fit using all inliers
    if ~isempty(inliers)
        [theta, s] = totalLeastSquares(x(inliers), y(inliers));
    end
end
