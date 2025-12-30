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

function [xc, yc, r] = circleHough(I, rmin, rmax, P, FS)
% CIRCLEHOUGH Find the best fit circle using the circle Hough transform.
%    [xc, yc, r] = circleHough(I, rmin, rmax, P, FS) detects circles in a 
%    binary image.
%
%    Inputs:
%        I - 2D binary image.
%        rmin - Smallest radius to be considered.
%        rmax - Largest radius to be considered.
%        P - Precision (suggested 1; >2 is slow).
%        FS - Gaussian filter size (suggested 5; >10 is slow).
%
%    Outputs:
%        xc - X-coordinate of detected center.
%        yc - Y-coordinate of detected center.
%        r - Detected radius.

    % Input checking
    xc = []; yc = []; r = [];
    if length(size(I)) > 2
        error('AGSM:circleHough:invalidInput', 'I must be a 2D binary image.');
    end
    if rmin >= rmax || rmin <= 0 || rmax <= 0
        error('AGSM:circleHough:invalidRange', 'Must ensure 0 < rmin < rmax.');
    end
    if rmin ~= round(rmin) || rmax ~= round(rmax)
        warning('AGSM:circleHough:nonIntegerRange', 'rmin and rmax should be integers.');
    end
    if P < 1 || round(P) ~= P || P > 100
        error('AGSM:circleHough:invalidP', 'P must be a small positive integer.');
    end
    if FS < 1 || round(FS) ~= FS || FS > 100
        error('AGSM:circleHough:invalidFS', 'FS must be a small positive integer.');
    end

    % Hough transform
    I = double(I);
    H = CircleHoughTransform(I, rmin, rmax, P);

    % Gaussian filtering
    h = gaussian3(FS);
    fprintf('Performing 3D Gaussian filtering ... \n');
    if is_octave()
        H = convn(H, h, 'same');
    else
        H = imfilter(H, h);
    end

    % Extract best parameters
    [~, index] = max(H(:));
    y_index = mod(index - 1, size(H, 1)) + 1;
    index = (index - y_index) / size(H, 1) + 1;
    x_index = mod(index - 1, size(H, 2)) + 1;
    r_index = (index - x_index) / size(H, 2) + 1;

    xc = (x_index - 1) / P + 1;
    yc = (y_index - 1) / P + 1;
    r = (r_index - 1) / P + rmin;
end