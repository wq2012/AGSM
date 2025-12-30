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

function [theta, s] = agsm_hough_line(I, P)
% AGSM_HOUGH_LINE Find the best fit line using the line Hough transform.
%    [theta, s] = agsm_hough_line(I, P) detects lines in a binary image
%    using a simple Hough transform approach for comparison.
%
%    Inputs:
%        I - 2D binary image.
%        P - Precision (suggested 1; >2 is slow).
%
%    Outputs:
%        theta - Angle of the detected line.
%        s - Offset of the detected line.

    [H, theta_Hough, s_Hough] = hough(I > 0);
    h = fspecial('gaussian', 10, 3);
    H2 = imfilter(H, h);
    [~, maxindex] = max(H2(:));
    
    theta = theta_Hough(floor((maxindex - 1) / size(H, 1)) + 1) / 180 * pi;
    s = s_Hough(mod(maxindex - 1, size(H, 1)) + 1);
end