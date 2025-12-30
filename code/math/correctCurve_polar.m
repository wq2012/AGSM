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

function R = correctCurve_polar(r, r1, r2, sigma, iter)
% CORRECTCURVE_POLAR Correct for the radius of curvature in polar system.
%    R = correctCurve_polar(r, r1, r2, sigma, iter) corrects the detected 
%    radius in a polar coordinate system using derivatives.
%
%    Inputs:
%        r - Detected radius.
%        r1 - 1st order derivative of r w.r.t. theta.
%        r2 - 2nd order derivative of r w.r.t. theta.
%        sigma - Standard deviation of Gaussian PSF.
%        iter - Number of iterations.
%
%    Outputs:
%        R - Corrected radius.

    kappa = (r^2 + 2*r1^2 - r*r2) / (r^2 + r1^2).^1.5;
    Kappa = 1 / correctCurve(abs(1 / kappa), sigma, iter);

    if kappa > 0
        R = r - 1 / abs(kappa) + 1 / Kappa;
    elseif kappa < 0
        R = r + 1 / abs(kappa) - 1 / Kappa;
    else
        R = r;
    end
end
