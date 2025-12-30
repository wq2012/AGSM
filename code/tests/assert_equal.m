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

function assert_equal(actual, expected, tol, name)
    if nargin < 3 || isempty(tol)
        tol = 1e-6;
    end
    if nargin < 4
        name = 'unnamed test';
    end
    
    if all(size(actual) == size(expected)) && all(abs(actual(:) - expected(:)) < tol)
        fprintf('PASS: %s\n', name);
    else
        fprintf('FAIL: %s\n', name);
        fprintf('  Actual: %s\n', mat2str(actual));
        fprintf('  Expected: %s\n', mat2str(expected));
        error('Assertion failed for %s', name);
    end
end
