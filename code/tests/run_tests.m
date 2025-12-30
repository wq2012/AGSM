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

function run_tests()
    fprintf('Starting AGSM Unit Test Suite...\n');
    
    if exist('OCTAVE_VERSION', 'builtin') > 0
        pkg load image;
    end
    
    try
        test_math();
        test_force_field();
        test_line();
        test_circle();
        test_ellipse();
        test_spline();
        
        fprintf('\n============================\n');
        fprintf('ALL TESTS PASSED SUCCESSFULLY\n');
        fprintf('============================\n');
        exit(0);
    catch err
        fprintf('\n============================\n');
        fprintf('TEST SUITE FAILED\n');
        fprintf('Error: %s\n', err.message);
        fprintf('============================\n');
        exit(1);
    end
end
