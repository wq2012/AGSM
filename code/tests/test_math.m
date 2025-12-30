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

function test_math()
    fprintf('\nTesting Math and Utility Functions...\n');
    addpath('..');
    addpath('../math');
    addpath('../force_field');
    addpath('../spline_contour_fitting');
    addpath('../circle_fitting');
    
    % test is_octave
    assert_equal(is_octave(), exist('OCTAVE_VERSION', 'builtin') > 0, [], 'is_octave()');
    
    % test gaussianBlur
    I = zeros(10, 10);
    I(5, 5) = 1;
    I2 = gaussianBlur(I, 1);
    assert_equal(size(I2), [10, 10], [], 'gaussianBlur size');
    assert_equal(sum(I2(:)) > 0, true, [], 'gaussianBlur non-zero');
    
    % test correctCurve
    r = 50;
    sig = 5;
    rc = correctCurve(r, sig, 100);
    assert_equal(rc > 0, true, [], 'correctCurve positive');
    
    % test correctCurve_polar
    rcp = correctCurve_polar(r, 0, 1/r, sig, 100);
    assert_equal(abs(rcp - rc) < 1e-4, true, [], 'correctCurve_polar consistency');
    
    % test myspline
    D0 = [50 60 70 80 90 100];
    theta = [0 pi/2 pi];
    DD = myspline(D0, theta);
    assert_equal(length(DD), 3, [], 'myspline output length');
    assert_equal(DD(1), D0(1), 1e-4, 'myspline first point');

    % test gaussian3
    h = gaussian3(5);
    assert_equal(size(h), [5 5 5], [], 'gaussian3 size');
    assert_equal(abs(sum(h(:)) - 1.0) < 1e-6, true, [], 'gaussian3 sum for normalization');
    
    fprintf('Math tests completed successfully.\n');
end
