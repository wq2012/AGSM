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

function test_circle()
    fprintf('\nTesting Circle Fitting Functions...\n');
    addpath('..');
    addpath('../force_field');
    addpath('../circle_fitting');
    
    % Data generation
    rows = 100; cols = 100;
    xc0 = 50; yc0 = 50; r0 = 30;
    I = zeros(rows, cols);
    theta = 0:0.1:2*pi;
    xi = round(xc0 + r0*cos(theta));
    yi = round(yc0 + r0*sin(theta));
    valid = xi >= 1 & xi <= cols & yi >= 1 & yi <= rows;
    for i=find(valid)
        I(yi(i), xi(i)) = 100;
    end
    [u, v] = GVF(gaussianBlur(I, 2), 1, 0.1, 10);

    % test circleHough
    [xc_h, yc_h, r_h] = circleHough(I, 20, 40, 1, 5);
    assert_equal(abs(xc_h - xc0) < 5, true, [], 'circleHough xc accuracy');
    assert_equal(abs(r_h - r0) < 5, true, [], 'circleHough r accuracy');

    % test CircleHoughTransform (MEX directly)
    H = CircleHoughTransform(double(I), 20, 40, 1);
    assert_equal(size(H), [rows, cols, 21], [], 'CircleHoughTransform size');

    % test fit_circle_force
    init = [xc0 + 2, yc0 - 2, r0 + 1];
    inc = [0.1, 0.1, 0.1];
    thr = [1e-6, 1e-6, 1e-6];
    [xc_f, yc_f, r_f] = fit_circle_force(init, inc, thr, u, v, 50);
    assert_equal(abs(xc_f - xc0) < 5, true, [], 'fit_circle_force xc accuracy');

    % test CircleNormalForce
    [xi_c, yi_c] = circle_in_image(rows, cols, xc0, yc0, r0);
    th_c = atan2(yi_c - yc0, xi_c - xc0);
    cnf = CircleNormalForce(u, v, xi_c, yi_c, th_c);
    assert_equal(isscalar(cnf), true, [], 'CircleNormalForce output type');

    % test circle_in_image
    [xi, yi] = circle_in_image(rows, cols, xc0, yc0, r0);
    assert_equal(length(xi) > 0, true, [], 'circle_in_image non-empty');

    fprintf('Circle fitting tests completed successfully.\n');
end
