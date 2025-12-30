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

function test_spline()
    fprintf('\nTesting Spline Fitting Functions...\n');
    addpath('..');
    addpath('../force_field');
    addpath('../spline_contour_fitting');
    
    % Data generation
    rows = 100; cols = 100;
    xc0 = 50; yc0 = 50; D0 = [30 35 30 35 30 35];
    I = zeros(rows, cols);
    [xi, yi] = spline_contour_in_image(rows, cols, xc0, yc0, D0);
    for i=1:length(xi)
        I(yi(i), xi(i)) = 100;
    end
    [u, v] = GVF(gaussianBlur(I, 2), 1, 0.1, 10);

    % test fit_spline_contour_force
    init = [xc0, yc0, D0];
    inc = [0.1, 0.1, 0.1];
    thr = [1e-6, 1e-6, 1e-6];
    bound = [10, 80];
    [xc_f, yc_f, D_f] = fit_spline_contour_force(init, inc, thr, bound, u, v, 50, 0);
    assert_equal(abs(xc_f - xc0) < 5, true, [], 'fit_spline_contour_force xc accuracy');

    % test SplineContourForce
    [xi_s, yi_s, th_s] = spline_contour_in_image(rows, cols, xc0, yc0, D0);
    scf = SplineContourForce(u, v, xi_s, yi_s, th_s);
    assert_equal(isscalar(scf), true, [], 'SplineContourForce output type');

    % test spline_contour_in_image
    [xi, yi] = spline_contour_in_image(rows, cols, xc0, yc0, D0);
    assert_equal(length(xi) > 0, true, [], 'spline_contour_in_image non-empty');

    fprintf('Spline fitting tests completed successfully.\n');
end
