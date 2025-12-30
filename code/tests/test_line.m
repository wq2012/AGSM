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

function test_line()
    fprintf('\nTesting Line Fitting Functions...\n');
    addpath('..');
    addpath('../force_field');
    addpath('../line_fitting');
    
    % Data generation
    rows = 100; cols = 100;
    theta0 = pi/4; s0 = 50;
    I = zeros(rows, cols);
    for x=1:cols
        y = round((s0 - x*cos(theta0))/sin(theta0));
        if y >= 1 && y <= rows
            I(y, x) = 100;
        end
    end
    [u, v] = GVF(gaussianBlur(I, 2), 1, 0.1, 10);
    
    % test totalLeastSquares
    [x_pts, y_pts] = find(I);
    [theta_tls, s_tls] = totalLeastSquares(y_pts, x_pts); % Note: flipped x,y for image coords
    assert_equal(abs(sin(theta_tls - theta0)) < 0.1, true, [], 'totalLeastSquares accuracy');

    % test ransac_TLS
    [theta_ransac, s_ransac] = ransac_TLS(y_pts, x_pts);
    assert_equal(abs(sin(theta_ransac - theta0)) < 0.1, true, [], 'ransac_TLS accuracy');

    % test agsm_hough_line
    [theta_h, s_h] = agsm_hough_line(I);
    assert_equal(abs(sin(theta_h - theta0)) < 0.1, true, [], 'agsm_hough_line accuracy');
    
    % test fit_line_force
    init = [theta0 + 0.1, s0 + 2];
    inc = [0.01, 0.1];
    thr = [1e-6, 1e-6];
    [theta_fit, s_fit] = fit_line_force(init, inc, thr, u, v, 50);
    assert_equal(abs(sin(theta_fit - theta0)) < 0.1, true, [], 'fit_line_force accuracy');

    % test LineInImage
    [xi, yi] = LineInImage(rows, cols, theta0, s0);
    assert_equal(length(xi) > 0, true, [], 'LineInImage non-empty');

    % test LineNormalForce
    lnf = LineNormalForce(u, v, xi, yi, theta0);
    assert_equal(isscalar(lnf), true, [], 'LineNormalForce output type');

    % test LineTorque
    ltq = LineTorque(u, v, xi, yi, theta0);
    assert_equal(isscalar(ltq), true, [], 'LineTorque output type');

    fprintf('Line fitting tests completed successfully.\n');
end
