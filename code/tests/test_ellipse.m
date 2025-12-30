function test_ellipse()
    fprintf('\nTesting Ellipse Fitting Functions...\n');
    addpath('..');
    addpath('../force_field');
    addpath('../arbitrary_ellipse_fitting');
    
    % Data generation
    rows = 100; cols = 100;
    xc0 = 50; yc0 = 50; a0 = 40; b0 = 20; phi0 = pi/6;
    I = zeros(rows, cols);
    theta = 0:0.1:2*pi;
    xi = round(xc0 + a0*cos(theta)*cos(phi0) - b0*sin(theta)*sin(phi0));
    yi = round(yc0 + a0*cos(theta)*sin(phi0) + b0*sin(theta)*cos(phi0));
    valid = xi >= 1 & xi <= cols & yi >= 1 & yi <= rows;
    for i=find(valid)
        I(yi(i), xi(i)) = 100;
    end
    [u, v] = GVF(gaussianBlur(I, 2), 1, 0.1, 10);

    % test fit_arb_ellipse_force
    init = [xc0, yc0, a0, b0, phi0];
    inc = [0.1, 0.1, 0.1, 0.1, 0.01];
    thr = ones(1, 5) * 1e-6;
    bound = [10, 100, 10, 100];
    [xc_f, yc_f, a_f, b_f, phi_f] = fit_arb_ellipse_force(init, inc, thr, bound, u, v, 50, 0);
    assert_equal(abs(xc_f - xc0) < 5, true, [], 'fit_arb_ellipse_force xc accuracy');

    % test EllipseTorque
    [xi_e, yi_e] = arb_ellipse_in_image(rows, cols, xc0, yc0, a0, b0, phi0);
    th_e = rand(size(xi_e)); % Dummy theta vector of same size
    etq = EllipseTorque(u, v, xi_e, yi_e, th_e, xc0, yc0, phi0);
    assert_equal(isscalar(etq), true, [], 'EllipseTorque output type');

    % test arb_ellipse_in_image
    [xi, yi] = arb_ellipse_in_image(rows, cols, xc0, yc0, a0, b0, phi0);
    assert_equal(length(xi) > 0, true, [], 'arb_ellipse_in_image non-empty');

    fprintf('Ellipse fitting tests completed successfully.\n');
end
