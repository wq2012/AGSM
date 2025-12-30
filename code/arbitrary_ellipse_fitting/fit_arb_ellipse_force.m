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

function [xc, yc, a, b, phi] = fit_arb_ellipse_force(init, increment, threshold, bound, field_x, field_y, iter, show_fitness)
% FIT_ARB_ELLIPSE_FORCE Fit an arbitrary ellipse in the force field.
%    [xc, yc, a, b, phi] = fit_arb_ellipse_force(init, increment, threshold, 
%    bound, field_x, field_y, iter, show_fitness) iteratively fits an 
%    arbitrary ellipse to a vector field.
%
%    Inputs:
%        init - [x0, y0, a0, b0, phi0] initial parameters.
%        increment - [d_xc, d_yc, d_a, d_b, d_phi] step sizes.
%        threshold - [t_xc, t_yc, t_a, t_b, t_phi] update thresholds.
%        bound - [a_low, a_up, b_low, b_up] constraints on axes.
%        field_x - The x component of force field.
%        field_y - The y component of force field.
%        iter - Number of iterations.
%        show_fitness - Flag to toggle fitness function tracking.
%
%    Outputs:
%        xc, yc - Final center coordinates.
%        a, b - Final semi-axes.
%        phi - Final orientation.

    x0 = init(1);
    y0 = init(2);
    a0 = init(3);
    b0 = init(4);
    phi0 = init(5);
    d_xc = increment(1);
    d_yc = increment(2);
    d_a = increment(3);
    d_b = increment(4);
    d_phi = increment(5);
    t_xc = threshold(1);
    t_yc = threshold(2);
    t_a = threshold(3);
    t_b = threshold(4);
    t_phi = threshold(5);
    a_low = bound(1);
    a_up = bound(2);
    b_low = bound(3);
    b_up = bound(4);

    % note: m rows and n columns, but x is column and y is row here
    [m, n] = size(field_x);

    xc = x0;
    yc = y0;
    a = a0;
    b = b0;
    phi = phi0;

    fit_save = zeros(iter, 1);

    for it = 1:iter
        [x, y, theta] = arb_ellipse_in_image(m, n, xc, yc, a, b, phi);
        
        % Calculate torque about center
        torque = EllipseTorque(field_x, field_y, x, y, theta, xc, yc, phi);

        % Update orientation (phi)
        if torque > t_phi
            phi = phi + d_phi;
        elseif torque < -t_phi
            phi = phi - d_phi;
        end
        
        % Calculate average force around contour
        F_round = sum(ContourForceArray(field_x, field_y, x, y));
        F_round = F_round / length(theta);
        
        % Quadrant inward forces
        idx_left = find(theta > pi * 3 / 4 & theta < pi * 5 / 4);
        F_left = ContourDirectionForce(field_x, field_y, x(idx_left), y(idx_left), [cos(phi), sin(phi)]);
        
        idx_right = find(theta < pi / 4 | theta > pi * 7 / 4);
        F_right = ContourDirectionForce(field_x, field_y, x(idx_right), y(idx_right), [-cos(phi), -sin(phi)]);
        
        idx_up = find(theta > pi / 4 & theta < pi * 3 / 4);
        F_up = ContourDirectionForce(field_x, field_y, x(idx_up), y(idx_up), [sin(phi), -cos(phi)]);
        
        idx_down = find(theta > pi * 5 / 4 & theta < pi * 7 / 4);
        F_down = ContourDirectionForce(field_x, field_y, x(idx_down), y(idx_down), [-sin(phi), cos(phi)]);
        
        % Update center (xc, yc)
        F_left_right = dot(F_round, [1, 0]);
        if F_left_right > t_xc
            xc = xc + d_xc;
        elseif F_left_right < -t_xc
            xc = xc - d_xc;
        end
        
        F_down_up = dot(F_round, [0, 1]);
        if F_down_up > t_yc
            yc = yc + d_yc;
        elseif F_down_up < -t_yc
            yc = yc - d_yc;
        end
        
        % Update center again using diagonal forces
        F_diag1 = dot(F_round, [0.7071, 0.7071]);
        if F_diag1 > t_xc + t_yc
            xc = xc + d_xc;
            yc = yc + d_yc;
        elseif F_diag1 < -t_xc - t_yc
            xc = xc - d_xc;
            yc = yc - d_yc;
        end
        
        F_diag2 = dot(F_round, [-0.7071, 0.7071]);
        if F_diag2 > t_xc + t_yc
            xc = xc - d_xc;
            yc = yc + d_yc;
        elseif F_diag2 < -t_xc - t_yc
            xc = xc + d_xc;
            yc = yc - d_yc;
        end
        
        % Update semi-axes (a, b)
        if F_left + F_right > t_a
            a = a - d_a;
        elseif F_left + F_right < -t_a
            a = a + d_a;
        end
        
        if F_up + F_down > t_b
            b = b - d_b;
        elseif F_up + F_down < -t_b
            b = b + d_b;
        end
        
        if b > a
            temp = a; a = b; b = temp;
            phi = mod(phi + pi / 2, pi);
        end
        
        % Constrain semi-axes
        a = max(a_low, min(a, a_up));
        b = max(b_low, min(b, b_up));
        
        % Performance/fitness tracking
        if show_fitness == 1
            beta = 0.9;
            [~, ~, tmp_theta] = arb_ellipse_in_image(m, n, xc, yc, a, b, phi);
            [~, ~, tmp_theta1] = arb_ellipse_in_image(m, n, xc, yc, a * beta, b * beta, phi);
            [~, ~, tmp_theta2] = arb_ellipse_in_image(m, n, xc, yc, a / beta, b / beta, phi);

            fit0 = sum(sqrt(sum(ContourForceArray(field_x, field_y, x, y).^2, 2))) / length(tmp_theta);
            
            fit1_x = xc + a * beta * cos(tmp_theta1) * cos(phi) - b * beta * sin(tmp_theta1) * sin(phi);
            fit1_y = yc + a * beta * cos(tmp_theta1) * sin(phi) + b * beta * sin(tmp_theta1) * cos(phi);
            fit1 = sum(sqrt(sum(ContourForceArray(field_x, field_y, fit1_x, fit1_y).^2, 2))) / length(tmp_theta1);
            
            fit2_x = xc + a / beta * cos(tmp_theta2) * cos(phi) - b / beta * sin(tmp_theta2) * sin(phi);
            fit2_y = yc + a / beta * cos(tmp_theta2) * sin(phi) + b / beta * sin(tmp_theta2) * cos(phi);
            fit2 = sum(sqrt(sum(ContourForceArray(field_x, field_y, fit2_x, fit2_y).^2, 2))) / length(tmp_theta2);
            
            fit_save(it) = fit0 - fit1 / 2 - fit2 / 2;
        end
    end

    if show_fitness == 1 && (~is_octave() || strcmp(get(0, 'defaultfigurevisible'), 'on'))
        figure; hold on;
        plot(1:iter, fit_save, 'b', 'LineWidth', 2);
        legend('fitness function');
        grid on;
        xlabel('iteration');
        ylabel('fitness function');
        title('fitness function in each iteration');
    end
end