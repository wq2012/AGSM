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

function [xc, yc, D, fit_save] = fit_spline_contour_force(init, increment, threshold, bound, field_x, field_y, iter, show_fitness)
% FIT_SPLINE_CONTOUR_FORCE Fit a cubic spline contour in the force field.
%    [xc, yc, D, fit_save] = fit_spline_contour_force(init, increment, 
%    threshold, bound, field_x, field_y, iter, show_fitness) iteratively 
%    fits a spline contour to a vector field.
%
%    Inputs:
%        init - [x0, y0, D0...] initial center and landmark distances.
%        increment - [d_xc, d_yc, d_D] step sizes for updates.
%        threshold - [t_xc, t_yc, t_D] update thresholds.
%        bound - [D_low, D_up] constraints on landmark distances.
%        field_x - The x component of force field.
%        field_y - The y component of force field.
%        iter - Number of iterations.
%        show_fitness - Flag to toggle fitness function tracking.
%
%    Outputs:
%        xc, yc - Final center coordinates.
%        D - Final landmark distances.
%        fit_save - Fitness values over iterations.

    x0 = init(1);
    y0 = init(2);
    D0 = init(3:end);
    d_xc = increment(1);
    d_yc = increment(2);
    d_D = increment(3);
    t_xc = threshold(1);
    t_yc = threshold(2);
    t_D = threshold(3);
    D_low = bound(1);
    D_up = bound(2);

    NN = length(D0);
    [m, n] = size(field_x);

    xc = x0;
    yc = y0;
    D = D0;

    fit_save = zeros(iter, 1);

    for it = 1:iter
        [x, y, theta] = spline_contour_in_image(m, n, xc, yc, D);
        
        % Average force around contour
        F_round = sum(ContourForceArray(field_x, field_y, x, y));
        F_round = F_round / length(theta);
        
        % Discrete forces on each landmark
        dt = 2 * pi / NN;
        F_k = zeros(NN, 1);
        for k = 1:NN
            if k == 1
                idx = find(theta < dt / 2 | theta > 2 * pi - dt / 2);
            else
                idx = find(theta > (k - 1.5) * dt & theta < (k - 0.5) * dt);
            end
            F_k(k) = SplineContourForce(field_x, field_y, x(idx), y(idx), theta(idx));
        end
        
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
        
        % Diagonal updates
        F_diag1 = dot(F_round, [0.7071, 0.7071]);
        if F_diag1 > t_xc + t_yc
            xc = xc + d_xc; yc = yc + d_yc;
        elseif F_diag1 < -t_xc - t_yc
            xc = xc - d_xc; yc = yc - d_yc;
        end
        
        F_diag2 = dot(F_round, [-0.7071, 0.7071]);
        if F_diag2 > t_xc + t_yc
            xc = xc - d_xc; yc = yc + d_yc;
        elseif F_diag2 < -t_xc - t_yc
            xc = xc + d_xc; yc = yc - d_yc;
        end
        
        % Update landmark distances
        for k = 1:NN
            if F_k(k) > t_D
                D(k) = D(k) + d_D;
            elseif F_k(k) < -t_D
                D(k) = D(k) - d_D;
            end
            D(k) = max(D_low, min(D(k), D_up));
        end
        
        % Performance/fitness tracking
        if show_fitness == 1
            beta = 0.9;
            [~, ~, tmp_theta] = spline_contour_in_image(m, n, xc, yc, D);
            [~, ~, tmp_theta1] = spline_contour_in_image(m, n, xc, yc, D * beta);
            [~, ~, tmp_theta2] = spline_contour_in_image(m, n, xc, yc, D / beta);
            
            fit0 = sum(sqrt(sum(ContourForceArray(field_x, field_y, x, y).^2, 2))) / length(tmp_theta);
            
            fit1_dd = myspline(D * beta, tmp_theta1);
            fit1_x = xc + fit1_dd .* cos(tmp_theta1);
            fit1_y = yc + fit1_dd .* sin(tmp_theta1);
            fit1 = sum(sqrt(sum(ContourForceArray(field_x, field_y, fit1_x, fit1_y).^2, 2))) / length(tmp_theta1);
            
            fit2_dd = myspline(D / beta, tmp_theta2);
            fit2_x = xc + fit2_dd .* cos(tmp_theta2);
            fit2_y = yc + fit2_dd .* sin(tmp_theta2);
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
