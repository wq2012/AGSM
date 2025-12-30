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

function [xc, yc, r] = fit_circle_force(init, increment, threshold, field_x, field_y, iter)
% FIT_CIRCLE_FORCE Fit a circle in the force field.
%    [xc, yc, r] = fit_circle_force(init, increment, threshold, field_x, 
%    field_y, iter) iteratively fits a circle to a vector field.
%
%    Inputs:
%        init - [x0, y0, r0] initial parameters.
%        increment - [d_xc, d_yc, d_r] step sizes.
%        threshold - [t_xc, t_yc, t_r] update thresholds.
%        field_x - The x component of the field.
%        field_y - The y component of the field.
%        iter - Number of iterations.
%
%    Outputs:
%        xc - Final center x-coordinate.
%        yc - Final center y-coordinate.
%        r - Final radius.

    x0 = init(1);
    y0 = init(2);
    r0 = init(3);
    d_xc = increment(1);
    d_yc = increment(2);
    d_r = increment(3);
    t_xc = threshold(1);
    t_yc = threshold(2);
    t_r = threshold(3);

    [m, n] = size(field_x);
    xc = x0;
    yc = y0;
    r = r0;

    xc_save = zeros(iter, 1);
    yc_save = zeros(iter, 1);
    r_save = zeros(iter, 1);
    F_in_save = zeros(iter, 1);
    F_left_right_save = zeros(iter, 1);
    F_down_up_save = zeros(iter, 1);
    F_diag1_save = zeros(iter, 1);
    F_diag2_save = zeros(iter, 1);
        
    for it = 1:iter
        [x, y, theta] = circle_in_image(m, n, xc, yc, r);
          
        % Calculate round force
        F_round = sum(ContourForceArray(field_x, field_y, x, y));
        F_round = F_round / length(theta);

        % Calculate inward force
        F_in = CircleNormalForce(field_x, field_y, x, y, theta);

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
        
        % Update center again with diagonal force
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
        
        % Update radius (r)
        if F_in > t_r
            r = r + d_r;
        elseif F_in < -t_r
            r = r - d_r;
        end
        
        xc_save(it) = xc;
        yc_save(it) = yc;
        r_save(it) = r;
        F_in_save(it) = abs(F_in);
        F_left_right_save(it) = abs(F_left_right);
        F_down_up_save(it) = abs(F_down_up);
        F_diag1_save(it) = abs(F_diag1);
        F_diag2_save(it) = abs(F_diag2);
    end

    % Visualization (only if figures are visible)
    if ~is_octave() || strcmp(get(0, 'defaultfigurevisible'), 'on')
        figure; hold on;
        plot(1:iter, F_in_save, 'b', 'LineWidth', 2);
        plot(1:iter, F_left_right_save, 'r', 'LineWidth', 2);
        plot(1:iter, F_down_up_save, 'g', 'LineWidth', 2);
        plot(1:iter, F_diag1_save, '--c', 'LineWidth', 2);
        plot(1:iter, F_diag2_save, '--m', 'LineWidth', 2);
        legend('F in', 'F left right', 'F down up', 'F diagonal', 'F anti-diagonal');
        title('Circle fitting forces in each iteration');
        xlabel('Iteration');
        ylabel('Force magnitude');
    end
end
