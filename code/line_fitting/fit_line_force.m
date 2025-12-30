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

function [theta, s] = fit_line_force(init, increment, threshold, field_x, field_y, iter)
% FIT_LINE_FORCE Fit a line in the force field.
%    [theta, s] = fit_line_force(init, increment, threshold, field_x, 
%    field_y, iter) iteratively fits a line to a vector field.
%    The line equation is: x * cos(theta) + y * sin(theta) = s
%
%    Inputs:
%        init - [theta0, s0] initial parameters.
%        increment - [d_theta, d_s] step sizes for updates.
%        threshold - [t_theta, t_s] thresholds for updates.
%        field_x - The x component of the vector field.
%        field_y - The y component of the vector field.
%        iter - Number of iterations.
%
%    Outputs:
%        theta - Final line angle.
%        s - Final line offset.

    theta0 = init(1);
    s0 = init(2);
    d_theta = increment(1);
    d_s = increment(2);
    t_theta = threshold(1);
    t_s = threshold(2);

    [m, n] = size(field_x);
    theta = theta0;
    s = s0;

    force_save = zeros(iter, 1);
    torque_save = zeros(iter, 1);
    theta_save = zeros(iter, 1);
    s_save = zeros(iter, 1);

    for it = 1:iter
        theta = mod(theta, 2 * pi);
        [x, y] = LineInImage(m, n, theta, s);
        [torque, k] = LineTorque(field_x, field_y, x, y, theta);

        % Update theta
        if mod(theta, 2 * pi) <= pi || x(1) == x(end)
            if torque > t_theta
                theta = theta - d_theta;
            elseif torque < -t_theta
                theta = theta + d_theta;
            end
        else
            if torque > t_theta
                theta = theta + d_theta;
            elseif torque < -t_theta
                theta = theta - d_theta;
            end
        end

        s = x(k) * cos(theta) + y(k) * sin(theta);
        torque_save(it) = abs(torque);

        % Update distance (s)
        [x, y] = LineInImage(m, n, theta, s);
        force_normal = LineNormalForce(field_x, field_y, x, y, theta);

        if force_normal > t_s
            s = s + d_s;
        elseif force_normal < -t_s
            s = s - d_s;
        end

        if theta > pi
            theta = theta - pi;
            s = -s;
        end

        force_save(it) = abs(force_normal);
        theta_save(it) = theta;
        s_save(it) = s;
    end

    % Visualization (only if figures are visible)
    if ~is_octave() || strcmp(get(0, 'defaultfigurevisible'), 'on')
        figure; hold on;
        plot(1:iter, torque_save, 'b', 'LineWidth', 2);
        plot(1:iter, force_save, '-.r', 'LineWidth', 2);
        grid on;
        legend('absolute value of torque', 'absolute value of force');
        xlabel('iteration');
        ylabel('torque and force');
        title('torque and force in each iteration');

        figure; hold on;
        plot(1:iter, theta_save, 'b', 'LineWidth', 2);
        legend('theta');
        ylabel('theta');
        xlabel('iteration');
        title('theta in each iteration');

        figure; hold on;
        plot(1:iter, s_save, 'r', 'LineWidth', 2);
        legend('s');
        ylabel('s');
        xlabel('iteration');
        title('s in each iteration');
    end
end
