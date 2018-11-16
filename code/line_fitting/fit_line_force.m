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

function [theta, s]=fit_line_force(init,increment,threshold,field_x,field_y,iter)

theta0=init(1);
s0=init(2);
d_theta=increment(1);
d_s=increment(2);
t_theta=threshold(1);
t_s=threshold(2);
%%  Fit a line in the force filed. Theta and s will be returned.
%   The function of the line: x * cos(theta) + y * sin(theta) = s
%   theta0: initial value of theta
%   s0: initial value of s
%   d_theta: increment of theta in each loop
%   d_s: increment of s in each loop
%   t_theta: threshold needed to update theta
%   t_s: threshold needed to update s
%   field_x: the x component of field
%   field_y: the y component of field
%   iter: number of iterations

[m,n]=size(field_x);
% m rows and n columns, but x is column and y is row here
theta=theta0;
s=s0;

force_save=zeros(iter,1);
torque_save=zeros(iter,1);
theta_save=zeros(iter,1);
s_save=zeros(iter,1);

for it=1:iter
    
    theta=mod(theta,2*pi);
    
    [x,y]=LineInImage(m,n,theta,s);

    [torque,k]=LineTorque(field_x,field_y,x,y,theta);

    
    % update theta and s
    if mod(theta,2*pi)<=pi || x(1)==x(end)
        if torque>t_theta
            theta=theta-d_theta;
            
        elseif torque<-t_theta
            theta=theta+d_theta;
            
        end
    else
        if torque>t_theta
            theta=theta+d_theta;
            
        elseif torque<-t_theta
            theta=theta-d_theta;
        end
    end
    
    s=x(k)*cos(theta)+y(k)*sin(theta);
    
    torque_save(it)=abs(torque);
    
    %% update distance
    [x,y]=LineInImage(m,n,theta,s);

    % the force on the line
    force_normal=LineNormalForce(field_x,field_y,x,y,theta);

    if force_normal>t_s
        s=s+d_s;
    elseif force_normal<-t_s
        s=s-d_s;
    end

    if theta>pi
        theta=theta-pi;
        s=-s;
    end
    
    force_save(it)=abs(force_normal);
    
    theta_save(it)=theta;
    s_save(it)=s;
end


figure;hold on;
plot(1:iter,torque_save,'b','LineWidth',2);
plot(1:iter,force_save,'-.r','LineWidth',2);
grid on;
legend('absolute value of torque','absulote value of force');
xlabel('iteration');
ylabel('torque and force');
title('torque and force in each iteration');

figure;hold on;
plot(1:iter,theta_save,'b','LineWidth',2);
legend('theta');
ylabel('theta');
xlabel('iteration');
title('theta in each iteration');
figure;hold on;
plot(1:iter,s_save,'r','LineWidth',2);
legend('s');
ylabel('s');
xlabel('iteration');
title('s in each iteration');
