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

function [xc, yc, D, fit_save]=fit_spline_contour_force(init,increment,threshold,bound,field_x,field_y,iter,show_fitness)

x0=init(1);
y0=init(2);
D0=init(3:end);
d_xc=increment(1);
d_yc=increment(2);
d_D=increment(3);
t_xc=threshold(1);
t_yc=threshold(2);
t_D=threshold(3);
D_low=bound(1);
D_up=bound(2);

%%  Fit a cubic spline contour in the force filed. xc, yc, and D will be returned.
%   The functions of the cubic spline contour:
%       x=xc+xspline(D(1),...,D(NN))
%       y=yc+yspline(D(1),...,D(NN))
%   x0, y0: initial position of center
%   D0: initial distances
%   d_xc, d_yc: increment of xc and yc in each loop
%   d_D: increment of D(i) in each loop
%   t_xc, t_yc: threshold needed to update xc and yc
%   t_D: threshold needed to update each D(i)
%   D_low, D_up: lower bound and upper bound of D(i)
%   field_x: the x component of force field
%   field_y: the y component of force field
%   iter: number of iterations
%   show_fitness: a flag of whether to show fitness function

NN=length(D0);

% m rows and n columns, but x is column and y is row here
[m,n]=size(field_x);

xc=x0;
yc=y0;
D=D0;

fit_save=zeros(iter,1);

for it=1:iter
    [x,y,theta]=spline_contour_in_image(m,n,xc,yc,D);
    
    %% F_around
    F_round=sum(ContourForceArray(field_x,field_y,x,y));
    F_round=F_round/length(theta);
    
    %% F_k, the force on the kth landmark
    dt=2*pi/NN; % dt: delta theta
    F_k=zeros(NN,1);
    for k=1:NN
        if k==1
            index=find(theta<dt/2 | theta>2*pi-dt/2);
        else
            index=find(theta>(k-1.5)*dt & theta<(k-0.5)*dt);
        end
        
        F_k(k)=SplineContourForce(field_x,field_y,x(index),y(index),theta(index));
        
    end
    
    %% update xc and yc
    F_left_right=dot(F_round,[1,0]);
    if F_left_right>t_xc
        xc=xc+d_xc;
    elseif F_left_right<-t_xc
        xc=xc-d_xc;
    end
    
    F_down_up=dot(F_round,[0,1]);
    if F_down_up>t_yc
        yc=yc+d_yc;
    elseif F_down_up<-t_yc
        yc=yc-d_yc;
    end
    
    %% update xc and yc again according to diagonal force
    F_diag1=dot(F_round,[0.7071,0.7071]);
    if F_diag1>t_xc+t_yc
        xc=xc+d_xc;
        yc=yc+d_yc;
    elseif F_diag1<-t_xc-t_yc
        xc=xc-d_xc;
        yc=yc-d_yc;
    end
    
    F_diag2=dot(F_round,[-0.7071,0.7071]);
    if F_diag2>t_xc+t_yc
        xc=xc-d_xc;
        yc=yc+d_yc;
    elseif F_diag2<-t_xc-t_yc
        xc=xc+d_xc;
        yc=yc-d_yc;
    end
    
    %% update and restrict D
    for k=1:NN
        if F_k(k)>t_D
            D(k)=D(k)+d_D;
        elseif F_k(k)<-t_D
            D(k)=D(k)-d_D;
        end
        
        if D(k)>D_up
            D(k)=D_up;
        end
        if D(k)<D_low
            D(k)=D_low;
        end
    end
    
    %% fitness function
    if show_fitness==1
        beta=0.9;
        
        [x,y,theta]=spline_contour_in_image(m,n,xc,yc,D);
        [x1,y1,theta1]=spline_contour_in_image(m,n,xc,yc,D*beta);
        [x2,y2,theta2]=spline_contour_in_image(m,n,xc,yc,D/beta);
        
        fit0=ContourForceArray(field_x,field_y,x,y);
        fit0=sum(sqrt(sum(fit0.*fit0,2)))/length(theta);
        
        fit1=ContourForceArray(field_x,field_y,x1,y1);
        fit1=sum(sqrt(sum(fit1.*fit1,2)))/length(theta1);
        
        fit2=ContourForceArray(field_x,field_y,x2,y2);
        fit2=sum(sqrt(sum(fit2.*fit2,2)))/length(theta2);
        
        fit_save(it)=fit0-fit1/2-fit2/2;
    end
    
end

if show_fitness==1
    figure;hold on;
    plot(1:iter,fit_save,'b','LineWidth',2);
    legend('fitness function');
    grid on;
    xlabel('iteration');
    ylabel('fitness function');
    title('fitness function in each iteration');
end
