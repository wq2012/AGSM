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

function [xc, yc, r]=fit_circle_force(init,increment,threshold,field_x,field_y,iter)

x0=init(1);
y0=init(2);
r0=init(3);
d_xc=increment(1);
d_yc=increment(2);
d_r=increment(3);
t_xc=threshold(1);
t_yc=threshold(2);
t_r=threshold(3);
%%  Fit a circle in the force filed. xc, yc and r will be returned.
%   The function of the circle:
%      x=xc+r*cos(theta)
%      y=yc+r*sin(theta)
%   x0, y0: initial position of center
%   r0: initial radius
%   d_xc, d_yc: increment of xc and yc in each loop
%   d_r: increment of r in each loop
%   t_xc, t_yc: threshold needed to update xc and yc
%   t_r: threshold needed to update r
%   field_x: the x component of field
%   field_y: the y component of field
%   iter: number of iterations

[m,n]=size(field_x);
% m rows and n columns, but x is column and y is row here
xc=x0;
yc=y0;
r=r0;

xc_save=zeros(iter,1);
yc_save=zeros(iter,1);
r_save=zeros(iter,1);
F_in_save=zeros(iter,1);
F_left_right_save=zeros(iter,1);
F_down_up_save=zeros(iter,1);
F_diag1_save=zeros(iter,1);
F_diag2_save=zeros(iter,1);
    
for it=1:iter
    [x,y,theta]=circle_in_image(m,n,xc,yc,r);
      
    %% F_round
    F_round=sum(ContourForceArray(field_x,field_y,x,y));
    F_round=F_round/length(theta);

    %% F_in, which is the inward force
    F_in=CircleNormalForce(field_x,field_y,x,y,theta);

    %% update xc and yc
    %F_left_right=F_left-F_right;
    F_left_right=dot(F_round,[1,0]);
    if F_left_right>t_xc
        xc=xc+d_xc;
    elseif F_left_right<-t_xc
        xc=xc-d_xc;
    end
    
    %F_down_up=F_down-F_up;
    F_down_up=dot(F_round,[0,1]);
    if F_down_up>t_yc
        yc=yc+d_yc;
    elseif F_down_up<-t_yc
        yc=yc-d_yc;
    end
    

    
    %% update xc and yc again according to diagonal force
    %F_diag1=F_downleft-F_upright;
    F_diag1=dot(F_round,[0.7071,0.7071]);
    if F_diag1>t_xc+t_yc
        xc=xc+d_xc;
        yc=yc+d_yc;
    elseif F_diag1<-t_xc-t_yc
        xc=xc-d_xc;
        yc=yc-d_yc;
    end
    
    %F_diag2=F_downright-F_upleft;
    F_diag2=dot(F_round,[-0.7071,0.7071]);
    if F_diag2>t_xc+t_yc
        xc=xc-d_xc;
        yc=yc+d_yc;
    elseif F_diag2<-t_xc-t_yc
        xc=xc+d_xc;
        yc=yc-d_yc;
    end
    
    %% update r
    if F_in>t_r
        r=r+d_r;
    elseif F_in<-t_r
        r=r-d_r;
    end
    

    
    xc_save(it)=xc;
    yc_save(it)=yc;
    r_save(it)=r;
    F_in_save(it)=abs(F_in);
    F_left_right_save(it)=abs(F_left_right);
    F_down_up_save(it)=abs(F_down_up);
    F_diag1_save(it)=abs(F_diag1);
    F_diag2_save(it)=abs(F_diag2);
end




figure;hold on;
plot(1:iter,F_in_save,'b','LineWidth',2);
plot(1:iter,F_left_right_save,'r','LineWidth',2);
plot(1:iter,F_down_up_save,'g','LineWidth',2);
plot(1:iter,F_diag1_save,'--c','LineWidth',2);
plot(1:iter,F_diag2_save,'--m','LineWidth',2);
legend('F in','F left right','F down up','F diagonal','F anti-diagonal');
