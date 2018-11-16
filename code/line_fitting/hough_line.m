function [theta_Hough,s_Hough]=hough_line(I)

[H, theta_Hough, s_Hough] = hough(I>0);
h=fspecial('gaussian',10,3);
H2=imfilter(H,h);
[~,maxindex]=max(H2(:));
theta_Hough=theta_Hough(floor((maxindex-1)/size(H,1))+1)/180*pi;
s_Hough=s_Hough(mod(maxindex-1,size(H,1))+1);