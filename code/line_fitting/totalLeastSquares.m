function [theta_TLS, s_TLS, error_TLS]=totalLeastSquares(x,y)

x_bar=mean(x);
y_bar=mean(y);
x2_bar=mean(x.^2);
y2_bar=mean(y.^2);
xy_bar=mean(x.*y);
X=[x2_bar-x_bar^2,xy_bar-x_bar*y_bar;xy_bar-x_bar*y_bar,y2_bar-y_bar^2];
[V,D] = eig(X);
a1=V(1,1);a2=V(2,1);
b1=V(2,1);b2=V(2,2);
c1=-a1*x_bar-b1*y_bar;c2=-a2*x_bar-b2*y_bar;
error_TLS1=mean((x.*a1+y.*b1+c1).^2);
error_TLS2=mean((x.*a2+y.*b2+c2).^2);

if error_TLS1<error_TLS2
    error_TLS=error_TLS1;
    a=a1;b=b1;c=c1;
else
    error_TLS=error_TLS2;
    a=a2;b=b2;c=c2;
end
theta_TLS=atan(b/a);
if a<0
    theta_TLS=theta_TLS+pi;
end
s_TLS=-c;

