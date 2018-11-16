function [theta,s]=ransac_TLS(x,y)

trial=100;
train=0.5;
accept=0.2;
reuse=0.5;

N=length(x);
Ntrain=round(N*train);

vote=zeros(N,1);

for t=1:trial
    rp=randperm(N);
    train_index=rp(1:Ntrain);
    test_index=rp(Ntrain+1:end);
    
    [theta_TLS, s_TLS]=totalLeastSquares(x(train_index),y(train_index));
    test_error=abs( x(test_index)*cos(theta_TLS)+y(test_index)*sin(theta_TLS)-s_TLS );
    [~,index]=sort(test_error);
    index=index(1:round(length(index)*accept));
    index=rp(index+Ntrain);
    vote(index)=vote(index)+1;
end

[~,index]=sort(vote,'descend');
index=index(1:round(N*reuse));
[theta, s]=totalLeastSquares(x(index),y(index));

