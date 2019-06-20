function [ data ] = tySampling( s, C, R)
% generates data set of length n and dimension nd from trial data
% A: state transition matrix: nd x nd*np
% C: measurement matrix
% Q: covariance matrix for states
% R: observation noise
% output: [signal measurement]

[nd,nt,n]=size(s);

try
    [~,~,Rt]=size(R);
    if nd>1, CholR=chol(R(:,:,1)); else  CholR=sqrt(R(:,:,1)); end;
catch
    Rt=0;
    if nd>1, CholR=chol(R); else  CholR=sqrt(R); end;
end;

y=zeros(nd,nt,n); % data matrix

% draw samples
l=1; % initial data point
if ~(Rt<l)
    if nd>1, CholR=chol(R(:,:,l)); else  CholR=sqrt(R(:,:,l)); end;
end;

for l=1:n
    if ~(Rt<l)
        if nd>1, CholR=chol(R(:,:,l)); else  CholR=sqrt(R(:,:,l)); end;
    end;
    for tr=1:nt
    y(:,tr,l)=C*s(:,tr,l) + CholR*randn(nd,1); % noisy observation
    end
end;

data=y;
end