function b=bAR2(mu,B,Omega,X,p) % mean, regression matrix, precision, data, model order
% calculates p(x(t)|x(t-1),s(t)) for all s(t)
% mu: [nd+nd*(np-1)] x ns = np*nd x ns
% B: (np*nd) x nd
% X: (np*nd) x n

[nd,ns]=size(mu);
n=size(X,2);
b=zeros(ns,n);
nd=floor(nd/p);
for a=1:ns
    % initialization
    x=X(:,1)-mu(:,a);
    b(a,1)=exp(-.5*x'*Omega(:,:,a)*x);
    for t=2:n
        x=X(:,t)-mu(:,a)-B(:,:,a)*X(:,t-1); % P(x(t)|x(t-1),s(t)=a)
        b(a,t)=exp(-.5*x'*Omega(:,:,a)*x); % for P(s(t)|X(t+1..T),M): X given
    end
    % normalization
    %b(a,:)=b(a,:)/sqrt((2*pi)^nd*det(S(1:nd,1:nd,a)));
    b(a,:)=b(a,:)/sqrt((2*pi)^nd)*sqrt(det(Omega(1:nd,1:nd,a)));
end
end