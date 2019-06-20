function [A,Q,G]=tmvar(PyArea,np,pt,s)
% computes MVAR parameters from trial data
%
% PyArea..........: trial data (nd x tr x nb)
% np..............: model order
% pt..............: quasi-stationary window size
% s...............: shift

if nargin<3, pt=size(PyArea,3); end
if nargin<4, s=floor(.5*pt); else, s=floor(s);end;

[nd,~,nb]=size(PyArea); % electr x trials x #dt
G=zeros(nd,nd,np+1); % initialize cross-covariances for lag zero to lag np
B=floor((nb-pt)/s)+1; % B=number of q.-stationary windows
A=zeros(nd,nd,np,B); % electr x electr x p x #dt
Q=zeros(nd,nd,B);

for b=1:B % qstationary windows
    if 0
    for p=1:np+1 % calculate cross covariances from lag 0 up to lag np
        G(:,:,p)=mean(ccov(PyArea(1:nd,:,s*(b-1)+1:s*(b-1)+pt),p-1),3); % nd x trials x pt -> mean(nd x nd x pt-p, 3)
        %G(:,:,p)=mean(ccov(dzscore(PyArea(1:nd,:,s*(b-1)+1:s*(b-1)+pt)),p-1),3); % nd x trials x pt -> mean(nd x nd x pt-p, 3)        
    end
    res=regress(G); % compute MVAR coefficient matrix A and covariance Q
    A(:,:,:,b)=res{1}; % set of A's for current time window
    Q(:,:,b)=res{3}; % MVAR noise covariance
    elseif 1
        [x,y]=memsTest(PyArea(1:nd,:,s*(b-1)+1:s*(b-1)+pt),np,0);
        A(:,:,:,b)=reshape(x{np},nd,nd,np);
        Q(:,:,b)=y;
    else
        [A(:,:,:,b),Q(:,:,b)]=regress2(PyArea(1:nd,:,s*(b-1)+1:s*(b-1)+pt),np);
    end
end

end
