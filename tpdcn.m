function xpdc=tpdcn(PyArea,np,pt,nf,s)
% computes partial directed coherence from trial data
%
% PyArea..........: trial data (nd x tr x nb)
% np..............: model order
% pt..............: quasi-stationary window size
% fmin, fmax......: min and max frequency
% nf..............: number of frequency intervals
% s...............: shift (default: pt/2)

if nargin<5, s=floor(pt/2);else, s=floor(s);end;

[nd,~,nb]=size(PyArea); % electr x trials x #dt
B=floor((nb-pt)/s)+1; % B=number of q.-stationary windows
if nf<0
    sym=-.5;
    nf=-nf;
    xpdc=zeros(nd,nd,nf+1,B);
else
    sym=0;
    xpdc=zeros(nd,nd,ceil(nf/2)+1,B);
end

A=tmvar(PyArea,np,pt,s); % estimate MVAR parameters from trial data

for b=1:B % q.-stationary windows
    i=0;
    for f=sym:1/nf:.5
        i=i+1;
        xpdc(:,:,i,b)=pdc(A(:,:,:,b),f,1);
    end
end

end
