function C=ccov(x,l)
% Cross Covariance for trial data: C=<x(t)x(t+l)>'=<x(t)x(t-l)>
% x: dimensions x trials x time
% l: lag (cov(x,-l)=cov(x,l)')

[c,tr,n]=size(x);
C=zeros(c,c,n-l);
for i=1:n-l
    for c1=1:c
        d1=squeeze(x(c1,:,i)-mean(x(c1,:,i),2));
        for c2=1:c
            C(c2,c1,i)=d1*squeeze(x(c2,:,i+l)-mean(x(c2,:,i+l),2))';
        end
    end
end
C=C/(tr-1);

end