function y=bin(raw,pt,shift)
% bin raw data in windows of size pt with translation shift

if nargin<3, shift=pt; end

[nd,nt,nb0]=size(raw);
if nb0>1
    nb=floor((nb0-pt)/shift)+1;
    y=zeros(nd,nt,nb);
    for k=1:nb
        d=shift*(k-1);
        y(:,:,k)=mean(raw(:,:,d+1:d+pt),3); % binning for each trial
    end
else
    nb=floor((nt-pt)/shift)+1;
    y=zeros(nd,nb);
    for k=1:nb
        d=shift*(k-1);
        y(:,k)=mean(raw(:,d+1:d+pt),2); % binning for each trial
    end
end

end