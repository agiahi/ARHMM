ntr=size(vi,1);
nt=size(vi,2);
cbar=0;
Res=zeros(nd,ntr,nt-1);
Res0=zeros(nd,ntr,nt-1);
for tr=1:ntr
    for t=2:nt
        Res(:,tr,t-1)=X(1:nd,t,tr)-B(1:nd,:,vi(tr,t))*X(:,t-1,tr);
        Res0(:,tr,t-1)=X(1:nd,t-1,tr);
    end;
end;

cR=zeros(nd,nd,nt-1,p);
cR0=zeros(nd,nd,nt-1,p);
cc=ccov(Res,0); cRab=diag(mean(cc,3)).^.5; cRab=cRab*cRab';
cc=ccov(Res0,0); cR0ab=diag(mean(cc,3)).^.5; cR0ab=cR0ab*cR0ab';
for l=1:p
    cR(:,:,1:nt-l-1,l)=ccov(Res,l); 
    cR0(:,:,1:nt-l-1,l)=ccov(Res0,l); 
end
if 0
    lim=1;
else
    lim=max(max(max(mean(cR0(nd:-1:1,1:nd,:,:)./cR0ab(nd:-1:1,:),3))));
end

if 0
figure;subplot(2,1,1)
axis off
mplot(mean(cR0(nd:-1:1,1:nd,:,:),3)./cR0ab(nd:-1:1,:),1,lim); if cbar,colorbar;end; title('raw temp correlation')
subplot(2,1,2)
axis off
mplot(mean(cR(nd:-1:1,1:nd,:,:),3)./cRab(nd:-1:1,:),1,lim); if cbar,colorbar;end; title('MVAR temp correlation')
end

if 1
figure;
for k=1:p
subplot(2,p,k)
axis off;
mplot(mean(cR0(nd:-1:1,1:nd,:,k),3)./cR0ab(nd:-1:1,:),1,lim); if cbar,colorbar;end; title(sprintf('data corr (lag=%i)',k))
subplot(2,p,k+p)
axis off;
mplot(mean(cR(nd:-1:1,1:nd,:,k),3)./cRab(nd:-1:1,:),1,lim); if cbar,colorbar;end; title(sprintf('residual corr (lag=%i)',k))
end
end