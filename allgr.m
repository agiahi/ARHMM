statecolor=[0.6000 0.1300 1.0000 0.2500]; % hues
cvector=1:ns; nd=6; % parameters
sq=ceil(sqrt(ns)); 
figure;
for k=1:ns
    ax1=subplot(sq,sq,k);
    axis([-1.25 1.25,-1 1]); axis off
    ppdc=zeros(nd,nd,ns);for pp=1:ns,for ff=0:1/100:.5,ppdc(:,:,pp)=ppdc(:,:,pp)+pdc(reshape(B(1:nd,:,pp),nd,nd,p),ff,1);end;end;
    ps=ppdc(:,:,seq(k))/ns;
    matrixToGraph2(ps,.7,max(ceil(20/sq),8),statecolor(cvector(k)),1)
    title(k)
end
axes(ax1);