% script for plotting viterbi traces

n=size(gamma,1); 
angl=zeros(ntr,n);
absz=zeros(ntr,n);
hsvimg=zeros(ntr,n,3);
vitimg=zeros(ntr,n,3);

angl2=zeros(ntr,n);
absz2=zeros(ntr,n);

delta=0;
vs=zeros(1,ns);
statecolor=zeros(1,ns);
if ~exist('showrxntimes','var')
    showrxntimes=1;
end;
if ~exist('cvector','var')
    cvector=1:ns;
elseif length(cvector)~=ns
    cvector=1:ns;
end

if ~exist('seq','var')
    mg=mean(gamma,3);
    [a,b]=find(mg(2:floor(7/10*end),:)==max(mg(2:floor(7/10*end),:)));[~,seq]=sort(a);
    sq=ceil(sqrt(ns));
end

if 0
    mplot(mean(gamma(:,seq,:),3)');
end

if 0
    for i=1:ns
        statecolor(i)=exp(2i*pi*(i+delta)/ns);
        statecolor(i)=(angle(statecolor(i))+pi)/2/pi;
    end
else
    switch ns
        case 2
            statecolor=[0.6000 0.2500];
        case 3
            statecolor=[0.60 0.250 1.000];
        case 4
            statecolor=[0.6000 0.1300 1.0000 0.2500];
        case 5
            statecolor=[0.7000 0.5000 0.1000 0.3000 0.9000];
        case 6
            %statecolor=[0.6667 0.8333  0.3333 0.1667 1.0000 0.5000];
            statecolor=[0.6667 0.8333 1.0000 0.1667 0.3333 0.5000];
    end
    statecolor=statecolor(cvector);
    for i=1:ns, vs(i)=find(seq==i);end
end

for tr=1:ntr
    for t=1:n % n: enable variable trial duration
        angl2(tr,t)=statecolor(vs(vi(tr,t))); % discrete angles from discrete color palette
        absz2(tr,t)=gamma(t,vi(tr,t),tr); % brightness prop to responsibility
        vitimg(tr,t,3)=gamma(t,vi(tr,t),tr);
        hsvimg(tr,t,3)=max(max(gamma(t,:,tr)));
        z=0;
        for j=1:ns
            z=z+gamma(t,j,tr)*exp(2i*pi*j/ns);
            angl(tr,t)=(angle(z)+pi)/2/pi;
            absz(tr,t)=abs(z);
        end
    end
    hsvimg(tr,:,1)=angl(tr,:);
    hsvimg(tr,:,2)=absz(tr,:);
    
    vitimg(tr,:,1)=angl2(tr,:);
    vitimg(tr,:,2)=absz2(tr,:);
end

if ecogdata&&showrxntimes
    for tr=1:ntr
        dt=ceil((rxntimes(tr)-offset+500)*f-.008*n):floor((rxntimes(tr)-offset+500)*f+.008*n);
        if dt>0
            hsvimg(tr,dt,1)=0;
            hsvimg(tr,dt,2)=0;
            hsvimg(tr,dt,3)=0;
            vitimg(tr,dt,1)=0;
            vitimg(tr,dt,2)=0;
            vitimg(tr,dt,3)=0;
            
            if exist('rxntimesoff','var')
                dt2=ceil((rxntimesoff(rxnfilter(tr))-offset+500)*f-.008*n):floor((rxntimesoff(rxnfilter(tr))-offset+500)*f+.008*n);
                hsvimg(tr,dt2,1)=0;
                hsvimg(tr,dt2,2)=0;
                hsvimg(tr,dt2,3)=1;
                vitimg(tr,dt2,1)=0;
                vitimg(tr,dt2,2)=0;
                vitimg(tr,dt2,3)=1;
            end
            
        end
    end
end
figure, imagesc(hsv2rgb(vitimg)), if scrambled, title(strcat(subject,' (scrambled)'));else, title(subject);end
if 0, figure, imagesc(hsv2rgb(hsvimg)), title('gamma'), end

if 0
    for tr=1:10:100
        figure;
        for k=1:ns
            plot(gamma(p:end-p,k,tr),'*'); title(tr)
            hold on
        end
        xlabel('time')
        ylabel('gamma')
        line(([rxntimes(tr),rxntimes(tr)]+500-offset)*f,[0,1],'LineWidth',2,'Color',[0 0 0]);
    end
end

if 0
    Spdc=reshape(tpdcn(x,p,pt,-nf),nd*nd,nf+1,[]); % calculate pdc from trial time series and reshape
    vecSpdc=squeeze(sum(Spdc,2)); % integrate over all frequencies
    e1=pcaplot3(vecSpdc',80,-1);
    title('PDC');
end