%
% Autoregressive Hidden Markov Model (aHMM) for ECoG data
% uses MVAR code (Morf et al.) written by Meagan Whaley for initialization
%
% created 2017-2018
%

time=cputime;
dplay=0; % show graphics
plotL=0; % plot likelihood over iterations

if crossval
    if nq==1
        rxntimes=rxntimes(trials_test);
        rxntimesoff=rxntimesoff(trials_test);
        yt=yt(:,trials_test,:);
    else
        rxntimes=rxntimes(trials_train);
        rxntimesoff=rxntimesoff(trials_train);
        yt=yt(:,trials_train,:);
    end
end

if ~0
    if ~exist('ecogdata','var'), ecogdata=1;end;
    if ecogdata
        Art=[583 1897];
        %Art=[599 2501];
        Art1=Art(1); Art2=Art(2);
        rxnfilter=find((Art1<rxntimes).*(rxntimes<Art2)==1); % trials with rxntimes between Art1 and Art2

        ff=D.sampHz/1000;
        pt=tau*ff; % 1pt=.5ms @ 2000Hz
        f=ff/pt; % effective sampling frequency (mHz)
        
        if 1
            yt=yt(:,rxnfilter,:);
            rxntimesoff=rxntimesoff(rxnfilter);
            [rxntimes,rxnfilter]=sort(rxntimes(rxnfilter)); % use only trials with rxntimes in specified range (Art=[...]) and sort trials by rxntime
        else
            w=get_words(paths,subject);
            trialind=trialind(rxnfilter);
            w=w(trialind);
            w=extractword(w);
            v=wordfreq(paths,w);
            rxnfilter=rxnfilter(v<97565);
            v=v(v<97565);
            yt=yt(:,rxnfilter,:);
            rxntimesoff=rxntimesoff(rxnfilter);
            rxntimes=rxntimes(rxnfilter);
            [v,rxnfilter]=sort(v);
        end
        y=bin(yt(:,rxnfilter,:),pt,pt);
        [nd,~,n]=size(y);
        x=y;
        
    else
        u=.1; z2=0;
        dw=10; % quasistationary window size for AMVAR
        pt=1; % bin size for data, pt/bin
        tau=1;
        fR=1/4; % observation noise 1/fR ~ SNR, fR 1/4 => SNR 1.45 (for 'data606')
        ns=3; % # of states
        p=2; % model order
        fprintf('\n simulated data')
        load 'data606'; % simulated data set used for the paper: 3 states with variable onset and model order 1
        %load 'dataej3r'; % 3 states with state sequence: 1 -> 2 -> 3 -> 1
        xt=xt(:,:,:);
        nd=size(xt,1); % # of channels
        x=tySampling(xt,eye(nd),fR*eye(nd)); % add white noise to data set
        %pause;
        rxntimes=sum(st==1,2); % define reaction times by duration of state 'st'
        [rxntimes,rxnfilter]=sort(rxntimes);
        x=x(:,rxnfilter,:); % sort trials by reaction time
        x=bin(x,pt,pt);
        n=size(x,3); % # of time bins
    end
% data preprocessing
if z2
    %x=zscore(x,0,3);
    x=zscore(x,0,2);
else
    x=x-mean(x,2);
end
else
    showrxntimes=0;
end

ndp=nd*p;
%pause;
% estimate A and Q from data
if ~(crossval&&(nq==1))
    fprintf('\n initializing dynamical parameters')
    if 1 % use AMVAR to initialize ARHMM (default)
        if ecogdata
            idx=find((600<rxntimes).*(rxntimes<1800)==1);
            %xx=bin(x(:,idx,:),ceil(100/tau),ceil(100/tau),0);
            %xx=xx-mean (xx,2);
            %dw=ceil(50/tau*D.sampHz);
        else
            scrambled=0; % no scrambled option for simulated data set
        end
        [AA,QQ]=tmvar(x,p,dw); AA=-AA; % estimate initial parameter values from trial ensemble
        if cluster % k-means clustering of AMVAR estimates
            nf=100;
            Spdc=reshape(tpdcn(x,p,dw,-nf),nd*nd,nf+1,[]); % calculate pdc from trial time series and reshape
            vecSpdc=squeeze(sum(Spdc,2)); % integrate over all frequencies
            [ikx]=kmeans(vecSpdc',ns); % k-means clustering with Euclidean distance metric
            
            kst=zeros(1,ns);for k=1:ns,kst(k)=((ikx==k)'*(1:size(ikx,1))')/sum(ikx==k);end
            [~,kseq]=sort(kst); kseq=[kseq(end) kseq(1:end-1)]'; ikx0=ikx;
            kor=zeros(1,ns);for k=1:ns, kor(k)=find(kseq==k);end; ikx=kor(ikx);
            
        end
    else % initialize ARHMM with random seeds
        AA=randn(nd,nd,p,50)/p/nd;
        QQ=zeros(nd,nd,1,50);QQ(:,:,1,:)=repmat(eye(nd)+.1,1,1,50);
    end
end
fprintf('\n inference ')

if ~(crossval&&(nq==1))
    % initialize Phi and pr
    Phi0=ones(ns)/ns;Phi=Phi0; % flat prior
    pr=ones(ns,1)/ns;
    
    % initialize state means with random seeds
    mu=zeros(ndp,ns); % augmented random seed mean vector
    Mu=zeros(nd,ns); % random seeds for state means 
    if 1
        for i=1:ns
            Mu(:,i)=randn(1,nd)*.2;
            mu(:,i)=[Mu(:,i);zeros(nd*(p-1),1)]; % random seed augmentation
        end
    end
    m=mu; % state mean, initialized with augmented random seed mean
    
    % initiate (augmented) B
    Baugm=[eye((p-1)*nd),zeros((p-1)*nd,nd)]; % augmentation to AR matrix
    S=zeros(ndp,ndp,ns); B=zeros(ndp,ndp,ns);
    nB=size(AA,4);
    for i=1:ns
        if ~cluster
            B(1:nd,:,i)=reshape(AA(:,:,1:p,ceil((i-.5)*nB/ns)),nd,ndp);
            S(1:nd,1:nd,i)=QQ(:,:,ceil((i-.5)*nB/ns));
        else
            B(1:nd,:,i)=reshape(mean(AA(:,:,1:p,ikx==i),4),nd,ndp);
            S(1:nd,1:nd,i)=mean(QQ(:,:,ikx==i),3);
        end
        if p>1
            B(1+nd:p*nd,1:p*nd,i)=Baugm; % augmentation
        end
    end
    B0=B; % initial estimates of dyn. parameters
    S0=S;
    
    Omega=zeros(ndp,ndp,ns);
    for i=1:ns
        Omega(1:nd,1:nd,i)=eye(nd)/(S(1:nd,1:nd,i));
    end
end

if 1
    % combine trials (individual E-steps)
    trs=1:size(x,2); % trials to be combined
    ntr=length(trs); % # of trials selected
    gamma=zeros(n,ns,ntr);
    xi=zeros(n-1,ns,ns,ntr);
    X=zeros(ndp,n,ntr);
    ttr=zeros(1,ntr);
    for trial=1:ntr
        if adaptcut
            ttr(trial)=ceil(min(max((rxntimes(trs(trial))+1000)/tau,.45*n),n));
        else
            ttr(trial)=n;
        end
        yy=x(:,trs(trial),1:ttr(trial));
        for t=p:size(yy,3)
            X(:,t-p+1,trial)=reshape(yy(:,1,t:-1:t-p+1),ndp,1);
        end
    end
else
    % concatenate trials (one E-step)
    trs=[1:size(x,2)]; % number of concatenated trials
    nt=length(trs); n=nt*n; n=n-p+1; ntr=1;
    yyt=zscore(x(:,1,:),0,3);
    for trial=1:nt
        x1=zscore(x(:,trs(trial),:),0,3);
        yyt=cat(3,yyt,x1);
    end
    gamma=zeros(n,ns,ntr);
    xi=zeros(n-1,ns,ns,ntr);
    X=zeros(ndp,n,ntr);
    for t=p:n
        X(:,t-p+1,1)=reshape(yyt(:,1,t:-1:t-p+1),ndp,1);
    end
    x=yyt;
end

state=ones(1,n-p+1);
vi=zeros(ntr,n);
vmax=zeros(ntr,n);
l=zeros(nq,ntr);

% EM
%###########################################################################
for q=1:nq %EM-loop with nq iterations
    
    % E-step
    % alpha-beta recursion
    %  [alpha,beta,gamma,xi,l(q,trial)]=alphabetan(pr, Phi, bAR(m, B, S, X,p));
    for trial=1:ntr
        %[gamma(:,:,trial),xi(:,:,:,trial),l(q,trial)]=alphabetan(pr, Phi, bAR2(m, B, Omega, X(:,:,trial),p));
        if 1
            if q<nq
                [gamma(1:ttr(trial),:,trial),xi(1:ttr(trial)-1,:,:,trial),l(q,trial)]=alphabetavg(pr, Phi, bAR2(m, B, Omega, X(:,1:ttr(trial),trial),p));
            else
                [gamma(1:ttr(trial),:,trial),xi(1:ttr(trial)-1,:,:,trial),l(q,trial),vi(trial,1:ttr(trial)),vmax(trial,1:ttr(trial))]=alphabetavg(pr, Phi, bAR2(m, B, Omega, X(:,1:ttr(trial),trial),p));
            end
        else
            [gamma(1:ttr(trial),:,trial),xi(1:ttr(trial)-1,:,:,trial),l(q,trial)]=alphabetav(pr, Phi, bAR2(m, B, Omega, X(:,1:ttr(trial),trial),p));
        end
    end
    
    % M-step
    % mean and ovariance update
    m=zeros(ndp,ns);
    mn=zeros(ndp,ns);
    mm=zeros(ndp,ns);
    if ntr==1
        Sg=sum(gamma);
    else
        Sg=sum(sum(gamma,3),1);
    end
    S2=zeros(ndp,ndp,ns);
    Omega=zeros(ndp,ndp,ns);
    
    for i=1:ns % states
        if 0
            % state mean estimate
            for trial=1:ntr
                for t=1:n % time
                    mn(:,i)=mn(:,i)+gamma(t,i,trial)*X(:,t,trial);
                    if t>1
                        mm(:,i)=mm(:,i)+gamma(t,i,trial)*X(:,t-1,trial);
                    end;
                end;
            end;
            m(:,i)=(mn(:,i)-B(:,:,i)*mm(:,i))/Sg(i);
        else
            % state mean estimate
            for trial=1:ntr
                for t=1:n % time
                    mn(1:ndp-1,i)=mn(1:ndp-1,i)+gamma(t,i,trial)*X(1:ndp-1,t,trial); % <x_t>_st
                    if t>1
                        mm(1:ndp-1,i)=mm(1:ndp-1,i)+gamma(t,i,trial)*X(1:ndp-1,t-1,trial); % <x_t-1>_st
                    end;
                end;
            end;
            m(1:nd,i)=(mm(1:nd,i)-B(1:nd,1:ndp-1,i)*mn(1:ndp-1,i))/Sg(i); % x_t = B_i * x_t-1 + m_i
        end
        
        % covariance estimate
        if 0
            for trial=1:ntr
                if ntr==1
                    S(:,:,i)=gamma(1,i,trial)*(X(:,1,trial)-m(:,i))*(X(:,1,trial)-m(:,i))';
                else
                    S(:,:,i)=S(:,:,i)+gamma(1,i,trial)*(X(:,1,trial)-m(:,i))*(X(:,1,trial)-m(:,i))';
                end
                for t=2:n
                    S(:,:,i)=S(:,:,i)+gamma(t,i,trial)*(X(:,t,trial)-m(:,i)-B(:,:,i)*X(:,t-1,trial))*(X(:,t,trial)-m(:,i)-B(:,:,i)*X(:,t-1,trial))';
                end
            end
            S(:,:,i)=S(:,:,i)/Sg(i);
        else
            mi=m(1:nd,i);
            for trial=1:ntr
                if ntr==1
                    S(1:nd,1:nd,i)=gamma(1,i,trial)*(X(1:nd,1,trial)-mi)*(X(1:nd,1,trial)-mi)';
                else
                    S(1:nd,1:nd,i)=S(1:nd,1:nd,i)+gamma(1,i,trial)*(X(1:nd,1,trial)-mi)*(X(1:nd,1,trial)-mi)';
                end
                for t=2:n
                    S(1:nd,1:nd,i)=S(1:nd,1:nd,i)+gamma(t,i,trial)*(X(1:nd,t,trial)-mi-B(1:nd,1:ndp,i)*X(1:ndp,t-1,trial))*(X(1:nd,t,trial)-mi-B(1:nd,1:ndp-1,i)*X(1:ndp-1,t-1,trial))';
                end
            end
            S(:,:,i)=S(:,:,i)/Sg(i);
        end
        
        
        % dynamics estimate
        B(:,:,i)=zeros(ndp,ndp);
        for trial=1:ntr
            for t=2:n
                S2(:,:,i)=S2(:,:,i)+gamma(t,i,trial)*X(:,t-1,trial)*X(:,t-1,trial)';
                B(:,:,i)=B(:,:,i)+gamma(t,i,trial)*(X(:,t,trial)-m(:,i))*X(:,t-1,trial)';
            end
        end
        if S2(:,:,i)==0
            B(:,:,i)=zeros(ndp);
        else
            B(:,:,i)=B(:,:,i)/S2(:,:,i);
        end
        
        Omega(1:nd,1:nd,i)=eye(nd)/(S(1:nd,1:nd,i));
        if p<0
            B(1+nd:ndp,1:ndp,i)=Baugm; % augmentation
        end
        
        
    end
    
    % re-estimating state transition probabilities and priors
    for j=1:ns % states
        for i=1:ns % states
            Phi(i,j)=sum(sum(xi(:,i,j,:),4),1)/sum(sum(sum(xi(:,i,:,:),1),3),4);
        end
        pr(j)=sum(gamma(1,j,:),3)/ntr; % update prior
    end
    Phi=(Phi+u*eye(ns))/(1+u); % pseudocounts for longer state duration
    
    fprintf('.')
end;
fprintf('\n')
%##########################################################################################

% full Likelihood
L=sum(l,2);

% plot results
L(nq)
k=(ns*(ns-1)+ns*nd^2*(p+1));
AICc=2*(k+(k^2+k)/(ntr*n-k-1)-L(nq))/ntr/n
BIC=(log(n*ntr)*k-2*L(nq))/ntr/n

AC(ns,p)=AICc;
BC(ns,p)=BIC;
LH(ns,p)=L(nq);

if plotL
    figure;
    plot(L(2:nq))
    title('log[p(Y|theta)]')
    xlabel('iteration')
end;

fprintf('\n %3.2fs\n',cputime-time)

graphics
%allgr