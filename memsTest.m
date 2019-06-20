function [ta,Sigma,ind] = memsTest(xt,order,order_test_flag)
order = order(end);

%the data needs to be dimension x number of time points x number of trials
[dim,ntr,nt] = size(xt);
%need to preprocess the data to put into format needed for the algorithm
if(nt > 1)
    xtOut = zeros(dim,nt,ntr); %need to transpose that dimension to output
    for k = 1:dim
        tmp = squeeze(xt(k,:,:));
        xtPre = zeros(ntr,nt); %reinitialize temporary variable
        
        if 1
        %ensemble pre-processing - gives each trial sample
        for tim = 1:nt
            xtPre(:,tim) = (tmp(:,tim) - mean(tmp(:,tim)))/std(tmp(:,tim));
        end
        end
        
        if 1
        %temporal preprocessing
        for trial = 1:ntr
            xtPre(trial,:) = (tmp(trial,:) - mean(tmp(trial,:)))/std(tmp(trial,:));
        end
        end
        
        %transpose before output
        xtOut(k,:,:) = transpose(xtPre);
    end
end
%xt = xtOut;
xt=permute(xt,[1 3 2]);
clear xtOut;

a00 = 0;
b00 = 0;
p0 = 0;

for r = 1:ntr
    %average over trials?
    a00 = a00 + xt(:,(2:nt),r) * xt(:,(2:nt),r)';
    b00 = b00 + xt(:,1:(nt-1),r) * xt(:,1:(nt-1),r)';
    p0 = p0 + xt(:,:,r) * xt(:,:,r)';
end

Pn.p0 = chol( (1/(ntr*nt)) * p0,'lower');
tr = Pn.p0;
 
%average over trials, take cholesky decomposition and invert
a.a00 = inv(chol( (1/ntr) * a00,'lower'));
b.b00 = inv(chol( (1/ntr) * b00,'lower'));

%initialize
AIC = zeros(1,order);
BIC = zeros(1,order);
Sigma = cell(1,order);
ta = cell(1,order);
for n = 0:(order-1)
    a.(['a' num2str(n) num2str(n+1)]) = 0;
    b.(['b' num2str(n) num2str(n+1)]) = 0;
    
    %initialize
    rn_eps = 0;
    rn_r= 0;
    rn_epsr = 0;
    %average over trials
    for r = 1:ntr
        axt = 0;
        bxt = 0;
        for k = 0:n
            axt = axt + a.(['a' num2str(n) num2str(k)]) * (xt(:,(n-k+2):(nt-k),r));
            bxt = bxt + b.(['b' num2str(n) num2str(n-k)]) * (xt(:,(n+1-k):(nt-1-k),r));
        end
        %%%%%%%%%%
        rn_eps = rn_eps + axt * axt';
        rn_r = rn_r + bxt * bxt';
        rn_epsr = rn_epsr + axt * bxt';
    end
    
    %overwrite (only need inverses of rn_eps and rn_r)
    rn_eps = chol(rn_eps,'lower');
    rn_r = chol(rn_r,'lower');
    rho_n = (rn_eps \ rn_epsr) / rn_r';
    
    Pn.(['p' num2str(n+1)]) = (chol(eye(dim) - rho_n*rho_n','lower'));
    tr = tr * Pn.(['p' num2str(n+1)]);
    Sigma{n+1} = tr*tr';
    AIC(n+1) = 2*log(det(Sigma{n+1})) + 2*(dim^2)*(n+1)/(nt*ntr);
    BIC(n+1) = 2*log(det(Sigma{n+1})) + 2*(dim^2)*(n+1)*log(nt*ntr) / (nt*ntr);
    
    
    Qn = (chol(eye(dim) - rho_n'*rho_n,'lower'));
    
    tmpa = [];
    for k = 0:(n+1)
        
        nk = [num2str(n) num2str(k)];
        nk_1 = [num2str(n) num2str(n-k+1)];
        a.(['a' num2str(n+1) num2str(k)]) = Pn.(['p' num2str(n+1)]) \ (a.(['a' nk]) - rho_n*b.(['b' nk_1]));
        b.(['b' num2str(n+1) num2str(k)]) = Qn \ (b.(['b' nk]) - rho_n'*a.(['a' nk_1]));
        
        if(k>0)
            tmpa = [tmpa, (a.(['a' num2str(n+1) '0'])) \ a.(['a' num2str(n+1) num2str(k)])];
        end
    end
    ta{n+1}=tmpa;
end

if(order_test_flag)
    figure
    plot(AIC,'o')
    hold on
    plot(BIC,'x')
    legend('AIC','BIC')
end

[minVal,ind] = min(AIC);
[minVal,ind] = min(BIC);

%parameter specified in papers (should be < 0.1)
paramOrd = dim * (ind+1) / (nt*ntr);
if(paramOrd > 0.1)
    disp('parameter failed test with aic order')
    disp(' ')
    paramOrd
    %    ind = floor(.1*nt*ntr / (dim^2));
end
Sigma = Sigma{ind};
return