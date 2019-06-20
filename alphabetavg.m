function [ gamma, xi, L, vi, vmax] = alphabetavg(pr,Phi,b)
% the states are static or LDSs (time varying m)
% state transition matrix..Phi: ns x ns
% m: state means.............m: nd x ns x n 
% state covariance...........Q: nd x nd x ns 
% prior state probability...pr: ns
% b(j,t)~P(Y(t)|s=j).........b: ns x n

vit=1-(nargout<4);

% initialization
[ns, n]=size(b); % determine time dimension
xi=zeros(n-1,ns,ns); % <S(t)S(t+1)>, needed for parameter estimation

% alpha-beta-recursion
alpha=zeros(n,ns);
alpha(1,:)=pr; % prior on initial state
beta=zeros(n,ns);
p=zeros(n,1); % auxiliary variable (for numerics): p(1)*p(2)*...p(t)=P(Y)

% forward (alpha) recursion: P(X(1..t),s(t)|M)
if vit % viterbi recursion
viterbi=zeros(n,ns);
vistate=zeros(n,ns);
vi=zeros(n,1);
vmax=zeros(n,1);
% forward recursion
for t=1:n % time
    for j=1:ns % states
        if t>1
            alphai=b(j,t)*Phi(:,j);
            alpha(t,j)=alpha(t-1,:)*alphai; % inner product: sum of all paths toward state j
            vaux=viterbi(t-1,:)'+log(alphai); % all paths towards state j
            viterbi(t,j)=max(vaux); % max path probability towards state j
            vistate(t,j)=min(find(viterbi(t,j)==vaux)); % find maximizing state i towards j
        else
            alpha(1,j)=alpha(1,j)*b(j,1);
            viterbi(1,j)=log(pr(j)*b(j,1));
        end
    end
    p(t)=sum(alpha(t,:));
    alpha(t,:)=alpha(t,:)/p(t);
end
% backtrace
vi(n)=find(max(viterbi(n,:)')==viterbi(n,:)');
vmax(n)=viterbi(n,vi(n));
for t=n:-1:2
    vi(t-1)=vistate(t,vi(t));
    vmax(t-1)=viterbi(t-1,vi(t-1));
end
else
    for t=1:n % time
    for j=1:ns % states
        if t>1
            alpha(t,j)=alpha(t-1,:)*Phi(:,j); % marginal over previous states
        end
        alpha(t,j)=alpha(t,j)*b(j,t);
    end
    p(t)=sum(alpha(t,:));
    alpha(t,:)=alpha(t,:)/p(t);
    end
end

% backward (beta) recursion: P(X(t+1..T)|s(t),M)
beta(n,:)=alpha(n,:);
for t=n-1:-1:1 % time
    for j=1:ns % states       
        for k=1:ns % states
            aq=alpha(t,:)*Phi(:,k);
            xi(t,j,k)=alpha(t,j)*beta(t+1,k)*Phi(j,k)/aq;
            beta(t,j)=beta(t,j)+xi(t,j,k);
        end
    end
end
gamma=beta;

L=sum(log(p)); % log likelihood P(Y|model): log(p1*p2*p3*...*pT)=sum(log(p(1..T)))
end