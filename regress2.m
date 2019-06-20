function [c,q]=regress2(y,np)
[nd,nt,n]=size(y);
yy=reshape(y,nd,nt*n);

% initial conditions
A=zeros(nd,nd,1,2); B=A;
Sy=yy*yy';
PP=chol((Sy-yy(:,1:nt)*yy(:,1:nt)')/(nt*n),'lower');
A(:,:,1,1)=eye(nd)/PP;
B(:,:,1,1)=eye(nd)/chol((Sy-yy(:,end-nt+1:end)*yy(:,end-nt+1:end)')/(nt*n),'lower');

% order update
yp=permute(y,[1 3 2]);
for p=1:np+1
    
    if 1
        % compute the normalized innvovations
        AA=reshape(A(:,:,p,1:p),nd,nd*p); BB=reshape(B(:,:,p,p:-1:1),nd,nd*p);
        e=zeros(nd,nt,n); r=zeros(nd,nt,n); 
        for t=p:n
            yr=reshape(yp(:,t:-1:t-p+1,:),nd*p,nt);
            e(:,:,t)=AA*yr;
            r(:,:,t)=BB*yr;           
        end
        ee=reshape(e(:,:,p+1:n),nd,[]);
        rr=reshape(r(:,:,p+1:n),nd,[]);
        er=reshape(r(:,:,p:n-1),nd,[]);
        Rne=ee*ee';
        Rnr=rr*rr';
        Rner=ee*er';
    else
        Rne=zeros(nd);Rnr=Rne;Rner=Rne;
        for r = 1:nt
            axt = 0;
            bxt = 0;
            for k = 1:p
                axt = axt + A(:,:,p,k) * squeeze(y(:,r,(p-k+2):(n-k)));
                bxt = bxt + B(:,:,p,p-k+1) * squeeze(y(:,r,(p+1-k):(n-1-k)));
            end
            Rne = Rne + axt * axt';
            Rnr = Rnr + bxt * bxt';
            Rner = Rner + axt * bxt';
        end
    end
    
    rho=chol(Rne,'lower')\Rner/chol(Rnr,'lower')';
    
    % update the normalized predictors
    P=chol(eye(nd)-rho*rho','lower');
    PP=PP*P;
    P=eye(nd)/P;
    Q=eye(nd)/chol((eye(nd)-rho'*rho),'lower');
    
    for k=1:p+1
        A(:,:,p+1,k)=P*(A(:,:,p,k)-rho*B(:,:,p,p-k+2));
        B(:,:,p+1,k)=Q*(B(:,:,p,k)-rho'*A(:,:,p,p-k+2));
    end
    A(:,:,p+1,p+2)=zeros(nd); B(:,:,p+1,p+2)=zeros(nd);
end

c=zeros(nd,nd,p-1);
for k=2:p
    c(:,:,k-1)=squeeze(A(:,:,p,1)\A(:,:,p,k));
end
q=PP*PP';