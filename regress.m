function y=regress(G)
% MVAR estimation algorithm (Whittle 1963)
% cross-covariances G........: n x nd x np+1

[n,~,np]=size(G);
np=np-1;

A=zeros(n,n,np,np);
Ab=zeros(n,n,np,np);
A(:,:,1,1)=-D(eye(n),G,0)/Vb(eye(n),G,0);
Ab(:,:,1,1)=-Db(eye(n),G,0)/V(eye(n),G,0);

for p=1:np-1
    A(:,:,p+1,p+1)=-D(A,G,p)/Vb(Ab,G,p);
    Ab(:,:,p+1,p+1)=-Db(A,G,p)/V(A,G,p);
    for k=1:p
        A(:,:,p+1,k)=A(:,:,p,k)+A(:,:,p+1,p+1)*Ab(:,:,p,p-k+1);
        Ab(:,:,p+1,k)=Ab(:,:,p,k)+Ab(:,:,p+1,p+1)*A(:,:,p,p-k+1);
    end
end

y={squeeze(A(:,:,np,:)),squeeze(Ab(:,:,np,:)), squeeze(V(A,G,p)), squeeze(Vb(A,G,p))};
return

function M=D(A,G,p)
M=G(:,:,p+2);
if p>0
    for k=1:p
        M=M+A(:,:,p,k)*G(:,:,p-k+2);
    end
end
return

function M=V(A,G,p)
M=G(:,:,1)';
if p>0
    for k=1:p
        M=M+A(:,:,p,k)*G(:,:,k+1)';
    end
end
return

function M=Db(A,G,p)
M=G(:,:,p+2)';
if p>0
    for k=1:p
        M=M+A(:,:,p,k)*G(:,:,p-k+2)';
    end
end
return

function M=Vb(A,G,p)
M=G(:,:,1);
if p>0
    for k=1:p
        M=M+A(:,:,p,k)*G(:,:,k+1);
    end
end
return