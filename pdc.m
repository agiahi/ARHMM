function y=pdc(A,f,dt)
% computes partial directed coherence

F=Af(A,f,dt);
%F=inv(Af(A,f,dt)); %DTF
[n,~]=size(F);
y=zeros(n);
for a=1:n
    for b=1:n
        y(a,b)=abs(F(a,b))/sqrt(F(:,b)'*F(:,b));
        %y(a,b)=abs(F(a,b))^2/(F(a,:)*F(a,:)'); % DTF
    end
end
y=y.^1;
end