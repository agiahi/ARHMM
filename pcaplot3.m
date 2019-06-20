function [explained,pc]=pcaplot3(vecStf,a,c)
% PCA in 3 dimensions
% vecStf...........: dim x time
[B,~]=size(vecStf);
[principal,~,~,~,explained]=pca(vecStf); 
crd3=vecStf*principal(:,1:3);
if a==-1, a=linspace(10,100,B); end;
if c==-1, c=linspace(1,10,B); end;
figure; hold on; 
scatter3(crd3(:,1),crd3(:,2),crd3(:,3),a,c,'filled');
size(crd3)
if nargout==2
    pc=principal(:,1:3);
end
end