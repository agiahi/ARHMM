function mplot(M,opt,c)
% mplot(M,[opt,c])
% 
% M.........: nd x (nd) x t
%
% opt.......: 1: hold on
%             2: rbw color scheme
%             3: flatten
%
% c.........: range ([min,max] or single value for symmetric range)

if nargin==1, opt=0;end 
if abs(opt)==1
    hold on;
else
    figure;
end
if nargin<3, c=0;end

M=squeeze(M);
if (size(M,1)-1)&&(size(M,2)-1)&&~(opt==3)
    if numel(c)>1, cmin=min(c);c=max(c);colormap(rwbmap0b(256));else
    if c==0, c=max(max(max(M))); end % greatest matrix element
    cmin=min(min(min(M))); % smallest matrix element
    if cmin<0||opt==2 % set color scheme according to lowest matrix element
        if c==0, c=max(c,abs(cmin));end
        colormap(rwbmap(256));
        cmin=-c;
    else
        cmin=0;
    end
    end
    imagesc(M(:,:)); caxis([cmin,c]);
else
    if opt==3
        plot(reshape(M,1,[]));
    else
        plot(M); 
    end
end
end