function matrixToGraph2(A,width,length,statecolor,b)
% convert a pdc connectivity matrix to a graph
% input: A ..............: pdc matrix
%        width ..........: arrow width
%        length .........: arrow length
%        statecolor .....: vector of colors encoding states
%        b ..............: brightness (optional) of colors (HSV system)

if nargin<4, statecolor=1;b=1;end

n = size(A,1);
d=2*pi/n/20;
if (nargin<6)
    xL = 2*pi/n * (1:n);
    vLp = [cos(xL+d);sin(xL+d)];
    vLm = [cos(xL-d);sin(xL-d)];
end
if (nargin<2)
    width = 1;
end
if (nargin<3)
    length = 15;
end

for i=1:n
    for j=1:n
        if i<j
            magAij=abs(A(i,j)); 
            magAji=abs(A(j,i)); 
            arrow(vLp(:,i),.95*vLm(:,j)+.05*vLp(:,i),7,'Length',length*magAij^.5,'BaseAngle',[],'TipAngle',[],'Width',width*magAij,'EdgeColor',colmapfunction(A(i,j),statecolor,b),'FaceColor',colmapfunction(A(i,j),statecolor,b));
            arrow(vLp(:,j),.95*vLm(:,i)+.05*vLp(:,j),7,'Length',length*magAji^.5,'BaseAngle',[],'TipAngle',[],'Width',width*magAji,'EdgeColor',colmapfunction(A(j,i),statecolor,b),'FaceColor',colmapfunction(A(j,i),statecolor,b));   
         end
    end
end