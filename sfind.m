function score=sfind(list,item)

%
% finds items in a list structure and returns row index
%
% Input  
%         list ......: cell array
%         item ......: query (in structure format: {'.','.',...})
%
% Output 
%         score .....: row index of item match, returns 0 if item not in list
%
% Aram Giahi Saravani Feb 15 2018
%

score=zeros(size(item,2),1);
for j=1:size(item,2)
    k=0;
    while k<size(list,1)&&score(j)<size(item{j},2)
        k=k+1;
        if size(list{k},2)==size(item{j},2)
            score(j)=sum(list{k}==item{j});
        end
    end
    if score<size(item{j},2)
        score(j)=0;
    else
        score(j)=k;
    end
end
end