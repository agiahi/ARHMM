function indx=wordfreq(paths,L)

%
% Finds trial indices of images (default to 'unscrambled') according to
% 'wordlist' (in cell format: {'dog', 'cat', ...} in trials specified
% by 'strials' (default to 'all') shown during experiment; returns 0 if no
% occurrence found.
%
% example:  find_word(paths,'ta505',{'dog','umbrella_open'},1, [1,20:50,100])
%
%           --> finds trial numbers of scrambled images of 'dog' and
%               'umbrella_open' for subject 'ta505' in specified trials
%
%
% Aram Giahi Saravani Feb 15 2018
%
indx=zeros(1,size(L,1));
f=sprintf('%s/WordFreqList.txt',paths.root);
W=readtable(f);
W=table2cell(W(1:end,1));

for k=1:size(L,1)
    j=1;
    while (j<size(W,1))&&(~strcmpi(L{k},W{j}))
        j=j+1;
    end   
    indx(k)=j;
end

end