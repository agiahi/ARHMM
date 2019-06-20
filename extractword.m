function word=extractword(L)
L=char(L);
for k=1:size(L,1)
    a=find(L(k,:)=='_')+1;
    e=find(L(k,:)=='.')-1;
    if size(a,2)>1
        word{k}=L(k,a(1):a(2)-2);
    else
        word{k}=L(k,a:e);
    end
    switch word{k}
        case 'brownpelican'
            word{k}='pelican';
        case 'well'
            word{k}='wells';
        case 'boot'
            word{k}='boots';
    end
end
word=word';
end