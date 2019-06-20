function L=get_words(paths,subject)
f=sprintf('%s%s_lfp_files/%s_common.csv',paths.vol.S1,subject,subject);
L=readtable(f);
L=table2cell(L(1:end,2));
end