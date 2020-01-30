function cvidx=CVindex_for_screening(runs,y,folds,dataname)
sname=sprintf('../cvindex/%dcv_%s_%druns.txt',folds,dataname,runs);
f=fopen(sname);
if (f<0)
    for r=1:runs
        indices = crossvalind('Kfold',y,folds);
        itable(r,:)=indices;
        cvidx{r}.i=indices;
    end
    dlmwrite(sname,itable);
else
    fclose(f);
    ptable=dlmread(sname);
    for r=1:runs
        cvidx{r}.i=ptable(r,:);
    end
end
