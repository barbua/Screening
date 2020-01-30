nruns=7;
folds=7;
verbos=1;
dataname='madelon';
allrank=get_cl_screen_rank(nruns,dataname,folds,verbos);
results.rank=allrank;
save(sprintf('../%s_result.mat',dataname),'results');
% load(sprintf('../%s_result.mat',dataname));
% allrank=results.rank;

if strcmp(dataname , 'CLI_85')
    scalen=1.78;
elseif strcmp(dataname,'dexter')
    scalen=1.78;
elseif strcmp(dataname , 'gisette')
    scalen=1.73;
elseif strcmp(dataname,'CAN_SMK_187')
    scalen=1.78;
elseif strcmp(dataname,'madelon')
    scalen=1.25;
end

for i=1:30
    if verbos, fprintf('out_loop_%d\n',i), end;   
    filter_auc=get_cl_filter_AUC(allrank,dataname,round(((i)*4).^scalen),nruns,folds);
    results.firstLayer.(sprintf('filter_select_%d',i)).numOFfeatures=round(((i)*4).^scalen);
    results.firstLayer.(sprintf('filter_select_%d',i)).AUC=filter_auc;
    save(sprintf('../%s_result.mat',dataname),'results');
    for j=1:i
        filter_fsa_auc=get_cl_filter_fsa_AUC(allrank,dataname,round(((j)*4).^scalen),round(((i)*4).^scalen),nruns,folds);
        results.firstLayer.(sprintf('filter_select_%d',i)).secondLayer.(sprintf('filter_select_%d',j)).numOFfeatures=round(((j)*4).^scalen);
        results.firstLayer.(sprintf('filter_select_%d',i)).secondLayer.(sprintf('filter_select_%d',j)).AUC=filter_fsa_auc;    
        save(sprintf('../%s_result.mat',dataname),'results');
    end
end

nofilter_auc=get_cl_nofilter_AUC(nruns,dataname,folds);
results.benchmark=nofilter_auc;
save(sprintf('../%s_result.mat',dataname),'results');
