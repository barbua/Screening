ncores=feature('numcores');
drawnow
parpool(ncores-1);
nruns=7;
folds=7;
verbos=1;
dataname='indoorloc';
fprintf('-------dataset: %s, used workers: %d--------\n',dataname,ncores-1)
drawnow

allrank=get_reg_screen_rank(nruns,dataname,folds,verbos);
results.rank=allrank;
save(sprintf('../%s_result.mat',dataname),'results');
% load(sprintf('../%s_result.mat',dataname));
% allrank=results.rank;

outter_loop=30;
if strcmp(dataname , 'BMI')
    scalen=1.825;
elseif strcmp(dataname,'Tumor')
    scalen=1.825;
elseif strcmp(dataname , 'CoEPrA20063')
    scalen=1.68;
elseif strcmp(dataname,'wikidataxyvggface')
    scalen=1.65;
    outter_loop=23;
elseif strcmp(dataname,'indoorloc')
    scalen=1.25;
end


for i=1:outter_loop
    if verbos, fprintf('out_loop_%d!!!!!!!!!\n',i),drawnow, end;
    filter_rsq=get_reg_filter_AUC(allrank,dataname,round(((i)*4).^scalen),nruns,folds,verbos);
    results.firstLayer.(sprintf('filter_select_%d',i)).numOFfeatures=round(((i)*4).^scalen);
	results.firstLayer.(sprintf('filter_select_%d',i)).rsq=filter_rsq;
	save(sprintf('../%s_result.mat',dataname),'results');
    for j=1:i
        if verbos ,fprintf('inner_loop_%d!!!!\n',j),drawnow, end;
        filter_fsa_rsq=get_reg_filter_fsa_AUC(allrank,dataname,round(((j)*4).^scalen),round(((i)*4).^scalen),nruns,folds,verbos);
        results.firstLayer.(sprintf('filter_select_%d',i)).secondLayer.(sprintf('filter_select_%d',j)).numOFfeatures=round(((j)*4).^scalen);
        results.firstLayer.(sprintf('filter_select_%d',i)).secondLayer.(sprintf('filter_select_%d',j)).rsq=filter_fsa_rsq;    
        save(sprintf('../%s_result.mat',dataname),'results');
    end
end

fprintf('benchmark-------------\n');
drawnow
nofilter_rsq=get_reg_nofilter_AUC(nruns,dataname,folds,verbos);
results.benchmark=nofilter_rsq;
save(sprintf('../%s_result.mat',dataname),'results');
fprintf('------all done-------\n');
drawnow
