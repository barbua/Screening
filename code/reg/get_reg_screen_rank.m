function allrank=get_reg_screen_rank(nruns,dataname,folds,verbos)
%--------------------load data
load(dataname);
sumx=var(X);
X=X(:,sumx~=0);
[n,p]=size(X);
%----------load cv index
cvindex=CVindex_for_screening(nruns,Y,folds,dataname);
%---------------methods
for r=1:nruns
    fprintf('run_%d\n',r);
	drawnow
    currentIdx=cvindex{r}.i;
    parfor c=1:folds
        tic
        [rankr(c,:),scorer(c,:)]=relieff(X(currentIdx~=c,:),Y(currentIdx~=c,:),10);%37;130/BMI;20/tumor;60/piazoall;5/2006;10/meat;75/ct;80/vggface;5/drive,5/indoor
        timer(c)=toc;
        if verbos ,fprintf('%s_fold%d\n','Relief',c),drawnow, end;
    end
    allrank.Relief.rank.(sprintf('run_%d',r)) = rankr;
    allrank.Relief.time.(sprintf('run_%d',r)) = timer;
    allrank.Relief.score.(sprintf('run_%d',r)) = scorer;
	fprintf('Relief finish\n');
    drawnow
    
    parfor c=1:folds
        tic
        [scoremu(c,:),rankmu(c,:)]=screen_reg_mutual(X(currentIdx~=c,:),Y(currentIdx~=c,:),10);%piazo 20;bmi 15;tumor 15;20face;15 16;6 2006;25 meat;2/ct;5/drive;10/indoor
        timemu(c)=toc;
        if verbos ,fprintf('%s_fold%d\n','MI',c),drawnow, end;
    end
    allrank.MI.rank.(sprintf('run_%d',r)) = rankmu;
    allrank.MI.time.(sprintf('run_%d',r)) = timemu;
    allrank.MI.score.(sprintf('run_%d',r)) = scoremu;
    fprintf('MI finish\n');
    drawnow
	
    parfor c=1:folds
        tic
        [scorec(c,:),rankc(c,:)]=screen_reg_corrB(X(currentIdx~=c,:),Y(currentIdx~=c,:));
        timec(c)=toc;
        if verbos ,fprintf('%s_fold%d\n','corr',c), drawnow,end;
    end
    allrank.corr.rank.(sprintf('run_%d',r)) = rankc;
    allrank.corr.time.(sprintf('run_%d',r)) = timec;
    allrank.corr.score.(sprintf('run_%d',r)) = scorec;
	fprintf('corr finish\n');
    drawnow
end