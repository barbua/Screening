function allrank=get_cl_screen_rank(nruns,dataname,folds,verbos)
load(dataname);
Y(Y~=1)=-1;
sumx=var(X);
X=X(:,sumx~=0);
num_col=size(X,2);
%----------load cv index
cvindex=CVindex_for_screening(nruns,Y,folds,dataname);
%---------------methods
%mrmr's parameter
if num_col<5500
    param.k=num_col;
else
    param.k=5500;
end
param.pool=size(X,2);
param.type=-1;
for r=1:nruns
    fprintf('run_%d\n',r);
    currentIdx=cvindex{r}.i;
    parfor c=1:folds
        %------relief
        tic
        [rankr(c,:),scorer(c,:)]=relieff(X(currentIdx~=c,:),Y(currentIdx~=c,:),10,'method','classification');%140/dexter;10/187;410/gisette;30/GLI85
        timer(c)=toc;
        if verbos ,fprintf('%s_fold%d\n','Relief',c), end;
    end
    allrank.Relief.rank.(sprintf('run_%d',r)) = rankr;
    allrank.Relief.time.(sprintf('run_%d',r)) = timer;
    allrank.Relief.score.(sprintf('run_%d',r)) = scorer;
    parfor c=1:folds
        %-----------ttest
        tic
        [rankt(c,:),scoret(c,:)]=rankfeatures(X(currentIdx~=c,:)',Y(currentIdx~=c,:));
        timet(c)=toc;
        if verbos ,fprintf('%s_fold%d\n','T_score',c), end;
    end
    allrank.Tscore.rank.(sprintf('run_%d',r)) = rankt;
    allrank.Tscore.time.(sprintf('run_%d',r)) = timet;
    allrank.Tscore.score.(sprintf('run_%d',r)) = scoret;
    parfor c=1:folds
        %-----------mutual
        tic
        [scoremu(c,:),rankmu(c,:)]=screen_cl_mutual(X(currentIdx~=c,:),Y(currentIdx~=c,:),25);
        timemu(c)=toc;
        if verbos ,fprintf('%s_fold%d\n','MI',c),end;
    end
    allrank.MI.rank.(sprintf('run_%d',r)) = rankmu;
    allrank.MI.time.(sprintf('run_%d',r)) = timemu;
    allrank.MI.score.(sprintf('run_%d',r)) = scoremu;
        %----------fisher score
    parfor c=1:folds
        tic;
        outfish = screen_cl_fsFisher(X(currentIdx~=c,:),Y(currentIdx~=c,:));
        timef(c)=toc;
        rankf(c,:)=outfish.fList;
        scoref(c,:)=outfish.W;
        if verbos ,fprintf('%s_fold%d\n','Fisher',c), end;
    end
    allrank.Fisher.rank.(sprintf('run_%d',r)) = rankf;
    allrank.Fisher.time.(sprintf('run_%d',r)) = timef;
    allrank.Fisher.score.(sprintf('run_%d',r)) = scoref;
        %--------chi square
    parfor c=1:folds
        tic;
        [rank2(c,:),score2(c,:)]=screen_cl_chi2(myQuantileDiscretize(X(currentIdx~=c,:)),Y(currentIdx~=c,:));
        time2(c)=toc;
        if verbos ,fprintf('%s_fold%d\n','Chi2',c), end;
    end
    allrank.Chi2.rank.(sprintf('run_%d',r)) = rank2;
    allrank.Chi2.time.(sprintf('run_%d',r)) = time2;
    allrank.Chi2.score.(sprintf('run_%d',r)) = score2;
        %--------------gini idx
    parfor c=1:folds
        tic;
        [scoreg(c,:),outgini]=screen_cl_fsGini(myQuantileDiscretize(X(currentIdx~=c,:)),Y(currentIdx~=c,:));
        rankg(c,:)=outgini.fList;
        timeg(c)=toc;
        if verbos ,fprintf('%s_fold%d\n','Gini',c), end;
    end
    allrank.Gini.rank.(sprintf('run_%d',r)) = rankg;
    allrank.Gini.time.(sprintf('run_%d',r)) = timeg;
    allrank.Gini.score.(sprintf('run_%d',r)) = scoreg;
    parfor c=1:folds
        %-------------mrmr
        tic
        out= screen_cl_MRMR(X(currentIdx~=c,:),Y(currentIdx~=c,:),param);
        rankmm(c,:)=out.fList;
        timemm(c)=toc;
        if verbos ,fprintf('%s_fold%d\n','MRMR',c), end;
    end
    allrank.MRMR.rank.(sprintf('run_%d',r)) = rankmm;
    allrank.MRMR.time.(sprintf('run_%d',r)) = timemm;
end
