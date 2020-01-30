function cvscore=get_reg_nofilter_AUC(nruns,dataname,folds,verbos)
%--------------load parameter------
if strcmp(dataname,'BMI')
    tdepth=2^0;
    ntrees=100;
    shrink=0.0001;
elseif strcmp(dataname,'Tumor')
    tdepth=2^0;
    ntrees=10;
    shrink=0.001;
elseif strcmp(dataname,'CoEPrA20063')
    tdepth=2^0;
    ntrees=10;
    shrink=0.9;
elseif strcmp(dataname,'indoorloc')
    tdepth=2^3;
    ntrees=500;
    shrink=0.0001;
elseif strcmp(dataname,'wikidataxyvggface')
    tdepth=2^0;
    ntrees=400;
    shrink=0.001;
end
%------------load data and standardize------
load(dataname);
sumx=var(X);
X=X(:,sumx~=0);

[n,p]=size(X);
mx=mean(X);
X=X-repmat(mx,n,1);
stdx=std(X);
X=X./repmat(stdx,n,1);
Y=Y-mean(Y);
%----------load cv index
cvindex=CVindex_for_screening(nruns,Y,folds,dataname);
%---------loop and method
ssy=sum((Y-mean(Y)).^2);
for r=1:nruns
    fprintf('run_%d\n',r);
    drawnow
    currentIdx=cvindex{r}.i;

    yp=[];
    yt=[];
    for c=1:folds
        [times(c),temp]=learner_reg_ridge(X(currentIdx==c,:),Y(currentIdx==c,:),X(currentIdx~=c,:),Y(currentIdx~=c,:),shrink);
        yp=[yp;temp];
        yt=[yt;Y(currentIdx==c,:)];
        if verbos ,fprintf('%s_fold%d\n','ridge',c),drawnow, end;
    end
    dif=yt-yp; er=sum(dif.^2); rsq=1-er/ssy; rsq(rsq<0)=0;
    cvscore.ridge.rsq.(sprintf('run_%d',r))= rsq;
    cvscore.ridge.time.(sprintf('run_%d',r))=times;
    cvscore.ridge.yp.(sprintf('run_%d',r))=yp;
    fprintf('ridge finish\n');
    drawnow

    yp=[];
    yt=[];
    parfor c=1:folds
        [times(c),temp]=learner_reg_boostT(X(currentIdx==c,:),Y(currentIdx==c,:),X(currentIdx~=c,:),Y(currentIdx~=c,:),tdepth,ntrees);
        yp=[yp;temp];
        yt=[yt;Y(currentIdx==c,:)];
        if verbos ,fprintf('%s_fold%d\n','BT',c),drawnow, end;
    end
    dif=yt-yp; er=sum(dif.^2); rsq=1-er/ssy; rsq(rsq<0)=0;
    cvscore.BT.rsq.(sprintf('run_%d',r))= rsq;
    cvscore.BT.time.(sprintf('run_%d',r))=times;
    cvscore.BT.yp.(sprintf('run_%d',r))=yp;
    fprintf('BT finish\n');
    drawnow
end

