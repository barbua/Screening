function cvscore=get_cl_nofilter_AUC(nruns,dataname,folds)
%--------------load parameter------
if strcmp(dataname,'madelon')
    tdepth=2^6;
    ntrees=1900;
elseif strcmp(dataname,'CAN_SMK_187')
    tdepth=2^0;
    ntrees=500;
elseif strcmp(dataname,'GLI_85')
    tdepth=2;
    ntrees=200;
else
    tdepth=2^2;
    ntrees=400;
end
%----------------------------------------------
load(dataname);
Y(Y~=1)=-1;
sumx=var(X);
X=X(:,sumx~=0);

[n,p]=size(X);
mx=mean(X);
X=X-repmat(mx,n,1);
stdx=std(X);
X=X./repmat(stdx,n,1);
%----------load cv index
cvindex=CVindex_for_screening(nruns,Y,folds,dataname);
%-------------------------
for r=1:nruns
    fprintf('run_%d\n',r);
    currentIdx=cvindex{r}.i;
    yp=[];
    yt=[];
    parfor c=1:folds
        [times(c),temp]=learner_cl_logist(X(currentIdx==c,:),Y(currentIdx==c,:),X(currentIdx~=c,:),Y(currentIdx~=c,:));
        yp=[yp;temp];
        yt=[yt;Y(currentIdx==c,:)];
    end
    [~,~,~,aucnow] = perfcurve(yt,yp,1);
    cvscore.logist.auc.(sprintf('run_%d',r))= aucnow;
    cvscore.logist.time.(sprintf('run_%d',r))=times;
    cvscore.logist.yp.(sprintf('run_%d',r))=yp;
    fprintf('logistic finish\n');
    yp=[];
    yt=[];
    parfor c=1:folds
        [times(c),temp]=learner_cl_nyb(X(currentIdx==c,:),Y(currentIdx==c,:),X(currentIdx~=c,:),Y(currentIdx~=c,:));
        yp=[yp;temp];
        yt=[yt; Y(currentIdx==c,:)];
    end
    [~,~,~,aucnow] = perfcurve(yt,yp,1);
    cvscore.nby.auc.(sprintf('run_%d',r))=aucnow;
    cvscore.nby.time.(sprintf('run_%d',r))=times;
    cvscore.nby.yp.(sprintf('run_%d',r))=yp;
    fprintf('nby finish\n');
    yp=[];
    yt=[];
    parfor c=1:folds
        [times(c),temp]=learner_cl_lsvm(X(currentIdx==c,:),Y(currentIdx==c,:),X(currentIdx~=c,:),Y(currentIdx~=c,:));
        yp=[yp; temp];
        yt=[yt; Y(currentIdx==c,:)];
    end
    [~,~,~,aucnow] = perfcurve(yt,yp,1);
    cvscore.svm.auc.(sprintf('run_%d',r))=aucnow;
    cvscore.svm.time.(sprintf('run_%d',r))=times;
    cvscore.svm.yp.(sprintf('run_%d',r))=yp;
    fprintf('SVM finish\n');
    yp=[];
    yt=[];
    parfor c=1:folds
        [times(c),temp]=learner_cl_boostT(X(currentIdx==c,:),Y(currentIdx==c,:),X(currentIdx~=c,:),Y(currentIdx~=c,:),tdepth,ntrees);
        yp=[yp; temp];
        yt=[yt; Y(currentIdx==c,:)];
    end
    [~,~,~,aucnow] = perfcurve(yt,yp,1);
    cvscore.BT.auc.(sprintf('run_%d',r))=aucnow;
    cvscore.BT.time.(sprintf('run_%d',r))=times;
    cvscore.BT.yp.(sprintf('run_%d',r))=yp;
    fprintf('BT finish\n');
end
