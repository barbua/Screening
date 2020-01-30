function cvscore=get_reg_filter_AUC(ranks,dataname,sselected,nruns,folds,verbos)
%--------------load parameter------ 
corefun=4;
if strcmp(dataname,'BMI')
    tdepth=2^0;
    ntrees=100;
    nEpochs=150;
    eta=0.00001;
    mu=800;
    minibatch=285;
    shrink=0.0001;
elseif strcmp(dataname,'Tumor')
    tdepth=2^0;
    ntrees=10;
    nEpochs=50;
    eta=0.0003;
    mu=30;
    minibatch=250;
    shrink=0.001;
elseif strcmp(dataname,'CoEPrA20063')
    tdepth=2^0;
    ntrees=10;
    nEpochs=101;
    eta=0.0001;
    mu=651;
    minibatch=15;
    shrink=0.9;
elseif strcmp(dataname,'indoorloc')
    tdepth=2^3;
    ntrees=500;
    nEpochs=250;
    eta=0.00001;
    mu=200;
    minibatch=30;
    shrink=0.0001;
elseif strcmp(dataname,'wikidataxyvggface')
    tdepth=2^0;
    ntrees=400;
    nEpochs=451;
    eta=0.001;
    mu=251;
    minibatch=151;
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
    yp_test=[]; yt_test=[];
    parfor c=1:folds
        [timefsatest(c),temptest]=learner_reg_fsa(X(currentIdx==c,:),Y(currentIdx==c,:),X(currentIdx~=c,:),Y(currentIdx~=c,:),corefun,sselected,shrink,eta,nEpochs,mu,minibatch);
        yp_test=[yp_test; temptest]; yt_test=[yt_test; Y(currentIdx==c,:)];
        if verbos ,fprintf('%s_fold%d\n','FSA',c),drawnow, end;
    end   
    dif=yt_test-yp_test; er=sum(dif.^2); rsq=1-er/ssy; rsq(rsq<0)=0;
    cvscore.fsatest.rsq.(sprintf('run_%d',r))=rsq;
    cvscore.fsatest.time.(sprintf('run_%d',r))=timefsatest;
    cvscore.fsatest.yp.(sprintf('run_%d',r))=yp_test;
    fprintf('FSA finish\n');
    drawnow    
    
    rankmu=ranks.MI.rank.(sprintf('run_%d',r));
    rankr=ranks.Relief.rank.(sprintf('run_%d',r));
    rankc=ranks.corr.rank.(sprintf('run_%d',r));
    
    yp_r=[]; yp_mu=[]; yp_c=[]; yt=[];
    parfor c=1:folds
        [timemu(c),tempmu]=learner_reg_ridge(X(currentIdx==c,rankmu(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankmu(c,1:sselected)),Y(currentIdx~=c,:),shrink);
        [timec(c),tempc]=learner_reg_ridge(X(currentIdx==c,rankc(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankc(c,1:sselected)),Y(currentIdx~=c,:),shrink);
        [timer(c),tempr]=learner_reg_ridge(X(currentIdx==c,rankr(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankr(c,1:sselected)),Y(currentIdx~=c,:),shrink);
        yt=[yt; Y(currentIdx==c,:)];
        yp_r=[yp_r; tempr]; yp_mu=[yp_mu; tempmu]; yp_c=[yp_c; tempc];
        if verbos ,fprintf('%s_fold%d\n','ridge',c),drawnow, end;
    end
    dif=yt-yp_r; er=sum(dif.^2); rsqr=1-er/ssy; rsqr(rsqr<0)=0;
    cvscore.ridge.Relief.rsq.(sprintf('run_%d',r))=rsqr;
    cvscore.ridge.Relief.time.(sprintf('run_%d',r))=timer;
    cvscore.ridge.Relief.yp.(sprintf('run_%d',r))=yp_r;
    dif=yt-yp_mu; er=sum(dif.^2); rsqmu=1-er/ssy; rsqmu(rsqmu<0)=0;
    cvscore.ridge.MI.rsq.(sprintf('run_%d',r))=rsqmu;
    cvscore.ridge.MI.time.(sprintf('run_%d',r))=timemu;
    cvscore.ridge.MI.yp.(sprintf('run_%d',r))=yp_mu;
    dif=yt-yp_c; er=sum(dif.^2); rsqc=1-er/ssy; rsqc(rsqc<0)=0;
    cvscore.ridge.corr.rsq.(sprintf('run_%d',r))=rsqc;
    cvscore.ridge.corr.time.(sprintf('run_%d',r))=timec;
    cvscore.ridge.corr.yp.(sprintf('run_%d',r))=yp_c;
    fprintf('ridge finish\n');
    drawnow

    yp_r=[]; yp_mu=[]; yp_c=[]; yt=[];
    parfor c=1:folds
        [timemu(c),tempmu]=learner_reg_boostT(X(currentIdx==c,rankmu(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankmu(c,1:sselected)),Y(currentIdx~=c,:),tdepth,ntrees);
        [timec(c),tempc]=learner_reg_boostT(X(currentIdx==c,rankc(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankc(c,1:sselected)),Y(currentIdx~=c,:),tdepth,ntrees);
        [timer(c),tempr]=learner_reg_boostT(X(currentIdx==c,rankr(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankr(c,1:sselected)),Y(currentIdx~=c,:),tdepth,ntrees);
        yt=[yt; Y(currentIdx==c,:)];
        yp_r=[yp_r; tempr]; yp_mu=[yp_mu; tempmu]; yp_c=[yp_c; tempc];
        if verbos ,fprintf('%s_fold%d\n','BT',c),drawnow, end;
    end
    dif=yt-yp_r; er=sum(dif.^2); rsqr=1-er/ssy; rsqr(rsqr<0)=0;
    cvscore.BT.Relief.rsq.(sprintf('run_%d',r))=rsqr;
    cvscore.BT.Relief.time.(sprintf('run_%d',r))=timer;
    cvscore.BT.Relief.yp.(sprintf('run_%d',r))=yp_r;
    dif=yt-yp_mu; er=sum(dif.^2); rsqmu=1-er/ssy; rsqmu(rsqmu<0)=0;
    cvscore.BT.MI.rsq.(sprintf('run_%d',r))=rsqmu;
    cvscore.BT.MI.time.(sprintf('run_%d',r))=timemu;
    cvscore.BT.MI.yp.(sprintf('run_%d',r))=yp_mu;
    dif=yt-yp_c; er=sum(dif.^2); rsqc=1-er/ssy; rsqc(rsqc<0)=0;
    cvscore.BT.corr.rsq.(sprintf('run_%d',r))=rsqc;
    cvscore.BT.corr.time.(sprintf('run_%d',r))=timec;
    cvscore.BT.corr.yp.(sprintf('run_%d',r))=yp_c;
    fprintf('BT finish\n');
    drawnow
end