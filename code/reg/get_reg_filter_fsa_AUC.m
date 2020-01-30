function cvscore=get_reg_filter_fsa_AUC(ranks,dataname,kstart,sselected,nruns,folds,verbos)
%--------------load parameter------
corefun=4;
if strcmp(dataname,'BMI')
    nEpochs=150;
    eta=0.00001;
    mu=800;
    minibatch=285;
    shrink=0.0001;
elseif strcmp(dataname,'Tumor')
    nEpochs=50;
    eta=0.0003;
    mu=30;
    minibatch=250;
    shrink=0.001;
elseif strcmp(dataname,'CoEPrA20063')
    nEpochs=101;
    eta=0.0001;
    mu=651;
    minibatch=15;
    shrink=0.9;
elseif strcmp(dataname,'indoorloc')
    nEpochs=250;
    eta=0.00001;
    mu=200;
    minibatch=30;
    shrink=0.0001;
elseif strcmp(dataname,'wikidataxyvggface')
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

    rankmu=ranks.MI.rank.(sprintf('run_%d',r));
    rankr=ranks.Relief.rank.(sprintf('run_%d',r));
    rankc=ranks.corr.rank.(sprintf('run_%d',r));

    yp_r=[]; yp_mu=[]; yp_c=[]; yt=[];
    parfor c=1:folds 
        [timefr(c),tempr]=learner_reg_fsa(X(currentIdx==c,rankr(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankr(c,1:sselected)),Y(currentIdx~=c,:),corefun,kstart,shrink,eta,nEpochs,mu,minibatch);
        [timefmu(c),tempmu]=learner_reg_fsa(X(currentIdx==c,rankmu(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankmu(c,1:sselected)),Y(currentIdx~=c,:),corefun,kstart,shrink,eta,nEpochs,mu,minibatch);
        [timefc(c),tempc]=learner_reg_fsa(X(currentIdx==c,rankc(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankc(c,1:sselected)),Y(currentIdx~=c,:),corefun,kstart,shrink,eta,nEpochs,mu,minibatch);
        yt=[yt; Y(currentIdx==c,:)];
        yp_r=[yp_r; tempr]; yp_mu=[yp_mu; tempmu]; yp_c=[yp_c; tempc];
        if verbos ,fprintf('%s_fold%d\n','innerFSA',c),drawnow, end;
    end
    dif=yt-yp_r; er=sum(dif.^2); rsqr=1-er/ssy; rsqr(rsqr<0)=0;
    cvscore.FSA.Relief.rsq.(sprintf('run_%d',r))=rsqr;
    cvscore.FSA.Relief.time.(sprintf('run_%d',r))=timefr;
    cvscore.FSA.Relief.yp.(sprintf('run_%d',r))=yp_r;
    dif=yt-yp_mu; er=sum(dif.^2); rsqmu=1-er/ssy; rsqmu(rsqmu<0)=0;
    cvscore.FSA.MI.rsq.(sprintf('run_%d',r))=rsqmu;
    cvscore.FSA.MI.time.(sprintf('run_%d',r))=timefmu;
    cvscore.FSA.MI.yp.(sprintf('run_%d',r))=yp_mu;
    dif=yt-yp_c; er=sum(dif.^2); rsqc=1-er/ssy; rsqc(rsqc<0)=0;
    cvscore.FSA.corr.rsq.(sprintf('run_%d',r))=rsqc;
    cvscore.FSA.corr.time.(sprintf('run_%d',r))=timefc;
    cvscore.FSA.corr.yp.(sprintf('run_%d',r))=yp_c;
end