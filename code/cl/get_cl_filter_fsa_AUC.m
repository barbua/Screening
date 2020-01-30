function cvscore=get_cl_filter_fsa_AUC(ranks,dataname,kstart,sselected,nruns,folds)
%"ranks" is the index of features sorted by importance accroding to their weights.
%--------------load parameter------
corefun=0;
if strcmp(dataname,'madelon')
    nEpochs=10;
    eta=0.0005;
    shrink=0.00001;
    minibatch=145;
    mu=40;
elseif strcmp(dataname,'CAN_SMK_187')
    nEpochs=500;
    eta=0.01;
    shrink=0.001;
    minibatch=145;
    mu=280;
elseif strcmp(datanam,'GLI_85')
    nEpochs=800;
    eta=0.1;
    mu=101;
    minibatch=101;
    shrink=0.005;
elseif strcmp(dataname,'dexter')
    nEpochs=300;
    eta=0.000001;
    shrink=0;
    minibatch=30;
    mu=100;
elseif strcmp(dataname,'gisette')
    nEpochs=60;
    eta=0.0001;
    shrink=0;
    minibatch=20;
    mu=600;
end

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
%------------
for r=1:nruns
    fprintf('run_%d\n',r);
    currentIdx=cvindex{r}.i;
    rankmu=ranks.MI.rank.(sprintf('run_%d',r));
    rankr=ranks.Relief.rank.(sprintf('run_%d',r));
    rankt=ranks.Tscore.rank.(sprintf('run_%d',r));
    rank2=ranks.Chi2.rank.(sprintf('run_%d',r));
    rankg=ranks.Gini.rank.(sprintf('run_%d',r));
    rankf=ranks.Fisher.rank.(sprintf('run_%d',r));
    rankmm=ranks.MRMR.rank.(sprintf('run_%d',r));
    
    yp_r=[]; yp_mu=[]; yp_mm=[]; yp_t=[];  yp_f=[]; yp_2=[];  yp_g=[]; yt=[];
    parfor c=1:folds
        [timefr(c),tempr]=learner_cl_fsa(X(currentIdx==c,rankr(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankr(c,1:sselected)),Y(currentIdx~=c,:),corefun,kstart,shrink,eta,nEpochs,mu,minibatch);
        [timefmu(c),tempmu]=learner_cl_fsa(X(currentIdx==c,rankr(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankr(c,1:sselected)),Y(currentIdx~=c,:),corefun,kstart,shrink,eta,nEpochs,mu,minibatch);
        [timefmm(c),tempmm]=learner_cl_fsa(X(currentIdx==c,rankr(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankr(c,1:sselected)),Y(currentIdx~=c,:),corefun,kstart,shrink,eta,nEpochs,mu,minibatch);
        [timeft(c),tempt]=learner_cl_fsa(X(currentIdx==c,rankr(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankr(c,1:sselected)),Y(currentIdx~=c,:),corefun,kstart,shrink,eta,nEpochs,mu,minibatch);
        [timeff(c),tempf]=learner_cl_fsa(X(currentIdx==c,rankr(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankr(c,1:sselected)),Y(currentIdx~=c,:),corefun,kstart,shrink,eta,nEpochs,mu,minibatch);
        [timef2(c),temp2]=learner_cl_fsa(X(currentIdx==c,rankr(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankr(c,1:sselected)),Y(currentIdx~=c,:),corefun,kstart,shrink,eta,nEpochs,mu,minibatch);
        [timefg(c),tempg]=learner_cl_fsa(X(currentIdx==c,rankr(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankr(c,1:sselected)),Y(currentIdx~=c,:),corefun,kstart,shrink,eta,nEpochs,mu,minibatch);
        yt=[yt; Y(currentIdx==c,:)];
        yp_r=[yp_r; tempr]; yp_mu=[yp_mu; tempmu]; yp_mm=[yp_mm; tempmm]; yp_t=[yp_t; tempt]; yp_f=[yp_f; tempf]; yp_2=[yp_2; temp2]; yp_g=[yp_g; tempg]; 
    end
    [~,~,~,aucr] = perfcurve(yt,yp_r,1);
    cvscore.FSA.Relief.auc.(sprintf('run_%d',r))=aucr;
    cvscore.FSA.Relief.time.(sprintf('run_%d',r))=timefr;
    cvscore.FSA.Relief.yp.(sprintf('run_%d',r))=yp_r;
    [~,~,~,aucmu] = perfcurve(yt,yp_mu,1);
    cvscore.FSA.MI.auc.(sprintf('run_%d',r))=aucmu;
    cvscore.FSA.MI.time.(sprintf('run_%d',r))=timefmu;
    cvscore.FSA.MI.yp.(sprintf('run_%d',r))=yp_mu;
    [~,~,~,aucmm] = perfcurve(yt,yp_mm,1);
    cvscore.FSA.MRMR.auc.(sprintf('run_%d',r))=aucmm;
    cvscore.FSA.MRMR.time.(sprintf('run_%d',r))=timefmm;
    cvscore.FSA.MRMR.yp.(sprintf('run_%d',r))=yp_mm;
    [~,~,~,auct] = perfcurve(yt,yp_t,1);
    cvscore.FSA.Tscore.auc.(sprintf('run_%d',r))=auct;
    cvscore.FSA.Tscore.time.(sprintf('run_%d',r))=timeft;
    cvscore.FSA.Tscore.yp.(sprintf('run_%d',r))=yp_t;
    [~,~,~,aucf] = perfcurve(yt,yp_f,1);
    cvscore.FSA.Fisher.auc.(sprintf('run_%d',r))=aucf;
    cvscore.FSA.Fisher.time.(sprintf('run_%d',r))=timeff;
    cvscore.FSA.Fisher.yp.(sprintf('run_%d',r))=yp_f;
    [~,~,~,auc2] = perfcurve(yt,yp_2,1);
    cvscore.FSA.Chi2.auc.(sprintf('run_%d',r))=auc2;
    cvscore.FSA.Chi2.time.(sprintf('run_%d',r))=timef2;
    cvscore.FSA.Chi2.yp.(sprintf('run_%d',r))=yp_2;
    [~,~,~,aucg] = perfcurve(yt,yp_g,1);
    cvscore.FSA.Gini.auc.(sprintf('run_%d',r))=aucg;
    cvscore.FSA.Gini.time.(sprintf('run_%d',r))=timefg;
    cvscore.FSA.Gini.yp.(sprintf('run_%d',r))=yp_g;
end