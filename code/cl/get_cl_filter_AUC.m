function cvscore=get_cl_filter_AUC(ranks,dataname,sselected,nruns,folds)
%--------------load parameter------
corefun=0;
if strcmp(dataname,'madelon')
    tdepth=2^6;
    ntrees=1900;
    nEpochs=10;
    eta=0.0005;
    shrink=0.00001;
    minibatch=145;
    mu=40;
elseif strcmp(dataname,'CAN_SMK_187')
    tdepth=2^0;
    ntrees=500;
    nEpochs=500;
    eta=0.01;
    shrink=0.001;
    minibatch=145;
    mu=280;
elseif strcmp(datanam,'GLI_85')
    tdepth=2;
    ntrees=200;
    nEpochs=800;
    eta=0.1;
    mu=101;
    minibatch=101;
    shrink=0.005;
elseif strcmp(dataname,'dexter')
    tdepth=2^2;
    ntrees=400;
    nEpochs=300;
    eta=0.000001;
    shrink=0;
    minibatch=30;
    mu=100;
elseif strcmp(dataname,'gisette')
    tdepth=2^2;
    ntrees=400;
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
%-------------------------
for r=1:nruns
    %fsa
    fprintf('run_%d\n',r);
    currentIdx=cvindex{r}.i;
    yp_test=[]; yt_test=[];
    parfor c=1:folds
        [timefsatest(c),temptest]=learner_cl_fsa(X(currentIdx==c,:),Y(currentIdx==c,:),X(currentIdx~=c,:),Y(currentIdx~=c,:),corefun,sselected,shrink,eta,nEpochs,mu,minibatch);
        yp_test=[yp_test; temptest]; yt_test=[yt_test; Y(currentIdx==c,:)];
    end
    [~,~,~,aucnow] = perfcurve(yt_test,yp_test,1);
    cvscore.fsatest.auc.(sprintf('run_%d',r))=aucnow;
    cvscore.fsatest.time.(sprintf('run_%d',r))=timefsatest;
    cvscore.fsatest.yp.(sprintf('run_%d',r))=yp_test;
    fprintf('FSA finish\n');
    
    rankmu=ranks.MI.rank.(sprintf('run_%d',r));
    rankr=ranks.Relief.rank.(sprintf('run_%d',r));
    rankt=ranks.Tscore.rank.(sprintf('run_%d',r));
    rank2=ranks.Chi2.rank.(sprintf('run_%d',r));
    rankg=ranks.Gini.rank.(sprintf('run_%d',r));
    rankf=ranks.Fisher.rank.(sprintf('run_%d',r));
    rankmm=ranks.MRMR.rank.(sprintf('run_%d',r));
    
    yp_r=[]; yp_mu=[]; yp_mm=[]; yp_t=[];  yp_f=[]; yp_2=[];  yp_g=[]; yt=[];
    parfor c=1:folds
        [timefr(c),tempr]=learner_cl_logist(X(currentIdx==c,rankr(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankr(c,1:sselected)),Y(currentIdx~=c,:));
        [timefmu(c),tempmu]=learner_cl_logist(X(currentIdx==c,rankmu(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankmu(c,1:sselected)),Y(currentIdx~=c,:));
        [timefmm(c),tempmm]=learner_cl_logist(X(currentIdx==c,rankmm(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankmm(c,1:sselected)),Y(currentIdx~=c,:));
        [timeft(c),tempt]=learner_cl_logist(X(currentIdx==c,rankt(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankt(c,1:sselected)),Y(currentIdx~=c,:));
        [timeff(c),tempf]=learner_cl_logist(X(currentIdx==c,rankf(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankf(c,1:sselected)),Y(currentIdx~=c,:));
        [timef2(c),temp2]=learner_cl_logist(X(currentIdx==c,rank2(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rank2(c,1:sselected)),Y(currentIdx~=c,:));
        [timefg(c),tempg]=learner_cl_logist(X(currentIdx==c,rankg(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankg(c,1:sselected)),Y(currentIdx~=c,:));
        yt=[yt; Y(currentIdx==c,:)];
        yp_r=[yp_r; tempr]; yp_mu=[yp_mu; tempmu]; yp_mm=[yp_mm; tempmm]; yp_t=[yp_t; tempt]; yp_f=[yp_f; tempf]; yp_2=[yp_2; temp2]; yp_g=[yp_g; tempg]; 
    end
    [~,~,~,aucr] = perfcurve(yt,yp_r,1);
    cvscore.logist.Relief.auc.(sprintf('run_%d',r))=aucr;
    cvscore.logist.Relief.time.(sprintf('run_%d',r))=timefr;
    cvscore.logist.Relief.yp.(sprintf('run_%d',r))=yp_r;
    [~,~,~,aucmu] = perfcurve(yt,yp_mu,1);
    cvscore.logist.MI.auc.(sprintf('run_%d',r))=aucmu;
    cvscore.logist.MI.time.(sprintf('run_%d',r))=timefmu;
    cvscore.logist.MI.yp.(sprintf('run_%d',r))=yp_mu;
    [~,~,~,aucmm] = perfcurve(yt,yp_mm,1);
    cvscore.logist.MRMR.auc.(sprintf('run_%d',r))=aucmm;
    cvscore.logist.MRMR.time.(sprintf('run_%d',r))=timefmm;
    cvscore.logist.MRMR.yp.(sprintf('run_%d',r))=yp_mm;
    [~,~,~,auct] = perfcurve(yt,yp_t,1);
    cvscore.logist.Tscore.auc.(sprintf('run_%d',r))=auct;
    cvscore.logist.Tscore.time.(sprintf('run_%d',r))=timeft;
    cvscore.logist.Tscore.yp.(sprintf('run_%d',r))=yp_t;
    [~,~,~,aucf] = perfcurve(yt,yp_f,1);
    cvscore.logist.Fisher.auc.(sprintf('run_%d',r))=aucf;
    cvscore.logist.Fisher.time.(sprintf('run_%d',r))=timeff;
    cvscore.logist.Fisher.yp.(sprintf('run_%d',r))=yp_f;
    [~,~,~,auc2] = perfcurve(yt,yp_2,1);
    cvscore.logist.Chi2.auc.(sprintf('run_%d',r))=auc2;
    cvscore.logist.Chi2.time.(sprintf('run_%d',r))=timef2;
    cvscore.logist.Chi2.yp.(sprintf('run_%d',r))=yp_2;
    [~,~,~,aucg] = perfcurve(yt,yp_g,1);
    cvscore.logist.Gini.auc.(sprintf('run_%d',r))=aucg;
    cvscore.logist.Gini.time.(sprintf('run_%d',r))=timefg;
    cvscore.logist.Gini.yp.(sprintf('run_%d',r))=yp_g;
    fprintf('logistic finish\n');
    yp_r=[]; yp_mu=[]; yp_mm=[]; yp_t=[];  yp_f=[]; yp_2=[];  yp_g=[]; yt=[];
    parfor c=1:folds
        [timefr(c),tempr]=learner_cl_nyb(X(currentIdx==c,rankr(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankr(c,1:sselected)),Y(currentIdx~=c,:));
        [timefmu(c),tempmu]=learner_cl_nyb(X(currentIdx==c,rankmu(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankmu(c,1:sselected)),Y(currentIdx~=c,:));
        [timefmm(c),tempmm]=learner_cl_nyb(X(currentIdx==c,rankmm(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankmm(c,1:sselected)),Y(currentIdx~=c,:));
        [timeft(c),tempt]=learner_cl_nyb(X(currentIdx==c,rankt(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankt(c,1:sselected)),Y(currentIdx~=c,:));
        [timeff(c),tempf]=learner_cl_nyb(X(currentIdx==c,rankf(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankf(c,1:sselected)),Y(currentIdx~=c,:));
        [timef2(c),temp2]=learner_cl_nyb(X(currentIdx==c,rank2(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rank2(c,1:sselected)),Y(currentIdx~=c,:));
        [timefg(c),tempg]=learner_cl_nyb(X(currentIdx==c,rankg(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankg(c,1:sselected)),Y(currentIdx~=c,:));
        yt=[yt; Y(currentIdx==c,:)];
        yp_r=[yp_r; tempr]; yp_mu=[yp_mu; tempmu]; yp_mm=[yp_mm; tempmm]; yp_t=[yp_t; tempt]; yp_f=[yp_f; tempf]; yp_2=[yp_2; temp2]; yp_g=[yp_g ;tempg]; 
    end
    [~,~,~,aucr] = perfcurve(yt,yp_r,1);
    cvscore.nby.Relief.auc.(sprintf('run_%d',r))=aucr;
    cvscore.nby.Relief.time.(sprintf('run_%d',r))=timefr;
    cvscore.nby.Relief.yp.(sprintf('run_%d',r))=yp_r;
    [~,~,~,aucmu] = perfcurve(yt,yp_mu,1);
    cvscore.nby.MI.auc.(sprintf('run_%d',r))=aucmu;
    cvscore.nby.MI.time.(sprintf('run_%d',r))=timefmu;
    cvscore.nby.MI.yp.(sprintf('run_%d',r))=yp_mu;
    [~,~,~,aucmm] = perfcurve(yt,yp_mm,1);
    cvscore.nby.MRMR.auc.(sprintf('run_%d',r))=aucmm;
    cvscore.nby.MRMR.time.(sprintf('run_%d',r))=timefmm;
    cvscore.nby.MRMR.yp.(sprintf('run_%d',r))=yp_mm;
    [~,~,~,auct] = perfcurve(yt,yp_t,1);
    cvscore.nby.Tscore.auc.(sprintf('run_%d',r))=auct;
    cvscore.nby.Tscore.time.(sprintf('run_%d',r))=timeft;
    cvscore.nby.Tscore.yp.(sprintf('run_%d',r))=yp_t;
    [~,~,~,aucf] = perfcurve(yt,yp_f,1);
    cvscore.nby.Fisher.auc.(sprintf('run_%d',r))=aucf;
    cvscore.nby.Fisher.time.(sprintf('run_%d',r))=timeff;
    cvscore.nby.Fisher.yp.(sprintf('run_%d',r))=yp_f;
    [~,~,~,auc2] = perfcurve(yt,yp_2,1);
    cvscore.nby.Chi2.auc.(sprintf('run_%d',r))=auc2;
    cvscore.nby.Chi2.time.(sprintf('run_%d',r))=timef2;
    cvscore.nby.Chi2.yp.(sprintf('run_%d',r))=yp_2;
    [~,~,~,aucg] = perfcurve(yt,yp_g,1);
    cvscore.nby.Gini.auc.(sprintf('run_%d',r))=aucg;
    cvscore.nby.Gini.time.(sprintf('run_%d',r))=timefg;
    cvscore.nby.Gini.yp.(sprintf('run_%d',r))=yp_g;
    fprintf('nby finish\n');
    yp_r=[]; yp_mu=[]; yp_mm=[]; yp_t=[];  yp_f=[]; yp_2=[];  yp_g=[]; yt=[];
    parfor c=1:folds
        [timefr(c),tempr]=learner_cl_lsvm(X(currentIdx==c,rankr(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankg(c,1:sselected)),Y(currentIdx~=c,:));
        [timefmu(c),tempmu]=learner_cl_lsvm(X(currentIdx==c,rankmu(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankmu(c,1:sselected)),Y(currentIdx~=c,:));
        [timefmm(c),tempmm]=learner_cl_lsvm(X(currentIdx==c,rankmm(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankmm(c,1:sselected)),Y(currentIdx~=c,:));
        [timeft(c),tempt]=learner_cl_lsvm(X(currentIdx==c,rankt(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankt(c,1:sselected)),Y(currentIdx~=c,:));
        [timeff(c),tempf]=learner_cl_lsvm(X(currentIdx==c,rankf(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankf(c,1:sselected)),Y(currentIdx~=c,:));
        [timef2(c),temp2]=learner_cl_lsvm(X(currentIdx==c,rank2(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rank2(c,1:sselected)),Y(currentIdx~=c,:));
        [timefg(c),tempg]=learner_cl_lsvm(X(currentIdx==c,rankg(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankg(c,1:sselected)),Y(currentIdx~=c,:));
        yt=[yt; Y(currentIdx==c,:)];
        yp_r=[yp_r; tempr]; yp_mu=[yp_mu; tempmu]; yp_mm=[yp_mm ;tempmm]; yp_t=[yp_t; tempt]; yp_f=[yp_f; tempf]; yp_2=[yp_2; temp2]; yp_g=[yp_g ;tempg]; 
    end
    [~,~,~,aucr] = perfcurve(yt,yp_r,1);
    cvscore.svm.Relief.auc.(sprintf('run_%d',r))=aucr;
    cvscore.svm.Relief.time.(sprintf('run_%d',r))=timefr;
    cvscore.svm.Relief.yp.(sprintf('run_%d',r))=yp_r;
    [~,~,~,aucmu] = perfcurve(yt,yp_mu,1);
    cvscore.svm.MI.auc.(sprintf('run_%d',r))=aucmu;
    cvscore.svm.MI.time.(sprintf('run_%d',r))=timefmu;
    cvscore.svm.MI.yp.(sprintf('run_%d',r))=yp_mu;
    [~,~,~,aucmm] = perfcurve(yt,yp_mm,1);
    cvscore.svm.MRMR.auc.(sprintf('run_%d',r))=aucmm;
    cvscore.svm.MRMR.time.(sprintf('run_%d',r))=timefmm;
    cvscore.svm.MRMR.yp.(sprintf('run_%d',r))=yp_mm;
    [~,~,~,auct] = perfcurve(yt,yp_t,1);
    cvscore.svm.Tscore.auc.(sprintf('run_%d',r))=auct;
    cvscore.svm.Tscore.time.(sprintf('run_%d',r))=timeft;
    cvscore.svm.Tscore.yp.(sprintf('run_%d',r))=yp_t;
    [~,~,~,aucf] = perfcurve(yt,yp_f,1);
    cvscore.svm.Fisher.auc.(sprintf('run_%d',r))=aucf;
    cvscore.svm.Fisher.time.(sprintf('run_%d',r))=timeff;
    cvscore.svm.Fisher.yp.(sprintf('run_%d',r))=yp_f;
    [~,~,~,auc2] = perfcurve(yt,yp_2,1);
    cvscore.svm.Chi2.auc.(sprintf('run_%d',r))=auc2;
    cvscore.svm.Chi2.time.(sprintf('run_%d',r))=timef2;
    cvscore.svm.Chi2.yp.(sprintf('run_%d',r))=yp_2;
    [~,~,~,aucg] = perfcurve(yt,yp_g,1);
    cvscore.svm.Gini.auc.(sprintf('run_%d',r))=aucg;
    cvscore.svm.Gini.time.(sprintf('run_%d',r))=timefg;
    cvscore.svm.Gini.yp.(sprintf('run_%d',r))=yp_g;
    fprintf('SVM finish\n');
    try
        yp_r=[]; yp_mu=[]; yp_mm=[]; yp_t=[];  yp_f=[]; yp_2=[];  yp_g=[]; yt=[];
        parfor c=1:folds
            fprintf('start_%d\n',c)
            drawnow
            [timefr(c),tempr]=learner_cl_boostT(X(currentIdx==c,rankr(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankr(c,1:sselected)),Y(currentIdx~=c,:),tdepth,ntrees);
            [timefmu(c),tempmu]=learner_cl_boostT(X(currentIdx==c,rankmu(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankmu(c,1:sselected)),Y(currentIdx~=c,:),tdepth,ntrees);
            [timefmm(c),tempmm]=learner_cl_boostT(X(currentIdx==c,rankmm(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankmm(c,1:sselected)),Y(currentIdx~=c,:),tdepth,ntrees);
            [timeft(c),tempt]=learner_cl_boostT(X(currentIdx==c,rankt(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankt(c,1:sselected)),Y(currentIdx~=c,:),tdepth,ntrees);
            [timeff(c),tempf]=learner_cl_boostT(X(currentIdx==c,rankf(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankf(c,1:sselected)),Y(currentIdx~=c,:),tdepth,ntrees);
            [timef2(c),temp2]=learner_cl_boostT(X(currentIdx==c,rank2(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rank2(c,1:sselected)),Y(currentIdx~=c,:),tdepth,ntrees);
            [timefg(c),tempg]=learner_cl_boostT(X(currentIdx==c,rankg(c,1:sselected)),Y(currentIdx==c,:),X(currentIdx~=c,rankg(c,1:sselected)),Y(currentIdx~=c,:),tdepth,ntrees);
            yt=[yt; Y(currentIdx==c,:)];
            yp_r=[yp_r; tempr]; yp_mu=[yp_mu; tempmu]; yp_mm=[yp_mm ;tempmm]; yp_t=[yp_t; tempt]; yp_f=[yp_f; tempf]; yp_2=[yp_2; temp2]; yp_g=[yp_g; tempg];
            fprintf('end_%d\n',c)
            drawnow
        end
    catch ME
        disp(ME.identifier)
        fprintf('run_%d\n',r)
        drawnow
        continue
    end
    [~,~,~,aucr] = perfcurve(yt,yp_r,1);
    cvscore.BT.Relief.auc.(sprintf('run_%d',r))=aucr;
    cvscore.BT.Relief.time.(sprintf('run_%d',r))=timefr;
    cvscore.BT.Relief.yp.(sprintf('run_%d',r))=yp_r;
    [~,~,~,aucmu] = perfcurve(yt,yp_mu,1);
    cvscore.BT.MI.auc.(sprintf('run_%d',r))=aucmu;
    cvscore.BT.MI.time.(sprintf('run_%d',r))=timefmu;
    cvscore.BT.MI.yp.(sprintf('run_%d',r))=yp_mu;
    [~,~,~,aucmm] = perfcurve(yt,yp_mm,1);
    cvscore.BT.MRMR.auc.(sprintf('run_%d',r))=aucmm;
    cvscore.BT.MRMR.time.(sprintf('run_%d',r))=timefmm;
    cvscore.BT.MRMR.yp.(sprintf('run_%d',r))=yp_mm;
    [~,~,~,auct] = perfcurve(yt,yp_t,1);
    cvscore.BT.Tscore.auc.(sprintf('run_%d',r))=auct;
    cvscore.BT.Tscore.time.(sprintf('run_%d',r))=timeft;
    cvscore.BT.Tscore.yp.(sprintf('run_%d',r))=yp_t;
    [~,~,~,aucf] = perfcurve(yt,yp_f,1);
    cvscore.BT.Fisher.auc.(sprintf('run_%d',r))=aucf;
    cvscore.BT.Fisher.time.(sprintf('run_%d',r))=timeff;
    cvscore.BT.Fisher.yp.(sprintf('run_%d',r))=yp_f;
    [~,~,~,auc2] = perfcurve(yt,yp_2,1);
    cvscore.BT.Chi2.auc.(sprintf('run_%d',r))=auc2;
    cvscore.BT.Chi2.time.(sprintf('run_%d',r))=timef2;
    cvscore.BT.Chi2.yp.(sprintf('run_%d',r))=yp_2;
    [~,~,~,aucg] = perfcurve(yt,yp_g,1);
    cvscore.BT.Gini.auc.(sprintf('run_%d',r))=aucg;
    cvscore.BT.Gini.time.(sprintf('run_%d',r))=timefg;
    cvscore.BT.Gini.yp.(sprintf('run_%d',r))=yp_g;
    fprintf('BT finish\n');
end
