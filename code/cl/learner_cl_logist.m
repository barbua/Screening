function [ftime,py]=learner_cl_logist(xt,yt,x,y)
tic;
mdl = fitclinear(x,y,'Learner','logistic','Regularization','ridge');
ftime=toc;
py=predict(mdl,xt);
[~,~,~,auc] = perfcurve(yt,py,1);
py=py.*yt;
er=mean(py<0);