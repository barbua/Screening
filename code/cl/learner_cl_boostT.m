function [ftime,py]=learner_cl_boostT(xt,yt,x,y,tdepth,ntrees)
tic;
t = templateTree('MaxNumSplits',tdepth);
mdl = fitensemble(x,y,'AdaBoostM1',ntrees,t);
ftime=toc;
[pyc,pro]=predict(mdl,xt);
[~,~,~,auc] = perfcurve(yt,pro(:,2),1);
py=pro(:,2);
pyc=pyc.*yt;
er=mean(pyc<0);