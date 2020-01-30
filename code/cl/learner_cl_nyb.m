function [mtime,py]=learner_cl_nyb(xt,yt,x,y)
tic;
mdl=fitcnb(x,y,'distributionnames','kernel');
mtime=toc;
[pyc,pro]=predict(mdl,xt);
[~,~,~,auc] = perfcurve(yt,pro(:,2),1);
py=pro(:,2);
pyc=pyc.*yt;
miner=mean(pyc<0);
