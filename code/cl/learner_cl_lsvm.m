function [ftime,py]=learner_cl_lsvm(xt,yt,x,y)
tic;
mdl = fitclinear(x,y,'Learner','svm');
ftime=toc;
[pyc,pro]=predict(mdl,xt);
[~,~,~,auc] = perfcurve(yt,pro(:,2),1);
py=pro(:,2);
pyc=pyc.*yt;
er=mean(pyc<0);