function [ftime,py]=learner_reg_boostT(xvalid,yvalid,xtrain,ytrain,tdepth,ntrees)
ssy=sum((yvalid-mean(yvalid)).^2);

tic
t = templateTree('MaxNumSplits',tdepth);
mdl=fitensemble(xtrain,ytrain,'LSBoost',ntrees,t);
ftime=toc;

py=predict(mdl,xvalid);
dif=yvalid-py;
er=sum(dif.^2);
rsq=1-er/ssy;
rsq(rsq<0)=0;
