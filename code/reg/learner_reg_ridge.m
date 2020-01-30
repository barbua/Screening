function [ftime,py]=learner_reg_ridge(xvalid,yvalid,xtrain,ytrain,shrink)
ssy=sum((yvalid-mean(yvalid)).^2);

tic
b = ridge(ytrain,xtrain,shrink,0);
ftime=toc;

py=[xvalid,ones(size(xvalid,1),1)]*b;
dif=yvalid-py;
er=sum(dif.^2);
rsq=1-er/ssy;
rsq(rsq<0)=0;
