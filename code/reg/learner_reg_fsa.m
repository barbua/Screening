function [rowtime,py]=learner_reg_fsa(xtest,ytest,xtrain,ytrain,corefun,kstart,shrink,eta,nEpochs,mu,minibatch)
ssy=sum((ytest-mean(ytest)).^2);    
tic

[b,b0,sel]=trainLinFSA(xtrain,ytrain,corefun,kstart,shrink,eta,nEpochs,mu,0.9,minibatch,0);

rowtime=toc;
py=xtest(:,sel)*b'+b0;
dif=ytest-py;
er=sum(dif.^2);

rsq=1-er/ssy;
rsq(rsq<0)=0;
