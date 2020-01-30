function [rowtime,py]=learner_cl_fsa(xtest,ytest,xtrain,ytrain,corefun,kstart,shrink,eta,nEpochs,mu,minibatch)

tic
[bb,b0,sel]=trainLinFSA(xtrain,ytrain,corefun,kstart,shrink,eta,nEpochs,mu,0.9,minibatch,0);
rowtime=toc;

py=xtest(:,sel)*bb'+b0;%S.Intercept;
py=exp(py)./(1+exp(py));
[~,~,~,auc] = perfcurve(ytest,py,1);

py=py.*ytest;
er=mean(py<0);
    
    
    

