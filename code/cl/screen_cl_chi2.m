function [rank2,chi2]=screen_cl_chi2(x,y)
for i=1:size(x,2)
    [~,chi2(i),~] = crosstab(y,x(:,i));
end
chi2(isnan(chi2)) = 0;
[~,rank2]=sort(chi2,'descend');