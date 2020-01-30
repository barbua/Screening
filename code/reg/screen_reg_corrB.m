function [abscorry,I]=screen_reg_corrB(x,y)

for i=1:size(x,2)
    corry(1,i)=corr(x(:,i),y);
    
end
abscorry=abs(corry);
abscorry(isnan(abscorry))=0;
[B,I]=sort(abscorry,'descend');
