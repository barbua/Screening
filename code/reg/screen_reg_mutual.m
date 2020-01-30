function [scores,rankidx]=screen_reg_mutual(X,y,bin)
parfor p=1:size(X,2)
    if sum((X(:,p)-mean(X(:,p))).^2)==0
        scores(p)=0;
    else
        scores(p)=muinf(y,X(:,p),bin);
    end
end
[B,I]=sort(scores,'descend');
rankidx=I;
end


function I=muinf(y,x,bin)
point_num=size(y,1);
maxx=max(x);
minx=min(x);
maxy=max(y);
miny=min(y);
xi=(maxx-minx)/bin;%interval
yi=(maxy-miny)/bin;

for b=1:bin+1
    pgrid(b,1)=miny+yi*(b-1);
    pgrid(b,2)=minx+xi*(b-1);
end

yx=[y,x];
for i=1:bin
    yxlist(yx(:,1)>=pgrid(i,1) & yx(:,1)<pgrid(i+1,1)+0.00001,1)=i;
    yxlist(yx(:,2)>=pgrid(i,2) & yx(:,2)<pgrid(i+1,2)+0.00001,2)=i;
end
[uqr,ia,ic]=unique(yxlist,'rows');
newyx=[yxlist,ic];
newuqr=unique(newyx,'rows');
for c=1:size(uqr,1)
    pyx=sum(newyx(:,3)==c)/point_num;
    iy=newuqr(newuqr(:,3)==c,1);
    ix=newuqr(newuqr(:,3)==c,2);
    py=sum(newyx(:,1)==iy)/point_num;
    px=sum(newyx(:,2)==ix)/point_num;
    ele(c)=pyx*log(pyx/(px*py));
end
I=sum(ele);
end