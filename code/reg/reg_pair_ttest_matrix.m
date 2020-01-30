function [mean_table,p_table]=reg_pair_ttest_matrix(maxvalue,namelist,rows)
%'maxvalue': column matrix,'namelist': cell array

k=size(maxvalue,2);
%order methods according to mean
meanmaxvalue=mean(maxvalue,1);
[sortv,maxidx]=sort(meanmaxvalue,'descend');
namelist=namelist(maxidx);
maxvalue=maxvalue(:,maxidx);
sortv=[sortv;std(maxvalue)./sqrt(rows)];

% -get matrix
for i=1:k
   for j=1:k
       [~,temp]=ttest(maxvalue(:,i),maxvalue(:,j));
       if temp<0.05
           p(i,j)=0;
       elseif isnan(temp)
           p(i,j)=1;
       elseif i==j
           p(i,j)=1;
       else
           p(i,j)=round(temp,4);
       end
   end
end

for i=1:k
   for j=1:k
       if p(i,j)==0
           ps(i,j)=string('<0.05');
       else
           ps(i,j)=string(p(i,j));
       end
   end
end

p_table=table(ps(:,1),ps(:,2),ps(:,3),ps(:,4),ps(:,5),ps(:,6),ps(:,7),ps(:,8),ps(:,9),ps(:,10),ps(:,11),ps(:,12),'RowNames',namelist,'VariableNames',namelist);%

mean_table=table(sortv(:,1),sortv(:,2),sortv(:,3),sortv(:,4),sortv(:,5),sortv(:,6),sortv(:,7),sortv(:,8),...
    sortv(:,9),sortv(:,10),sortv(:,11),sortv(:,12),'VariableNames',namelist);%