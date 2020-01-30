function [mean_table,p_table]=cl_pair_ttest_matrix(maxvalue,namelist,rows)
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

p_table=table(ps(:,1),ps(:,2),ps(:,3),ps(:,4),ps(:,5),ps(:,6),ps(:,7),ps(:,8),...
    ps(:,9),ps(:,10),ps(:,11),ps(:,12),ps(:,13),ps(:,14),ps(:,15),...
    ps(:,16),ps(:,17),ps(:,18),ps(:,19),ps(:,20),ps(:,21),ps(:,22),...
    ps(:,23),ps(:,24),ps(:,25),ps(:,26),ps(:,27),ps(:,28),ps(:,29),...
    ps(:,30),ps(:,31),ps(:,32),ps(:,33),ps(:,34),ps(:,35),ps(:,36),...
    ps(:,37),ps(:,38),ps(:,39),ps(:,40),...
    'RowNames',namelist,'VariableNames',namelist);

mean_table=table(sortv(:,1),sortv(:,2),sortv(:,3),sortv(:,4),sortv(:,5),sortv(:,6),sortv(:,7),sortv(:,8),...
    sortv(:,9),sortv(:,10),sortv(:,11),sortv(:,12),sortv(:,13),sortv(:,14),sortv(:,15),...
    sortv(:,16),sortv(:,17),sortv(:,18),sortv(:,19),sortv(:,20),sortv(:,21),sortv(:,22),...
    sortv(:,23),sortv(:,24),sortv(:,25),sortv(:,26),sortv(:,27),sortv(:,28),sortv(:,29),...
    sortv(:,30),sortv(:,31),sortv(:,32),sortv(:,33),sortv(:,34),sortv(:,35),sortv(:,36),...
    sortv(:,37),sortv(:,38),sortv(:,39),sortv(:,40),...
    'VariableNames',namelist);