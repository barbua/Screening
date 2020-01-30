function pertable=reg_percentable(names,values)
pt=prctile(values,90,1);
%order
meanvalue=mean(values,1);
[~,idx]=sort(meanvalue,'descend');
names=names(idx);
pt=pt(idx);

pertable=table(pt(:,1),pt(:,2),pt(:,3),pt(:,4),pt(:,5),pt(:,6),pt(:,7),pt(:,8),pt(:,9),pt(:,10),pt(:,11),pt(:,12),'VariableNames',names);%