function pertable=cl_percentable(names,values)
pt=prctile(values,90,1);
%order
meanvalue=mean(values,1);
[~,idx]=sort(meanvalue,'descend');
names=names(idx);
pt=pt(idx);

pertable=table(pt(:,1),pt(:,2),pt(:,3),pt(:,4),pt(:,5),pt(:,6),pt(:,7),pt(:,8),pt(:,9),pt(:,10),...
    pt(:,11),pt(:,12),pt(:,13),pt(:,14),pt(:,15),pt(:,16),pt(:,17),pt(:,18),pt(:,19),pt(:,20),...
    pt(:,21),pt(:,22),pt(:,23),pt(:,24),pt(:,25),pt(:,26),pt(:,27),pt(:,28),pt(:,29),pt(:,30),...
    pt(:,31),pt(:,32),pt(:,33),pt(:,34),pt(:,35),pt(:,36),pt(:,37),pt(:,38),pt(:,39),pt(:,40),'VariableNames',names);