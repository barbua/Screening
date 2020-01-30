clear
nruns=7;
dataname='indoorloc';
load(sprintf('../%s_result.mat',dataname));

outter_loop=30;
if strcmp(dataname , 'BMI')
    scalen=1.825;
elseif strcmp(dataname,'Tumor')
    scalen=1.825;
elseif strcmp(dataname , 'CoEPrA20063')
    scalen=1.68;
elseif strcmp(dataname,'wikidataxyvggface')
    scalen=1.65;
    outter_loop=23;
elseif strcmp(dataname,'indoorloc')
    scalen=1.25;
end


for i=1:outter_loop
    %bt+relief
    for r=1:nruns
      rsqrb(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).rsq.BT.Relief.rsq.(sprintf('run_%d',r));
    end    
    mrb(i)=mean(rsqrb(i,:));
    srb(i)=std(rsqrb(i,:))/sqrt(nruns);
    %bt+corr
    for r=1:nruns
      rsqcb(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).rsq.BT.corr.rsq.(sprintf('run_%d',r));
    end 
    mcb(i)=mean(rsqcb(i,:));
    scb(i)=std(rsqcb(i,:))/sqrt(nruns);
    %bt+mi
    for r=1:nruns
      rsqmub(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).rsq.BT.MI.rsq.(sprintf('run_%d',r));
    end 
    mmub(i)=mean(rsqmub(i,:));
    smub(i)=std(rsqmub(i,:))/sqrt(nruns);
    %ridge+relief
    for r=1:nruns
      rsqrr(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).rsq.ridge.Relief.rsq.(sprintf('run_%d',r));
    end 
    mrr(i)=mean(rsqrr(i,:));
    srr(i)=std(rsqrr(i,:))/sqrt(nruns);
    %ridge+corr
    for r=1:nruns
      rsqcr(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).rsq.ridge.corr.rsq.(sprintf('run_%d',r));
    end
    mcr(i)=mean(rsqcr(i,:));
    scr(i)=std(rsqcr(i,:))/sqrt(nruns);
    %ridge+mi
    for r=1:nruns
      rsqmur(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).rsq.ridge.MI.rsq.(sprintf('run_%d',r));
    end
    mmur(i)=mean(rsqmur(i,:));
    smur(i)=std(rsqmur(i,:))/sqrt(nruns);
    %fsa
    for r=1:nruns
        rsqfsatest(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).rsq.fsatest.rsq.(sprintf('run_%d',r));
    end
    mfsa(i)=mean(rsqfsatest(i,:));
    sfsa(i)=std(rsqfsatest(i,:))/sqrt(nruns);
    

    %fsa+relief
    rsqfsa_r=[];
    for j=1:i
        for r=1:nruns
            rsqfsa_r(j,r)=results.firstLayer.(sprintf('filter_select_%d',i)).secondLayer.(sprintf('filter_select_%d',j)).rsq.FSA.Relief.rsq.(sprintf('run_%d',r));
        end
    end
    inner_meanr=mean(rsqfsa_r,2);
    [~,mir(i)]=max(inner_meanr);
    maxrsqr(i,:)=rsqfsa_r(mir(i),:);
    mr(i)=mean(maxrsqr(i,:));
    sr(i)=std(maxrsqr(i,:))/sqrt(nruns);
    %fsa+mi
    rsqfsa_mu=[];
    for j=1:i
        for r=1:nruns
            rsqfsa_mu(j,r)=results.firstLayer.(sprintf('filter_select_%d',i)).secondLayer.(sprintf('filter_select_%d',j)).rsq.FSA.MI.rsq.(sprintf('run_%d',r));
        end
    end
    inner_meanmu=mean(rsqfsa_mu,2);
    [~,mimu(i)]=max(inner_meanmu);
    maxrsqmu(i,:)=rsqfsa_mu(mimu(i),:);
    mmu(i)=mean(maxrsqmu(i,:));
    smu(i)=std(maxrsqmu(i,:))/sqrt(nruns);
    %fsa+mi
    rsqfsa_c=[];
    for j=1:i
        for r=1:nruns
            rsqfsa_c(j,r)=results.firstLayer.(sprintf('filter_select_%d',i)).secondLayer.(sprintf('filter_select_%d',j)).rsq.FSA.corr.rsq.(sprintf('run_%d',r));
        end
    end
    inner_meanc=mean(rsqfsa_c,2);
    [~,mic(i)]=max(inner_meanc);
    maxrsqc(i,:)=rsqfsa_c(mic(i),:);
    mc(i)=mean(maxrsqc(i,:));
    sc(i)=std(maxrsqc(i,:))/sqrt(nruns);
end
%-----------------optimal in each method------
[hc,ic]=max(mc);
cor_fsamax=maxrsqc(ic,:);
[hr,ir]=max(mr);
rel_fsamax=maxrsqr(ir,:);
[hp,ip]=max(mmu);
mutual_fsamax=maxrsqmu(ip,:);
[hrr,irr]=max(mrr);
rel_ridmax=rsqrr(irr,:);
[hpr,ipr]=max(mmur);
mutual_ridmax=rsqmur(ipr,:);
[hcr,icr]=max(mcr);
cor_ridmax=rsqcr(icr,:);
[hrb,irb]=max(mrb);
rel_boomax=rsqrb(irb,:);
[hpb,ipb]=max(mmub);
mutual_boomax=rsqmub(ipb,:);
[hcb,icb]=max(mcb);
cor_boomax=rsqcb(icb,:);
[hff,iff]=max(mfsa);
fsamax=rsqfsatest(iff,:);
%%P-value table
%%-------------put max into matrix--------
for r=1:nruns
    rsqridge(r)=results.benchmark.ridge.rsq.(sprintf('run_%d',r));
    rsqbt(r)=results.benchmark.BT.rsq.(sprintf('run_%d',r));
end
namelist={'Fsa','CorrelationFsa','ReliefFsa','MutualInfoFsa','Ridge','ReliefRidge','MutualInfoRidge','CorrelationRidge','BoostedRegTree','ReliefBoostedTree','MutualInfoBoostedTree','CorrelationBoostedTree'};
maxvalue=[fsamax',cor_fsamax',rel_fsamax',mutual_fsamax',rsqridge',rel_ridmax',mutual_ridmax',cor_ridmax',rsqbt',rel_boomax',mutual_boomax',cor_boomax'];
%-------------ttest--------------------
load(sprintf('../data/%s.mat',dataname));
rows=size(X,1);
[mean_table,p_table]=reg_pair_ttest_matrix(maxvalue,namelist,rows);
writetable(p_table,sprintf('../%s_ptable.csv',dataname),'WriteRowNames',true);
writetable(mean_table,sprintf('../%s_meantable.csv',dataname),'WriteRowNames',true);
%--------------percentil----------------
pertable=reg_percentable(namelist,maxvalue);
writetable(pertable,sprintf('../%s_pertable.csv',dataname),'WriteRowNames',true);


relief=[mrb;mrr;mr];
mi=[mmub;mmur;mmu];
cor=[mcb;mcr;mc];
[nofs,nofsidx]=max([mean(rsqbt),mean(rsqridge)]);
nofsidx
relief=max(relief);
mi=max(mi);
cor=max(cor);
[reliefdot,reliefidx]=max(relief);
[midot,miidx]=max(mi);
[cordot,coridx]=max(cor);
%%
%------------comprehensive plot
h=figure;
t=1:outter_loop;
l=round(((t)*4).^scalen);

plot(t,relief,'-','color',[255 0 0]./255,'Linewidth',2,'markersize',13);
hold on;
plot(t,mi,'-','color',[0 255 0]./255,'Linewidth',2,'markersize',13);
hold on;
plot(t,cor,'-','color',[0,0,255]./255,'Linewidth',2,'markersize',13);
plot(t,mfsa,'--','color','black','Linewidth',2,'markersize',10);
plot(t,mean(rsqridge)*ones(size(t)),'--','color','red','Linewidth',2,'markersize',10);
plot(t,mean(rsqbt)*ones(size(t)),'--','color','green','Linewidth',2,'markersize',10);

scatter(reliefidx,reliefdot,1000,'.','MarkerEdgeColor',[255 0 0]./255,'Linewidth',2);
scatter(miidx,midot,1000,'.','MarkerEdgeColor',[0 255 0]./255,'Linewidth',2);
scatter(coridx,cordot,1000,'.','MarkerEdgeColor',[0,0,255]./255,'Linewidth',2);
scatter(iff,mean(fsamax),1000,'.','MarkerEdgeColor','black','Linewidth',2);

h1=legend('RReliefF','Mutual Inf','Correlation','FSA','No FS(Ridge)','No FS(BT)','location','SouthEast');
xticks(t);
xticklabels(string(l));
xtickangle(45)
set(gca,'fontWeight','bold', 'FontSize', 18);
set( h1, 'FontSize', 20 ,'color','none');
xlabel('Number of Features Selected', 'FontSize', 18);
ylabel('R-square', 'FontSize', 18);
tname=dataname;
if strcmp(dataname,'wikidataxyvggface')
    tname='Wikiface';
end
title(tname,'FontSize', 25);
set(gcf, 'Position', [0,0,1200,700]); % Maximize figure.

axis([0,30,0.8,1]);
saveTightFigure(h,sprintf('../%s_complex',dataname),'png');

%-------------------**************figure for best learner*************-----------------------------
nofname={'Fsa','Ridge','BoostedRegTree'};
nofmax=[fsamax',rsqridge',rsqbt'];
[~,nofid]=max(mean(nofmax));

h=figure;
t=1:outter_loop;
l=round(((t)*4).^scalen);

if nofid==1
    plot(t,mr,'color',[255 0 0]./255,'Linewidth',2,'markersize',13);
    hold on;
    plot(t,mmu,'color',[0 255 0]./255,'Linewidth',2,'markersize',25);
    hold on;
    plot(t,mc,'color',[0,0,255]./255,'Linewidth',2,'markersize',8);
    hold on;
    plot(t,mfsa,'--','color','black','Linewidth',2,'markersize',10);
    hold on;
    scatter(ir,mr(ir),1000,'.','MarkerEdgeColor',[255 0 0]./255,'Linewidth',2);
    hold on;
    scatter(ip,mmu(ip),1000,'.','MarkerEdgeColor',[0 255 0]./255,'Linewidth',2);
    hold on;
    scatter(ic,mc(ic),1000,'.','MarkerEdgeColor',[0,0,255]./255,'Linewidth',2);
    hold on;
    scatter(iff,mean(fsamax),1000,'.','MarkerEdgeColor','black','Linewidth',2);
    h1=legend('RReliefF/FSA','Mutual Inf/FSA','Correlation/FSA','FSA','location','SouthEast');
    xticks(t);
    xticklabels(string(l));
    xtickangle(45);
    set(gca,'fontWeight','bold', 'FontSize', 18);
    set( h1, 'FontSize', 20 ,'color','none');
    xlabel('Number of Features Selected', 'FontSize', 18);
    ylabel('R-square', 'FontSize', 18);
    
    if strcmp(dataname,'wikidataxyvggface')
        tname=sprintf('%s_%s','Wikiface',string(nofname(1)));
    else
        tname=sprintf('%s_%s',dataname,string(nofname(1)));
    end
    title(tname,'FontSize', 25,'Interpreter', 'none');
    set(gcf, 'Position', [0,0,1200,700]);       
elseif nofid==2
    plot(t,mrr,'color',[255 0 0]./255,'Linewidth',2,'markersize',13);
    hold on;
    plot(t,mmur,'color',[0 255 0]./255,'Linewidth',2,'markersize',25);
    hold on;
    plot(t,mcr,'color',[0,0,255]./255,'Linewidth',2,'markersize',8);
    hold on;
    plot(t,mean(rsqridge)*ones(size(t)),'--','color','black','Linewidth',2,'markersize',8);
    hold on;
    scatter(irr,mean(rel_ridmax),1000,'.','MarkerEdgeColor',[255 0 0]./255,'Linewidth',2);
    hold on;
    scatter(ipr,mean(mutual_ridmax),1000,'.','MarkerEdgeColor',[0 255 0]./255,'Linewidth',2);
    hold on;
    scatter(icr,mean(cor_ridmax),1000,'.','MarkerEdgeColor',[0,0,255]./255,'Linewidth',2);
    h1=legend('RReliefF/Ridge','Mutual Inf/Ridge','Correlation/Ridge','Ridge','location','SouthEast');
    xticks(t);
    xticklabels(string(l));
    xtickangle(45);
    set(gca,'fontWeight','bold', 'FontSize', 18);
    set( h1, 'FontSize', 20 ,'color','none');
    xlabel('Number of Features Selected', 'FontSize', 18);
    ylabel('R-square', 'FontSize', 18);
    
    if strcmp(dataname,'wikidataxyvggface')
        tname=sprintf('%s_%s','Wikiface',string(nofname(2)));
    else
        tname=sprintf('%s_%s',dataname,string(nofname(2)));
    end
    title(tname,'FontSize', 25,'Interpreter', 'none');
    set(gcf, 'Position', [0,0,1200,700]);
elseif nofid==3
    plot(t,mrb,'color',[255 0 0]./255,'Linewidth',2,'markersize',13);
    hold on;
    plot(t,mmub,'color',[0 255 0]./255,'Linewidth',2,'markersize',25);
    hold on;
    plot(t,mcb,'color',[0,0,255]./255,'Linewidth',2,'markersize',8);
    hold on;
    plot(t,mean(rsqbt)*ones(size(t)),'--','color','black','Linewidth',2,'markersize',8);
    hold on;
    scatter(irb,mean(rel_boomax),1000,'.','MarkerEdgeColor',[255 0 0]./255,'Linewidth',2);
    hold on;
    scatter(ipb,mean(mutual_boomax),1000,'.','MarkerEdgeColor',[0 255 0]./255,'Linewidth',2);
    hold on;
    scatter(icb,mean(cor_boomax),1000,'.','MarkerEdgeColor',[0,0,255]./255,'Linewidth',2);
    h1=legend('RReliefF/BT','Mutual Inf/BT','Correlation/BT','BoostRegTree','location','SouthEast');
    xticks(t);
    xticklabels(string(l));
    xtickangle(45);
    set(gca,'fontWeight','bold', 'FontSize', 18);
    set( h1, 'FontSize', 20 ,'color','none');
    xlabel('Number of Features Selected', 'FontSize', 18);
    ylabel('R-square', 'FontSize', 18);
    
    if strcmp(dataname,'wikidataxyvggface')
        tname=sprintf('%s_%s','Wikiface',string(nofname(3)));
    else
        tname=sprintf('%s_%s',dataname,string(nofname(3)));
    end
    title(tname,'FontSize', 25,'Interpreter', 'none');
    set(gcf, 'Position', [0,0,1200,700]);
end

axis([0,30,0.8,1]);
saveTightFigure(h,sprintf('../%s_detail.png',dataname),'png');

