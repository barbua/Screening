clear
nruns=7;
dataname='GLI_85';
load(sprintf('../%s_result.mat',dataname));

if strcmp(dataname , 'GLI_85')
    scalen=1.78;
elseif strcmp(dataname,'dexter')
    scalen=1.78;
elseif strcmp(dataname , 'gisette')
    scalen=1.73;
elseif strcmp(dataname,'SMK_CAN_187')
    scalen=1.78;
elseif strcmp(dataname,'madelon')
    scalen=1.25;
end

for i=1:30   %features selected by filter
      %fsa+relief
      aucfsa_r=[];
      for j=1:i
          for r=1:nruns
            aucfsa_r(j,r)=results.firstLayer.(sprintf('filter_select_%d',i)).secondLayer.(sprintf('filter_select_%d',j)).AUC.FSA.Relief.auc.(sprintf('run_%d',r));
          end
      end
      inner_meanr=mean(aucfsa_r,2);
      [~,mir(i)]=max(inner_meanr);  
      maxaucr(i,:)=aucfsa_r(mir(i),:);
      mr(i)=mean(maxaucr(i,:));
      sr(i)=std(maxaucr(i,:))/sqrt(nruns);
      %fsa+MI
      aucfsa_mu=[];
      for j=1:i
          for r=1:nruns
            aucfsa_mu(j,r)=results.firstLayer.(sprintf('filter_select_%d',i)).secondLayer.(sprintf('filter_select_%d',j)).AUC.FSA.MI.auc.(sprintf('run_%d',r));
          end
      end
      inner_meanmu=mean(aucfsa_mu,2);
      [~,mimu(i)]=max(inner_meanmu);  
      maxaucmu(i,:)=aucfsa_mu(mimu(i),:);
      mmu(i)=mean(maxaucmu(i,:));
      smu(i)=std(maxaucmu(i,:))/sqrt(nruns);
      %fsa+mrmr
      aucfsa_mm=[];
      for j=1:i
          for r=1:nruns
            aucfsa_mm(j,r)=results.firstLayer.(sprintf('filter_select_%d',i)).secondLayer.(sprintf('filter_select_%d',j)).AUC.FSA.MRMR.auc.(sprintf('run_%d',r));
          end
      end
      inner_meanmm=mean(aucfsa_mm,2);
      [~,mimm(i)]=max(inner_meanmm);  
      maxaucmm(i,:)=aucfsa_mm(mimm(i),:);
      mmm(i)=mean(maxaucmm(i,:));
      smm(i)=std(maxaucmm(i,:))/sqrt(nruns);
      %fsa+tscore
      aucfsa_t=[];
      for j=1:i
          for r=1:nruns
            aucfsa_t(j,r)=results.firstLayer.(sprintf('filter_select_%d',i)).secondLayer.(sprintf('filter_select_%d',j)).AUC.FSA.Tscore.auc.(sprintf('run_%d',r));
          end
      end
      inner_meant=mean(aucfsa_t,2);
      [~,mit(i)]=max(inner_meant);  
      maxauct(i,:)=aucfsa_t(mit(i),:);
      mt(i)=mean(maxauct(i,:));
      st(i)=std(maxauct(i,:))/sqrt(nruns);
      %fsa+chi2
      aucfsa_2=[];
      for j=1:i
          for r=1:nruns
            aucfsa_2(j,r)=results.firstLayer.(sprintf('filter_select_%d',i)).secondLayer.(sprintf('filter_select_%d',j)).AUC.FSA.Chi2.auc.(sprintf('run_%d',r));
          end
      end
      inner_mean2=mean(aucfsa_2,2);
      [~,mi2(i)]=max(inner_mean2);  
      maxauc2(i,:)=aucfsa_2(mi2(i),:);
      m2(i)=mean(maxauc2(i,:));
      s2(i)=std(maxauc2(i,:))/sqrt(nruns);
      %fsa+fisher
      aucfsa_fish=[];
      for j=1:i
          for r=1:nruns
            aucfsa_fish(j,r)=results.firstLayer.(sprintf('filter_select_%d',i)).secondLayer.(sprintf('filter_select_%d',j)).AUC.FSA.Fisher.auc.(sprintf('run_%d',r));
          end
      end
      inner_meanfish=mean(aucfsa_fish,2);
      [~,mifish(i)]=max(inner_meanfish);  
      maxaucfish(i,:)=aucfsa_fish(mifish(i),:);
      mfish(i)=mean(maxaucfish(i,:));
      sfish(i)=std(maxaucfish(i,:))/sqrt(nruns);
      %fsa+gini
      aucfsa_gini=[];
      for j=1:i
          for r=1:nruns
            aucfsa_gini(j,r)=results.firstLayer.(sprintf('filter_select_%d',i)).secondLayer.(sprintf('filter_select_%d',j)).AUC.FSA.Gini.auc.(sprintf('run_%d',r));
          end
      end
      inner_meangini=mean(aucfsa_gini,2);
      [~,migini(i)]=max(inner_meangini);  
      maxaucgini(i,:)=aucfsa_gini(migini(i),:);
      mgini(i)=mean(maxaucgini(i,:));
      sgini(i)=std(maxaucgini(i,:))/sqrt(nruns);
      
      %fsa
      for r=1:nruns
          aucfsatest(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.fsatest.auc.(sprintf('run_%d',r));
      end
      mfsa(i)=mean(aucfsatest(i,:));
      sfsa(i)=std(aucfsatest(i,:))/sqrt(nruns);
      %logistic+relief
      for r=1:nruns
          aucrl(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.logist.Relief.auc.(sprintf('run_%d',r));
      end
      mrl(i)=mean(aucrl(i,:));
      srl(i)=std(aucrl(i,:))/sqrt(nruns);
      %logistic+MI
      for r=1:nruns
          aucmul(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.logist.MI.auc.(sprintf('run_%d',r));
      end
      mmul(i)=mean(aucmul(i,:));
      smul(i)=std(aucmul(i,:))/sqrt(nruns);
      %logistic+mrmr
      for r=1:nruns
          aucmml(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.logist.MRMR.auc.(sprintf('run_%d',r));
      end
      mmml(i)=mean(aucmml(i,:));
      smml(i)=std(aucmml(i,:))/sqrt(nruns);
      %logistic+tscore
      for r=1:nruns
          auctl(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.logist.Tscore.auc.(sprintf('run_%d',r));
      end
      mtl(i)=mean(auctl(i,:));
      stl(i)=std(auctl(i,:))/sqrt(nruns);
      %logistic+fisher
      for r=1:nruns
          aucfishl(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.logist.Fisher.auc.(sprintf('run_%d',r));
      end
      mfishl(i)=mean(aucfishl(i,:));
      sfishl(i)=std(aucfishl(i,:))/sqrt(nruns);
      %logistic+chi2
      for r=1:nruns
          auc2l(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.logist.Chi2.auc.(sprintf('run_%d',r));
      end
      m2l(i)=mean(auc2l(i,:));
      s2l(i)=std(auc2l(i,:))/sqrt(nruns);
      %logistic+gini
      for r=1:nruns
          aucginil(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.logist.Gini.auc.(sprintf('run_%d',r));
      end
      mginil(i)=mean(aucginil(i,:));
      sginil(i)=std(aucginil(i,:))/sqrt(nruns);
      %-----------------------------------------------------
      %nby+relief
      for r=1:nruns
          aucrn(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.nby.Relief.auc.(sprintf('run_%d',r));
      end
      mrn(i)=mean(aucrn(i,:));
      srn(i)=std(aucrn(i,:))/sqrt(nruns);
      %nby+mi
      for r=1:nruns
          aucmun(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.nby.MI.auc.(sprintf('run_%d',r));
      end
      mmun(i)=mean(aucmun(i,:));
      smun(i)=std(aucmun(i,:))/sqrt(nruns);
      %nby+mrmr
      for r=1:nruns
          aucmmn(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.nby.MRMR.auc.(sprintf('run_%d',r));
      end
      mmmn(i)=mean(aucmmn(i,:));
      smmn(i)=std(aucmmn(i,:))/sqrt(nruns);
      %nby+tscore
      for r=1:nruns
          auctn(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.nby.Tscore.auc.(sprintf('run_%d',r));
      end
      mtn(i)=mean(auctn(i,:));
      stn(i)=std(auctn(i,:))/sqrt(nruns);
      %nby+fish
      for r=1:nruns
          aucfishn(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.nby.Fisher.auc.(sprintf('run_%d',r));
      end
      mfishn(i)=mean(aucfishn(i,:));
      sfishn(i)=std(aucfishn(i,:))/sqrt(nruns);
      %nby+chi2
      for r=1:nruns
          auc2n(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.nby.Chi2.auc.(sprintf('run_%d',r));
      end
      m2n(i)=mean(auc2n(i,:));
      s2n(i)=std(auc2n(i,:))/sqrt(nruns);
      %nby+gini
      for r=1:nruns
          aucginin(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.nby.Gini.auc.(sprintf('run_%d',r));
      end
      mginin(i)=mean(aucginin(i,:));
      sginin(i)=std(aucginin(i,:))/sqrt(nruns);
      %------------------------------------------
      %svm+relief
      for r=1:nruns
          aucrs(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.svm.Relief.auc.(sprintf('run_%d',r));
      end
      mrs(i)=mean(aucrs(i,:));
      srs(i)=std(aucrs(i,:))/sqrt(nruns);
      %svm+mi
      for r=1:nruns
          aucmus(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.svm.MI.auc.(sprintf('run_%d',r));
      end
      mmus(i)=mean(aucmus(i,:));
      smus(i)=std(aucmus(i,:))/sqrt(nruns);
      %svm+mrmr
      for r=1:nruns
          aucmms(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.svm.MRMR.auc.(sprintf('run_%d',r));
      end
      mmms(i)=mean(aucmms(i,:));
      smms(i)=std(aucmms(i,:))/sqrt(nruns);
      %svm+tscore
      for r=1:nruns
          aucts(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.svm.Tscore.auc.(sprintf('run_%d',r));
      end
      mts(i)=mean(aucts(i,:));
      sts(i)=std(aucts(i,:))/sqrt(nruns);
      %svm+fish
      for r=1:nruns
          aucfishs(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.svm.Fisher.auc.(sprintf('run_%d',r));
      end
      mfishs(i)=mean(aucfishs(i,:));
      sfishs(i)=std(aucfishs(i,:))/sqrt(nruns);
      %svm+chi2
      for r=1:nruns
          auc2s(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.svm.Chi2.auc.(sprintf('run_%d',r));
      end
      m2s(i)=mean(auc2s(i,:));
      s2s(i)=std(auc2s(i,:))/sqrt(nruns);
      %svm+gini
      for r=1:nruns
          aucginis(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.svm.Gini.auc.(sprintf('run_%d',r));
      end
      mginis(i)=mean(aucginis(i,:));
      sginis(i)=std(aucginis(i,:))/sqrt(nruns);
      %-----------------------------------------
      %bt+relief
      for r=1:nruns
          aucrb(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.BT.Relief.auc.(sprintf('run_%d',r));
      end
      mrb(i)=mean(aucrb(i,:));
      srb(i)=std(aucrb(i,:))/sqrt(nruns);
      %bt+mi
      for r=1:nruns
          aucmub(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.BT.MI.auc.(sprintf('run_%d',r));
      end
      mmub(i)=mean(aucmub(i,:));
      smub(i)=std(aucmub(i,:))/sqrt(nruns);
      %bt+mrmr
      for r=1:nruns
          aucmmb(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.BT.MRMR.auc.(sprintf('run_%d',r));
      end
      mmmb(i)=mean(aucmmb(i,:));
      smmb(i)=std(aucmmb(i,:))/sqrt(nruns);
      %bt+tscore
      for r=1:nruns
          auctb(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.BT.Tscore.auc.(sprintf('run_%d',r));
      end
      mtb(i)=mean(auctb(i,:));
      stb(i)=std(auctb(i,:))/sqrt(nruns);
      %bt+fish
      for r=1:nruns
          aucfishb(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.BT.Fisher.auc.(sprintf('run_%d',r));
      end
      mfishb(i)=mean(aucfishb(i,:));
      sfishb(i)=std(aucfishb(i,:))/sqrt(nruns);
      %bt+chi2
      for r=1:nruns
          auc2b(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.BT.Chi2.auc.(sprintf('run_%d',r));
      end
      m2b(i)=mean(auc2b(i,:));
      s2b(i)=std(auc2b(i,:))/sqrt(nruns);
      %bt+gini
      for r=1:nruns
          aucginib(i,r)=results.firstLayer.(sprintf('filter_select_%d',i)).AUC.BT.Gini.auc.(sprintf('run_%d',r));
      end
      mginib(i)=mean(aucginib(i,:));
      sginib(i)=std(aucginib(i,:))/sqrt(nruns);
        
end    
  
%-----------------optimal in each method for comparison------
[hm,im]=max(mmm);
mrmr_fsamax=maxaucmm(im,:);
[hr,ir]=max(mr);
rel_fsamax=maxaucr(ir,:);
[hp,ip]=max(mmu);
mutual_fsamax=maxaucmu(ip,:);
[ht,it]=max(mt);
ttest_fsamax=maxauct(it,:);
[hfish,ifish]=max(mfish);
fisher_fsamax=maxaucfish(ifish,:);
[h2,i2]=max(m2);
chi2_fsamax=maxauc2(i2,:);
[hgini,igini]=max(mgini);
gini_fsamax=maxaucgini(igini,:);


[hrl,irl]=max(mrl);
rel_logmax=aucrl(irl,:);
[hpl,ipl]=max(mmul);
mutual_logmax=aucmul(ipl,:);
[hml,iml]=max(mmml);
mrmr_logmax=aucmml(iml,:);
[htl,itl]=max(mtl);
ttest_logmax=auctl(itl,:);
[hfishl,ifishl]=max(mfishl);
fisher_logmax=aucfishl(ifishl,:);
[h2l,i2l]=max(m2l);
chi2_logmax=auc2l(i2l,:);
[hginil,iginil]=max(mginil);
gini_logmax=aucginil(iginil,:);
%-------------------------------------------2
[hrn,irn]=max(mrn);
rel_nybmax=aucrn(irn,:);
[hpn,ipn]=max(mmun);
mutual_nybmax=aucmun(ipn,:);
[hmn,imn]=max(mmmn);
mrmr_nybmax=aucmmn(imn,:);
[htn,itn]=max(mtn);
ttest_nybmax=auctn(itn,:);
[hfishn,ifishn]=max(mfishn);
fisher_nybmax=aucfishn(ifishn,:);
[h2n,i2n]=max(m2n);
chi2_nybmax=auc2n(i2n,:);
[hginin,iginin]=max(mginin);
gini_nybmax=aucginin(iginin,:);
%-----------------------------------------3
[hrs,irs]=max(mrs);
rel_svmmax=aucrs(irs,:);
[hps,ips]=max(mmus);
mutual_svmmax=aucmus(ips,:);
[hms,ims]=max(mmms);
mrmr_svmmax=aucmms(ims,:);
[hts,its]=max(mts);
ttest_svmmax=aucts(its,:);
[hfishs,ifishs]=max(mfishs);
fisher_svmmax=aucfishs(ifishs,:);
[h2s,i2s]=max(m2s);
chi2_svmmax=auc2s(i2s,:);
[hginis,iginis]=max(mginis);
gini_svmmax=aucginis(iginis,:);
%------------------------------------------4
[hrb,irb]=max(mrb);
rel_boostmax=aucrb(irb,:);
[hpb,ipb]=max(mmub);
mutual_boostmax=aucmub(ipb,:);
[hmb,imb]=max(mmmb);
mrmr_boostmax=aucmmb(imb,:);
[htb,itb]=max(mtb);
ttest_boostmax=auctb(itb,:);
[hfishb,ifishb]=max(mfishb);
fisher_boostmax=aucfishb(ifishb,:);
[h2b,i2b]=max(m2b);
chi2_boostmax=auc2b(i2b,:);
[hginib,iginib]=max(mginib);
gini_boostmax=aucginib(iginib,:);

[hff,iff]=max(mfsa);
fsamax=aucfsatest(iff,:);

for r=1:nruns
    aucl(r)=results.benchmark.logist.auc.(sprintf('run_%d',r));
    aucn(r)=results.benchmark.nby.auc.(sprintf('run_%d',r));
    aucs(r)=results.benchmark.svm.auc.(sprintf('run_%d',r));
    aucb(r)=results.benchmark.BT.auc.(sprintf('run_%d',r));
end
%%P-value table
%-------------put max into matrix--------
namelist={'Fsa','MutualFsa','ReliefFsa','MrmrFsa','TscoreFsa','FisherFsa','Chi2Fsa','GiniFsa',...
    'LogisticReg','MutualLogist','ReliefLogist','MrmrLogist','TscoreLogist','FisherLogist','Chi2Logist','GiniLogist',...
    'NaiveBayes','MutualNaiveB','ReliefNaiveB','MrmrNaiveB','TscoreNaiveB','FisherNaiveB','Chi2NaiveB','GiniNaiveB',...
    'SVM','MutualSVM','ReliefSVM','MrmrSVM','TscoreSVM','FisherSVM','Chi2SVM','GiniSVM',...
    'BoostedDT','MutualBoostedDT','ReliefBoostedDT','MrmrBoostedDT','TscoreBoostedDT','FisherBoostedDT','Chi2BoostedDT','GiniBoostedDT'};
maxvalue=[fsamax',mutual_fsamax',rel_fsamax',mrmr_fsamax',ttest_fsamax',fisher_fsamax',chi2_fsamax',gini_fsamax',...
    aucl',mutual_logmax',rel_logmax',mrmr_logmax',ttest_logmax',fisher_logmax',chi2_logmax',gini_logmax',...
    aucn',mutual_nybmax',rel_nybmax',mrmr_nybmax',ttest_nybmax',fisher_nybmax',chi2_nybmax',gini_nybmax',...
    aucs',mutual_svmmax',rel_svmmax',mrmr_svmmax',ttest_svmmax',fisher_svmmax',chi2_svmmax',gini_svmmax',...
    aucb',mutual_boostmax',rel_boostmax',mrmr_boostmax',ttest_boostmax',fisher_boostmax',chi2_boostmax',gini_boostmax'];
%-------------ttest--------------------
load(sprintf('../data/%s.mat',dataname));
rows=size(X,1);
[mean_table,p_table]=cl_pair_ttest_matrix(maxvalue,namelist,rows);
writetable(p_table,sprintf('../%s_ptable.csv',dataname),'WriteRowNames',true);
writetable(mean_table,sprintf('../%s_meantable.csv',dataname),'WriteRowNames',true);
%--------------percentil----------------
pertable=cl_percentable(namelist,maxvalue,rows);
writetable(pertable,sprintf('../%s_pertable.csv',dataname),'WriteRowNames',true);
%%
relief=[mrb;mrl;mr;mrn;mrs];
mi=[mmub;mmul;mmu;mmun;mmus];
tscore=[mtb;mtl;mt;mtn;mts];
mrmr=[mmmb;mmml;mmm;mmmn;mmms];
fisher=[mfishb;mfishl;mfish;mfishn;mfishs];
chi2=[m2b;m2l;m2;m2n;m2s];
GI=[mginib;mginil;mgini;mginin;mginis];
[nofs,nofsidx]=max([mean(mean(aucn,2)),mean(mean(aucs,2)),mean(mean(aucl,2)),mean(mean(aucb,2))]);
nofsidx
relief=max(relief);
mi=max(mi);
tscore=max(tscore);
mrmr=max(mrmr);
fisher=max(fisher);
chi2=max(chi2);
GI=max(GI);
[reliefdot,reliefidx]=max(relief);
[midot,miidx]=max(mi);
[tscoredot,tscoreidx]=max(tscore);
[mrmrdot,mrmridx]=max(mrmr);
[fisherdot,fisheridx]=max(fisher);
[chi2dot,chi2idx]=max(chi2);
[ginidot,giniidx]=max(GI);
%%
%------------comprehensive plot
h=figure;
t=1:30;
l=round(((t)*4).^scalen);

plot(t,relief,'-','color',[255 0 0]./255,'Linewidth',2,'markersize',13);
hold on;
plot(t,mi,'-','color',[0 255 0]./255,'Linewidth',2,'markersize',13);
plot(t,tscore,'-','color',[0,0,255]./255,'Linewidth',2,'markersize',13);
plot(t,mrmr,'-','color',[255 204 0]./255,'LineWidth',2,'MarkerSize',13);
plot(t,fisher,'-','color',[255 128 0]./255,'LineWidth',2,'MarkerSize',13);
plot(t,chi2,'-','color',[51 153 255]./255,'LineWidth',2,'MarkerSize',13);
plot(t,GI,'-','color',[255 102 178]./255,'LineWidth',2,'MarkerSize',13);
plot(t,mfsa,'--','color','black','Linewidth',2,'markersize',10);
plot(t,mean(aucn)*ones(size(t)),'--','color','green','Linewidth',2,'markersize',10);
plot(t,mean(aucb)*ones(size(t)),'--','color','red','Linewidth',2,'markersize',10);
plot(t,mean(aucl)*ones(size(t)),'--','color','yellow','Linewidth',2,'markersize',10);
plot(t,mean(aucs)*ones(size(t)),'--','color','blue','Linewidth',2,'markersize',10);

scatter(reliefidx,reliefdot,1000,'.','MarkerEdgeColor',[255 0 0]./255,'Linewidth',2);
scatter(miidx,midot,1000,'.','MarkerEdgeColor',[0 255 0]./255,'Linewidth',2);
scatter(tscoreidx,tscoredot,1000,'.','MarkerEdgeColor',[0,0,255]./255,'Linewidth',2);
scatter(mrmridx,mrmrdot,1000,'.','MarkerEdgeColor',[255 204 0]./255,'LineWidth',2);
scatter(fisheridx,fisherdot,1000,'.','MarkerEdgeColor',[255 128 0]./255,'LineWidth',2);
scatter(chi2idx,chi2dot,1000,'.','MarkerEdgeColor',[51 153 255]./255,'LineWidth',2);
scatter(giniidx,ginidot,1000,'.','MarkerEdgeColor',[255 102 178]./255,'LineWidth',2);
scatter(iff,mean(fsamax),1000,'.','MarkerEdgeColor','black','Linewidth',2);
h1=legend('RReliefF','Mutual Inf','T-score','MRMR','Fisher Score','Chi-square','Gini Index','FSA','No FS(NaiveB)','No FS(BT)','No FS(Logistic)','No FS(SVM)','location','SouthEast');%BoostDT,
% h1=legend('RReliefF','Mutual Inf','T-score','MRMR','Fisher Score','Chi-square','Gini Index','FSA','No FS(NaiveB)','No FS(BT)','location','SouthEast');%BoostDT
% h1=legend('RReliefF','Mutual Inf','T-score','MRMR','Fisher Score','Chi-square','Gini Index','FSA','No FS(BT)','No FS(SVM)','location','SouthEast');%,'No FS(Logistic)''No FS(NaiveB)',BoostDT
xticks(t);
xticklabels(string(l));
xtickangle(45);
set(gca,'fontWeight','bold', 'FontSize', 18);
set( h1, 'FontSize', 20 ,'color','none');
xlabel('Number of Features Selected', 'FontSize', 18);
ylabel('AUC', 'FontSize', 18);
title(dataname,'FontSize', 25,'Interpreter', 'none');
set(gcf, 'Position', [0,0,1200,700]); % Maximize figure.

axis([0,30,0.45,0.97]);
saveTightFigure(h,sprintf('../%s_complex',dataname),'png');

%-------------------**************figure for best learner*************-----------------------------
nofname={'Fsa','LogisticReg','NaiveBayes','SVM','BoostedDT'};
nofmax=[fsamax',aucl',aucn',aucs', aucb'];
[~,nofid]=max(mean(nofmax));

h=figure;
t=1:30;
l=round(((t)*4).^scalen);
if nofid==1
    plot(t,mr,'Color',[255 0 0]./255,'LineWidth',2,'MarkerSize',10);%'d-',
    hold on;
    plot(t,mmu,'Color',[0 255 0]./255,'LineWidth',2,'MarkerSize',10);%'o-',
    hold on;
    plot(t,mt,'Color',[0,0,255]./255,'LineWidth',2,'MarkerSize',18);%'s-',
    hold on;
    plot(t,mmm,'Color',[255 204 0]./255,'LineWidth',2,'MarkerSize',10);%'p-',
    hold on;
    plot(t,mfish,'Color',[255 128 0]./255,'LineWidth',2,'MarkerSize',10);%'^-',
    hold on;
    plot(t,m2,'Color',[51 153 255]./255,'LineWidth',2,'MarkerSize',18);%'h-',
    hold on;
    plot(t,mgini,'Color',[255 102 178]./255,'LineWidth',2,'MarkerSize',18);%'x-',
    hold on;
    plot(t,mfsa,'--','Color','black','MarkerSize',10,'LineWidth',2);
    hold on;
    scatter(ir,mr(ir),1000,'.','MarkerEdgeColor',[255 0 0]./255,'Linewidth',2);
    hold on;
    scatter(ip,mmu(ip),1000,'.','MarkerEdgeColor',[0 255 0]./255,'Linewidth',2);
    hold on;
    scatter(it,mt(it),1000,'.','MarkerEdgeColor',[0,0,255]./255,'Linewidth',2);
    hold on;
    scatter(im,mmm(im),1000,'.','MarkerEdgeColor',[255 204 0]./255,'Linewidth',2);
    hold on;
    scatter(ifish,mfish(ifish),1000,'.','MarkerEdgeColor',[255 128 0]./255,'Linewidth',2);
    hold on;
    scatter(i2,m2(i2),1000,'.','MarkerEdgeColor',[51 153 255]./255,'Linewidth',2);
    hold on;
    scatter(igini,mgini(igini),1000,'.','MarkerEdgeColor',[255 102 178]./255,'Linewidth',2);
    hold on;
    scatter(iff,mean(fsamax),1000,'.','MarkerEdgeColor','black','Linewidth',2);
    h1=legend('Relief/FSA','Mutual Inf./FSA','T-score/FSA','MRMR/FSA','Fisher Score/FSA','Chi-square/FSA','Gini Index/FSA','FSA','Location','southeast');
    xticks(t);
    xticklabels(string(l));
    xtickangle(45);
    set(gca,'fontWeight','bold', 'FontSize', 18);
    set( h1, 'FontSize', 20 ,'color','none');
    xlabel('Number of Features Selected', 'FontSize', 18);
    ylabel('AUC', 'FontSize', 18);
    tname=sprintf('%s_%s',dataname,string(nofname(1)));
    title(tname,'FontSize', 25,'Interpreter', 'none');
    set(gcf, 'Position', [0,0,1200,700]);
elseif nofid==2
    plot(t,mrl,'Color',[255 0 0]./255,'LineWidth',2,'MarkerSize',10);%'d-',
    hold on;
    plot(t,mmul,'Color',[0 255 0]./255,'LineWidth',2,'MarkerSize',10);%'o-',
    hold on;
    plot(t,mtl,'Color',[0,0,255]./255,'LineWidth',2,'MarkerSize',18);%'s-',
    hold on;
    plot(t,mmml,'Color',[255 204 0]./255,'LineWidth',2,'MarkerSize',10);%'p-',
    hold on;
    plot(t,mfishl,'Color',[255 128 0]./255,'LineWidth',2,'MarkerSize',10);%'^-',
    hold on;
    plot(t,m2l,'Color',[51 153 255]./255,'LineWidth',2,'MarkerSize',18);%'h-',
    hold on;
    plot(t,mginil,'Color',[255 102 178]./255,'LineWidth',2,'MarkerSize',18);%'x-',
    hold on;
    plot(t,mean(aucl)*ones(size(t)),'--','color','black','Linewidth',2,'markersize',10);
    hold on;
    scatter(irl,mrl(irl),1000,'.','MarkerEdgeColor',[255 0 0]./255,'Linewidth',2);
    hold on;
    scatter(ipl,mmul(ipl),1000,'.','MarkerEdgeColor',[0 255 0]./255,'Linewidth',2);
    hold on;
    scatter(itl,mtl(itl),1000,'.','MarkerEdgeColor',[0,0,255]./255,'Linewidth',2);
    hold on;
    scatter(iml,mmml(iml),1000,'.','MarkerEdgeColor',[255 204 0]./255,'Linewidth',2);
    hold on;
    scatter(ifishl,mfishl(ifishl),1000,'.','MarkerEdgeColor',[255 128 0]./255,'Linewidth',2);
    hold on;
    scatter(i2l,m2l(i2l),1000,'.','MarkerEdgeColor',[51 153 255]./255,'Linewidth',2);
    hold on;
    scatter(iginil,mginil(iginil),1000,'.','MarkerEdgeColor',[255 102 178]./255,'Linewidth',2);
    hold on;
    h1=legend('Relief/LogistReg','Mutual Inf./LogistReg','T-score/LogistReg','MRMR/LogistReg','Fisher Score/LogistReg','Chi-square/LogistReg','Gini Index/LogistReg','LogistReg','Location','southeast');
    xticks(t);
    xticklabels(string(l));
    xtickangle(45);
    set(gca,'fontWeight','bold', 'FontSize', 18);
    set( h1, 'FontSize', 20 ,'color','none');
    xlabel('Number of Features Selected', 'FontSize', 18);
    ylabel('AUC', 'FontSize', 18);
    tname=sprintf('%s_%s',dataname,string(nofname(2)));
    title(tname,'FontSize', 25,'Interpreter', 'none');
    set(gcf, 'Position', [0,0,1200,700]);
elseif nofid==3
    plot(t,mrn,'Color',[255 0 0]./255,'LineWidth',2,'MarkerSize',10);%'d-',
    hold on;
    plot(t,mmun,'Color',[0 255 0]./255,'LineWidth',2,'MarkerSize',10);%'o-',
    hold on;
    plot(t,mtn,'Color',[0,0,255]./255,'LineWidth',2,'MarkerSize',18);%'s-',
    hold on;
    plot(t,mmmn,'Color',[255 204 0]./255,'LineWidth',2,'MarkerSize',10);%'p-',
    hold on;
    plot(t,mfishn,'Color',[255 128 0]./255,'LineWidth',2,'MarkerSize',10);%'^-',
    hold on;
    plot(t,m2n,'Color',[51 153 255]./255,'LineWidth',2,'MarkerSize',18);%'h-',
    hold on;
    plot(t,mginin,'Color',[255 102 178]./255,'LineWidth',2,'MarkerSize',18);%'x-',
    hold on;
    plot(t,mean(aucn)*ones(size(t)),'--','color','black','Linewidth',2,'markersize',10);
    hold on;
    scatter(irn,mrn(irn),1000,'.','MarkerEdgeColor',[255 0 0]./255,'Linewidth',2);
    hold on;
    scatter(ipn,mmun(ipn),1000,'.','MarkerEdgeColor',[0 255 0]./255,'Linewidth',2);
    hold on;
    scatter(itn,mtn(itn),1000,'.','MarkerEdgeColor',[0,0,255]./255,'Linewidth',2);
    hold on;
    scatter(imn,mmmn(imn),1000,'.','MarkerEdgeColor',[255 204 0]./255,'Linewidth',2);
    hold on;
    scatter(ifishn,mfishn(ifishn),1000,'.','MarkerEdgeColor',[255 128 0]./255,'Linewidth',2);
    hold on;
    scatter(i2n,m2n(i2n),1000,'.','MarkerEdgeColor',[51 153 255]./255,'Linewidth',2);
    hold on;
    scatter(iginin,mginin(iginin),1000,'.','MarkerEdgeColor',[255 102 178]./255,'Linewidth',2);
    hold on;
    h1=legend('Relief/NaiveBy','Mutual Inf./NaiveBy','T-score/NaiveBy','MRMR/NaiveBy','Fisher Score/NaiveBy','Chi-square/NaiveBy','Gini Index/NaiveBy','NaiveBy','Location','southeast');
    xticks(t);
    xticklabels(string(l));
    xtickangle(45);
    set(gca,'fontWeight','bold', 'FontSize', 18);
    set( h1, 'FontSize', 20 ,'color','none');
    xlabel('Number of Features Selected', 'FontSize', 18);
    ylabel('AUC', 'FontSize', 18);
    tname=sprintf('%s_%s',dataname,string(nofname(3)));
    title(tname,'FontSize', 25,'Interpreter', 'none');
    set(gcf, 'Position', [0,0,1200,700]);
elseif nofid==4
    plot(t,mrs,'Color',[255 0 0]./255,'LineWidth',2,'MarkerSize',10);%'d-',
    hold on;
    plot(t,mmus,'Color',[0 255 0]./255,'LineWidth',2,'MarkerSize',10);%'o-',
    hold on;
    plot(t,mts,'Color',[0,0,255]./255,'LineWidth',2,'MarkerSize',18);%'s-',
    hold on;
    plot(t,mmms,'Color',[255 204 0]./255,'LineWidth',2,'MarkerSize',10);%'p-',
    hold on;
    plot(t,mfishs,'Color',[255 128 0]./255,'LineWidth',2,'MarkerSize',10);%'^-',
    hold on;
    plot(t,m2s,'Color',[51 153 255]./255,'LineWidth',2,'MarkerSize',18);%'h-',
    hold on;
    plot(t,mginis,'Color',[255 102 178]./255,'LineWidth',2,'MarkerSize',18);%'x-',
    hold on;
    plot(t,mean(aucs)*ones(size(t)),'--','color','black','Linewidth',2,'markersize',10);
    hold on;
    scatter(irs,mrs(irs),1000,'.','MarkerEdgeColor',[255 0 0]./255,'Linewidth',2);
    hold on;
    scatter(ips,mmus(ips),1000,'.','MarkerEdgeColor',[0 255 0]./255,'Linewidth',2);
    hold on;
    scatter(its,mts(its),1000,'.','MarkerEdgeColor',[0,0,255]./255,'Linewidth',2);
    hold on;
    scatter(ims,mmms(ims),1000,'.','MarkerEdgeColor',[255 204 0]./255,'Linewidth',2);
    hold on;
    scatter(ifishs,mfishs(ifishs),1000,'.','MarkerEdgeColor',[255 128 0]./255,'Linewidth',2);
    hold on;
    scatter(i2s,m2s(i2s),1000,'.','MarkerEdgeColor',[51 153 255]./255,'Linewidth',2);
    hold on;
    scatter(iginis,mginis(iginis),1000,'.','MarkerEdgeColor',[255 102 178]./255,'Linewidth',2);
    hold on;
    h1=legend('Relief/SVM','Mutual Inf./SVM','T-score/SVM','MRMR/SVM','Fisher Score/SVM','Chi-square/SVM','Gini Index/SVM','SVM','Location','southeast');
    xticks(t);
    xticklabels(string(l));
    xtickangle(45);
    set(gca,'fontWeight','bold', 'FontSize', 18);
    set( h1, 'FontSize', 20 ,'color','none');
    xlabel('Number of Features Selected', 'FontSize', 18);
    ylabel('AUC', 'FontSize', 18);
    tname=sprintf('%s_%s',dataname,string(nofname(4)));
    title(tname,'FontSize', 25,'Interpreter', 'none');
    set(gcf, 'Position', [0,0,1200,700]);
elseif nofid==5
    plot(t,mrb,'Color',[255 0 0]./255,'LineWidth',2,'MarkerSize',10);%'d-',
    hold on;
    plot(t,mmub,'Color',[0 255 0]./255,'LineWidth',2,'MarkerSize',10);%'o-',
    hold on;
    plot(t,mtb,'Color',[0,0,255]./255,'LineWidth',2,'MarkerSize',18);%'s-',
    hold on;
    plot(t,mmmb,'Color',[255 204 0]./255,'LineWidth',2,'MarkerSize',10);%'p-',
    hold on;
    plot(t,mfishb,'Color',[255 128 0]./255,'LineWidth',2,'MarkerSize',10);%'^-',
    hold on;
    plot(t,m2b,'Color',[51 153 255]./255,'LineWidth',2,'MarkerSize',18);%'h-',
    hold on;
    plot(t,mginib,'Color',[255 102 178]./255,'LineWidth',2,'MarkerSize',18);%'x-',
    hold on;
    plot(t,mean(aucb)*ones(size(t)),'--','color','black','Linewidth',2,'markersize',10);
    hold on;
    scatter(irb,mrb(irb),1000,'.','MarkerEdgeColor',[255 0 0]./255,'Linewidth',2);
    hold on;
    scatter(ipb,mmub(ipb),1000,'.','MarkerEdgeColor',[0 255 0]./255,'Linewidth',2);
    hold on;
    scatter(itb,mtb(itb),1000,'.','MarkerEdgeColor',[0,0,255]./255,'Linewidth',2);
    hold on;
    scatter(imb,mmmb(imb),1000,'.','MarkerEdgeColor',[255 204 0]./255,'Linewidth',2);
    hold on;
    scatter(ifishb,mfishb(ifishb),1000,'.','MarkerEdgeColor',[255 128 0]./255,'Linewidth',2);
    hold on;
    scatter(i2b,m2b(i2b),1000,'.','MarkerEdgeColor',[51 153 255]./255,'Linewidth',2);
    hold on;
    scatter(iginib,mginib(iginib),1000,'.','MarkerEdgeColor',[255 102 178]./255,'Linewidth',2);
    hold on;
    h1=legend('Relief/BoostedDT','Mutual Inf./BoostedDT','T-score/BoostedDT','MRMR/BoostedDT','Fisher Score/BoostedDT','Chi-square/BoostedDT','Gini Index/BoostedDT','BoostedDT','Location','southeast');
    xticks(t);
    xticklabels(string(l));
    xtickangle(45);
    set(gca,'fontWeight','bold', 'FontSize', 18);
    set( h1, 'FontSize', 20 ,'color','none');
    xlabel('Number of Features Selected', 'FontSize', 18);
    ylabel('AUC', 'FontSize', 18);
    tname=sprintf('%s_%s',dataname,string(nofname(5)));
    title(tname,'FontSize', 25,'Interpreter', 'none');
    set(gcf, 'Position', [0,0,1200,700]);
end
axis([0,30,0.45,0.97]);
saveTightFigure(h,sprintf('../%s_detail.png',dataname),'png');
