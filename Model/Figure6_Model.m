

%folderModelFiles = './DataFiles_Model/';

%% ---------------  read data from file -----------------------------------
T = readtable('OptoBehavioralTable.csv');

% this generates a table with 4 columns, in this order:
% T.RSPTIME=Response Time for animal, in seconds
% T.LASER = Laser (Whether the laser is inhibiting D1/D2 MSNs or not)
% T.ANIMALID = animal ID
% T.GENOTYPE = either D1-cre or D2-cre or control , i.e.
% here only: OptoD1Halo = D1-cre,  OptoD2Halo = D2-cre, 
%            OptoD2CtrlHalo = Control mice but perturbations to D2 MSNs (?)

%% ---------------- keep only the D2-cre mice ---------------------------
D2creMice = (strcmp(T.GENOTYPE,'OptoD2Halo'));
T = T(D2creMice,:);

% if you want to exclude any animal from the analysis, add their ID here
animalsToExclude = [];   % e.g.,   animalsToExclude = {'IOR3','SL10'}; 
%
trialType = {'LaserOff', 'LaserOn'};
colorA={'c','m'};
colorB={'b','r'};

%% calculations
for j=1:length(trialType)
    if strcmp(trialType(j),'LaserOff')
        LaserStatus = 'ctrl';  
        [dataLaserOff,m_LaserOff,md_LaserOff,cv_LaserOff] = AnimalData(T,LaserStatus);
        load([ 'D2_MSNs_LaserOff_ModelExample1.mat']);
        %
        TrialDur_Off =TrialDur; Ymat_Off =Ymat;  samp_path_Off = samp_path;
        baseline_Off = baseline; FR_target_Off = FR_target;  D_Off = D;
        sigma_Off = sigma; t0_Off = t0; tf_Off = tf; h_Off = h;
        clear Ymat samp_path baseline FR_target D sigma t0 tf h
       %-----------------------------------------------
    elseif strcmp(trialType(j),'LaserOn')
        LaserStatus = 'inhibition';    
        [dataLaserOn,m_LaserOn,md_LaserOn, cv_LaserOn] = AnimalData(T,LaserStatus);
        load([ 'D2_MSNs_LaserOn_ModelExample1.mat']);
        %
        TrialDur_On =TrialDur; Ymat_On =Ymat;  samp_path_On = samp_path;
        baseline_On = baseline; FR_target_On = FR_target;  D_On = D;
        sigma_On = sigma; t0_On = t0; tf_On = tf; h_On = h;
        clear Ymat samp_path baseline FR_target D sigma t0 tf h
        % --------------------------------------------
    else
        error('Wrong D2-cre-OptoInhib TrialType')
    end
end


%% ---- test for goodness of  model data with behavioral data -------------
%
% two-sample Kolmogorov-Smirnov test
% h = kstest2(x1,x2) returns a test decision for the null hypothesis that
% the data in vectors x1 and x2 are from the same continuous distribution,
% using the two-sample Kolmogorov-Smirnov test. The alternative hypothesis
% is that x1 and x2 are from different continuous distributions. The result
% h is 1 if the test rejects the null hypothesis at the 5% significance 
% level, and 0 otherwise.
[hvalOff,pvalOff] = kstest2(dataLaserOff,TrialDur_Off,'Alpha',0.05);
disp(['D2 Cre h for ' char(trialType(1)) ' = ' num2str(hvalOff)])
disp(['p for ' char(trialType(1)) ' = ' num2str(pvalOff)])

[hvalOn,pvalOn] = kstest2(dataLaserOn,TrialDur_On,'Alpha',0.05);
disp(['D2 MSN Disrupted h for ' char(trialType(2)) ' = ' num2str(hvalOn)])
disp(['p for ' char(trialType(2)) ' = ' num2str(pvalOn)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------------plot results for behavior ----------------------- %%
xx=linspace(0,20);
xxlabel = 'durations';
ymax=0.15;
nbins=30;

successfulTrials_Off = size(TrialDur_Off(~isnan(TrialDur_Off)),1); 
successfulTrials_On = size(TrialDur_On(~isnan(TrialDur_On)),1); 
set(0,'DefaultFigureWindowStyle','docked');
figure(51)
%----------------
ax1=subplot(3,2,1); 
h=histogram(ax1,dataLaserOff,nbins); 
h.FaceColor = colorA{1};
h.Normalization = 'pdf';
hold on
%
% try fitting with Gamma distrib
pdL1_Off = fitdist(dataLaserOff,'Gamma');
% ---- in shape-rate form -----------------------------------------------
aa=pdL1_Off.a; bb = 1/pdL1_Off.b;
zz =  (bb^aa/ gamma(aa)) .* xx.^(aa-1) .* exp(-bb*xx) ; 
%-------------------------------------------------------------------------
plot(xx,zz,':','color',colorB{1},'Linewidth',2.5)
%
title(['DATA:  D2cre-Opto (N=' num2str(length(dataLaserOff)) ', '...
    num2str(nbins) ' bins): ' char(trialType(1))])
legend('Empirical', '\bf Gamma fit',...
    'Location','northeast','Interpreter','Latex')
legend boxoff
box off
box off
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20])
ylim([0 ymax])
%xlabel(xxlabel)
ylabel('PDF')
text(15.5, 0.07, ['\bf shape $\alpha$=' num2str(pdL1_Off.a,'%1.2f')],'Interpreter','Latex')
text(15.5,0.05,['\bf rate $\beta$=' num2str(1/pdL1_Off.b,'%1.2f')],'Interpreter','Latex')
clear mu sigma h yy zz aa bb 

%-------------%
ax2=subplot(3,2,3); 
h=histogram(ax2,TrialDur_Off,nbins); 
h.FaceColor = [0.7 0.7 0.7] ;
h.Normalization = 'pdf';
hold on
%
% try fitting with Gamma distrib
pdL2_Off = fitdist(TrialDur_Off,'Gamma');
% ---- in shape-rate form -----------------------------------------------
aa=pdL2_Off.a; bb = 1/pdL2_Off.b;
zz =  (bb^aa/ gamma(aa)) .* xx.^(aa-1) .* exp(-bb*xx) ; 
%-------------------------------------------------------------------------
plot(xx,zz,'--k','Linewidth',2)
%
title(['MODEL SIMULATIONS (N=' num2str(successfulTrials_Off) ', '...
    num2str(nbins) ' bins): ' char(trialType(1))])
legend('Empirical','\bf Gamma fit','Location','northeast','Interpreter','Latex')
legend boxoff
box off
box off
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20])
ylim([0 ymax])
%xlabel(xxlabel)
ylabel('PDF')
text(15.5, 0.07, ['\bf shape $\alpha$=' num2str(pdL2_Off.a,'%1.2f')],'Interpreter','Latex')
text(15.5,0.05,['\bf rate $\beta$=' num2str(1/pdL2_Off.b,'%1.2f')],'Interpreter','Latex')
clear mu sigma h yy zz  aa bb

%----------------
ax1=subplot(3,2,2); 
h=histogram(ax1,dataLaserOn,nbins); 
h.FaceColor = colorA{2};
h.Normalization = 'pdf';
hold on
%
% try fitting with Gamma distrib
pdL1_On = fitdist(dataLaserOn,'Gamma');
% ---- in shape-rate form -----------------------------------------------
aa=pdL1_On.a; bb = 1/pdL1_On.b;
zz =  (bb^aa/ gamma(aa)) .* xx.^(aa-1) .* exp(-bb*xx) ; 
%-------------------------------------------------------------------------
plot(xx,zz,':','color',colorB{2},'Linewidth',2.5)
%
title(['DATA:  D2cre-Opto (N=' num2str(length(dataLaserOn)) ', '...
    num2str(nbins) ' bins): ' char(trialType(2))])
legend('Empirical', '\bf Gamma fit',...
    'Location','northeast','Interpreter','Latex')
legend boxoff
box off
box off
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20])
ylim([0 ymax])
%xlabel(xxlabel)
ylabel('PDF')
text(15.5, 0.07, ['\bf shape $\alpha$=' num2str(pdL1_On.a,'%1.2f')],'Interpreter','Latex')
text(15.5, 0.05,['\bf rate $\beta$=' num2str(1/pdL1_On.b,'%1.2f')],'Interpreter','Latex')
clear mu sigma h yy zz aa bb 

%-------------%
ax2=subplot(3,2,4); 
h=histogram(ax2,TrialDur_On,nbins); 
h.FaceColor = [0.7 0.7 0.7] ;
h.Normalization = 'pdf';
hold on
%
% try fitting with Gamma distrib
pdL2_On = fitdist(TrialDur_On,'Gamma');
% ---- in shape-rate form -----------------------------------------------
aa=pdL2_On.a; bb = 1/pdL2_On.b;
zz =  (bb^aa/ gamma(aa)) .* xx.^(aa-1) .* exp(-bb*xx) ; 
%-------------------------------------------------------------------------
plot(xx,zz,'--k','Linewidth',2)
%
title(['MODEL SIMULATIONS (N=' num2str(successfulTrials_On) ', '...
    num2str(nbins) ' bins): ' char(trialType(2))])
legend('Empirical','\bf Gamma fit','Location','northeast','Interpreter','Latex')
legend boxoff
box off
box off
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20])
ylim([0 ymax])
%xlabel(xxlabel)
ylabel('PDF')
text(15.5, 0.07, ['\bf shape $\alpha$=' num2str(pdL2_On.a,'%1.2f')],'Interpreter','Latex')
text(15.5,0.05,['\bf rate $\beta$=' num2str(1/pdL2_On.b,'%1.2f')],'Interpreter','Latex')
clear mu sigma h yy zz  aa bb

%------------
subplot(3,2,5)
h=cdfplot(dataLaserOff);
h.Color= colorA{1};
h.LineWidth=2;
h.LineStyle = '-';
hold on
plot(xx,cdf(pdL1_Off,xx),':','color',colorB{1},'Linewidth',2) % for gamma data
%
plot(xx,cdf(pdL2_Off,xx),'--','color', 'k','Linewidth',1) % for gamma model
title(['DATA/MODEL Comparison - D2cre-Opto: ' char(trialType(1))])
legend('Empirical','\bf Gamma fit', '\bf Model Fit',...
    'Location','northwest','Interpreter','Latex')
box off
legend boxoff
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20])
ylim([-0.05 1.05])
xlabel(xxlabel)
ylabel('CDF')
text(12, 0.60, ['Data: $\alpha$=' num2str(pdL1_Off.a,'%1.2f') ...
    ', $\beta$=' num2str(1/pdL1_Off.b,'%1.2f')],'Color',colorB{1},'Interpreter','Latex')
text(12, 0.45, ['Model: $\alpha$=' num2str(pdL2_Off.a,'%1.2f') ...
    ', $\beta$=' num2str(1/pdL2_Off.b,'%1.2f')],'Interpreter','Latex')
% ----------
text(12, 0.25, ['Data mean/CV: ' num2str(m_LaserOff ,'%1.2f') '; '  ...
    num2str(cv_LaserOff,'%1.2f')],'Color',colorB{1},'Interpreter','Latex')
text(12, 0.10, ['Model mean/CV: ' num2str(pdL2_Off.a*pdL2_Off.b ,'%1.2f') '; '  ...
    num2str(1/sqrt(pdL2_Off.a),'%1.2f')],'Color','k','Interpreter','Latex')

%------------
subplot(3,2,6)
h=cdfplot(dataLaserOn);
h.Color= colorA{2};
h.LineWidth=2;
h.LineStyle = '-';
hold on
plot(xx,cdf(pdL1_On,xx),':','color',colorB{2},'Linewidth',2) % for gamma data
%
plot(xx,cdf(pdL2_On,xx),'--','color', 'k','Linewidth',1) % for gamma model
title(['DATA/MODEL Comparison - D2cre-Opto: ' char(trialType(2))])
legend('Empirical','\bf Gamma fit', '\bf Model Fit',...
    'Location','northwest','Interpreter','Latex')
box off
legend boxoff
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20])
ylim([-0.05 1.05])
xlabel(xxlabel)
ylabel('CDF')
text(12, 0.60, ['Data: $\alpha$=' num2str(pdL1_On.a,'%1.2f') ...
    ', $\beta$=' num2str(1/pdL1_On.b,'%1.2f')],'Color',colorB{2},'Interpreter','Latex')
text(12, 0.45, ['Model: $\alpha$=' num2str(pdL2_On.a,'%1.2f') ...
    ', $\beta$=' num2str(1/pdL2_On.b,'%1.2f')],'Interpreter','Latex')
% ----------
text(12, 0.25, ['Data mean/CV: ' num2str(m_LaserOn ,'%1.2f') '; '  ...
    num2str(cv_LaserOn,'%1.2f')],'Color',colorB{2},'Interpreter','Latex')
text(12, 0.10, ['Model mean/CV: ' num2str(pdL2_On.a*pdL2_On.b ,'%1.2f') '; '  ...
    num2str(1/sqrt(pdL2_On.a),'%1.2f')],'Color','k','Interpreter','Latex')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------- plot results for model dynamics --------------------- %%
 set(0,'DefaultFigureWindowStyle','docked');
figure(5); 
%------------
subplot(2,4,3)
h=cdfplot(dataLaserOff);
colormap = [0 0 0; .9 .1 .1; .1 .1 .75; .5 .5 .5; .4 .4 .4];
h.Color= colormap(2,:);
h.LineWidth=2;
h.LineStyle = '-';
hold on
plot(xx,cdf(pdL1_Off,xx),':','color','r','Linewidth',2) % for gamma data
%
plot(xx,cdf(pdL2_Off,xx),'--','color', 'r','Linewidth',1) % for gamma model
title([char(trialType(1))])
box off
grid off
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20])
ylim([-0.05 1.05])
ylabel('CDF')
xlim([0 20]); set(gca, 'xtick', [0 6 18]); xlabel("Time (In Seconds)");

disp('D2-Off')
disp(['Data: Alpha' num2str(pdL1_Off.a), 'Beta' num2str(1/pdL1_Off.b)])
disp(['Model: Alpha' num2str(pdL2_Off.a), 'Beta' num2str(1/pdL2_Off.b)])
[dataCDF] = cdfcalc(dataLaserOff); xx1 = linspace(0,1,size(dataCDF,1));
gammaCDF = cdf(pdL1_Off,xx1); modelCDF = cdf(pdL2_Off,xx1);
%dataVsGamma = fitlm(dataCDF, gammaCDF); dataVsGamma.Rsquared.Adjusted
dataVsModel = fitlm(dataCDF, modelCDF); 
disp(['D2 MSN  data vs Model  ' num2str(dataVsModel.Rsquared.Adjusted)])
% GammaVsModel = fitlm(gammaCDF, modelCDF)
% GammaVsModel.Rsquared.Adjusted

%------------
subplot(2,4,4)
h=cdfplot(dataLaserOn);
h.Color= colormap(2,:);
h.LineWidth=2;
h.LineStyle = '-';
hold on
plot(xx,cdf(pdL1_On,xx),':','color','r','Linewidth',2) % for gamma data
%
plot(xx,cdf(pdL2_On,xx),'--','color', 'r','Linewidth',1) % for gamma model
title([ char(trialType(2))])
legend('Empirical','\bf Gamma fit', '\bf Model Fit',...
    'Location','northwest','Interpreter','Latex')
box off
grid off
legend boxoff
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20]); set(gca, 'xtick', [0 6 18]); xlabel("Time (In Seconds)");
ylim([-0.05 1.05])
ylabel('CDF')

disp('D2-On')
disp(['Data: Alpha' num2str(pdL1_On.a), 'Beta' num2str(1/pdL1_On.b)])
disp(['Model: Alpha' num2str(pdL2_On.a), 'Beta' num2str(1/pdL2_On.b)])
[dataCDF] = cdfcalc(dataLaserOn); xx1 = linspace(0,1,size(dataCDF,1));
gammaCDF = cdf(pdL1_On,xx1); modelCDF = cdf(pdL2_On,xx1);
%dataVsGamma = fitlm(dataCDF, gammaCDF); dataVsGamma.Rsquared.Adjusted
dataVsModel = fitlm(dataCDF, modelCDF); 
disp(['D2 MSN Disrupted data vs Model  ' num2str(dataVsModel.Rsquared.Adjusted)])
%GammaVsModel = fitlm(gammaCDF, modelCDF)
%GammaVsModel.Rsquared.Adjusted

%------------
subplot(2,4,1)
TT = Ymat_Off(samp_path_Off,:);
TT(isnan(TT)) = FR_target_Off;
y_min = min(min(TT));
y_max = FR_target_Off;

yLIMS = [y_min y_max] + 0.1*[-1 1] .* (y_max - y_min);
%
t= t0_Off:h_Off:tf_Off;
%
plot(t, Ymat_Off(samp_path_Off,:), 'LineWidth', 2, 'Color', [0.4 0.1 0.1])
hold on
line([t0_Off tf_Off], [baseline_Off baseline_Off], 'LineWidth', 2, 'Color', 'k', ...
    'LineStyle', ':')
text(14.5, 0.95*baseline_Off, ['baseline $b$=' num2str(baseline_Off,'%1.2f')],...
    'Interpreter','LaTex','FontSize', 12)
line([t0_Off tf_Off], [FR_target_Off FR_target_Off], 'LineWidth', 2, 'Color', 'k', ...
    'LineStyle', '--')
text(15, 0.97*FR_target_Off,['target $T$=' num2str(FR_target_Off,'%1.2f') ],...
    'Interpreter','LaTex','FontSize', 12)
text(9, 0.7,['Drift rate $D$=' num2str(D_Off,'%1.4f')],'Interpreter','LaTex','FontSize', 12)
text(9, 0.66,['Noise $\sigma$=' num2str(sigma_Off,'%1.4f')],'Interpreter','LaTex','FontSize', 12)
ylim(yLIMS)
xlim([0 20]); set(gca, 'xtick', [0 6 18]); xlabel("Time (In Seconds)");
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
box on
xlim([0 20]); set(gca, 'xtick', [0 6 18]); xlabel("Time (In Seconds)"); ylabel('$r$, Firing Rate',  'Interpreter', 'LaTex')
% title([trialType ': $dr = (1-r)\,D\, dt + \sigma \, d \, \xi(t)$. If $r=T$ then reset to $r=b$. (Here $\sigma$=' num2str(sigma_Off,'%1.4f') ', $D$=' num2str(D_Off,'%1.4f') ')'], ...
%     'Interpreter', 'latex', 'FontSize', 18)
title([ char(trialType(1))])
box off
clear TT y_min y_max yLIMS t

       
%------------
subplot(2,4,2)
TT = Ymat_On(samp_path_On,:);
TT(isnan(TT)) = FR_target_On;
y_min = min(min(TT));
y_max = FR_target_On;

yLIMS = [y_min y_max] + 0.1*[-1 1] .* (y_max - y_min);
%
t= t0_On:h_On:tf_On;
%
plot(t, Ymat_On(samp_path_On,:), 'LineWidth', 2, 'Color', [0.9 0.1 0.1])
hold on
line([t0_On tf_On], [baseline_On baseline_On], 'LineWidth', 2, 'Color', 'k', ...
    'LineStyle', ':')
text(14.5, 0.95*baseline_On, ['baseline $b$=' num2str(baseline_On,'%1.2f')],...
    'Interpreter','LaTex','FontSize', 12)
line([t0_On tf_On], [FR_target_On FR_target_On], 'LineWidth', 2, 'Color', 'k', ...
    'LineStyle', '--')
text(15, 0.97*FR_target_On,['target $T$=' num2str(FR_target_On,'%1.2f') ],...
    'Interpreter','LaTex','FontSize', 12)
text(10, 0.7,['Drift rate $D$=' num2str(D_On,'%1.4f')],'Interpreter','LaTex','FontSize', 12)
text(10, 0.66,['Noise $\sigma$=' num2str(sigma_On,'%1.4f')],'Interpreter','LaTex','FontSize', 12)
ylim(yLIMS)
xlim([0 20]); set(gca, 'xtick', [0 6 18]); xlabel("Time (In Seconds)");
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
box on
xlim([0 20]); set(gca, 'xtick', [0 6 18]); xlabel("Time (In Seconds)"); ylabel('$r$, Firing Rate',  'Interpreter', 'LaTex')
% title([trialType ': $dr = (1-r)\,D\, dt + \sigma \, d \, \xi(t)$. If $r=T$ then reset to $r=b$. (Here $\sigma$=' num2str(sigma_Off,'%1.4f') ', $D$=' num2str(D_Off,'%1.4f') ')'], ...
%     'Interpreter', 'latex', 'FontSize', 18)
title([char(trialType(2))])
box off
clear TT y_min y_max yLIMS t


%% ---------------- keep only the D1-cre mice ---------------------------
clear
T = readtable('OptoBehavioralTable.csv');
D1creMice = (strcmp(T.GENOTYPE,'OptoD1Halo'));
T = T(D1creMice,:);

% if you want to exclude any animal from the analysis, add their ID here
animalsToExclude = [];   % e.g.,   animalsToExclude = {'IOR3','SL10'}; 
%
trialType = {'LaserOff', 'LaserOn'};
colorA={'c','m'};
colorB={'b','r'};

%% calculations
for j=1:length(trialType)
    if strcmp(trialType(j),'LaserOff')
        LaserStatus = 'ctrl';  
        [dataLaserOff,m_LaserOff,md_LaserOff,cv_LaserOff] = AnimalData(T,LaserStatus);
        load(['D1_MSNs_LaserOff_ModelExample1.mat']);
        %
        TrialDur_Off =TrialDur; Ymat_Off =Ymat;  samp_path_Off = samp_path;
        baseline_Off = baseline; FR_target_Off = FR_target;  D_Off = D;
        sigma_Off = sigma; t0_Off = t0; tf_Off = tf; h_Off = h;
        clear Ymat samp_path baseline FR_target D sigma t0 tf h
       %-----------------------------------------------
    elseif strcmp(trialType(j),'LaserOn')
        LaserStatus = 'inhibition';    
        [dataLaserOn,m_LaserOn,md_LaserOn, cv_LaserOn] = AnimalData(T,LaserStatus);
        load(['D1_MSNs_LaserOn_ModelExample1.mat']);
        %
        TrialDur_On =TrialDur; Ymat_On =Ymat;  samp_path_On = samp_path;
        baseline_On = baseline; FR_target_On = FR_target;  D_On = D;
        sigma_On = sigma; t0_On = t0; tf_On = tf; h_On = h;
        clear Ymat samp_path baseline FR_target D sigma t0 tf h
        % --------------------------------------------
    else
        error('Wrong D2-cre-OptoInhib TrialType')
    end
end


%% ---- test for goodness of  model data with behavioral data -------------
%
% two-sample Kolmogorov-Smirnov test
% h = kstest2(x1,x2) returns a test decision for the null hypothesis that
% the data in vectors x1 and x2 are from the same continuous distribution,
% using the two-sample Kolmogorov-Smirnov test. The alternative hypothesis
% is that x1 and x2 are from different continuous distributions. The result
% h is 1 if the test rejects the null hypothesis at the 5% significance 
% level, and 0 otherwise.
[hvalOff,pvalOff] = kstest2(dataLaserOff,TrialDur_Off,'Alpha',0.05);
disp(['D1 Cre h for ' char(trialType(1)) ' = ' num2str(hvalOff)])
disp(['p for ' char(trialType(1)) ' = ' num2str(pvalOff)])

[hvalOn,pvalOn] = kstest2(dataLaserOn,TrialDur_On,'Alpha',0.05);
disp(['D1 MSN Disrupted h for ' char(trialType(2)) ' = ' num2str(hvalOn)])
disp(['p for ' char(trialType(2)) ' = ' num2str(pvalOn)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------------plot results for behavior ----------------------- %%
xx=linspace(0,20);
xxlabel = 'durations';
ymax=0.15;
nbins=30;

successfulTrials_Off = size(TrialDur_Off(~isnan(TrialDur_Off)),1); 
successfulTrials_On = size(TrialDur_On(~isnan(TrialDur_On)),1); 
figure(52); 
%----------------
ax1=subplot(3,2,1); 
h=histogram(ax1,dataLaserOff,nbins); 
h.FaceColor = colorA{1};
h.Normalization = 'pdf';
hold on
%
% try fitting with Gamma distrib
pdL1_Off = fitdist(dataLaserOff,'Gamma');
% ---- in shape-rate form -----------------------------------------------
aa=pdL1_Off.a; bb = 1/pdL1_Off.b;
zz =  (bb^aa/ gamma(aa)) .* xx.^(aa-1) .* exp(-bb*xx) ; 
%-------------------------------------------------------------------------
plot(xx,zz,':','color',colorB{1},'Linewidth',2.5)
%
title(['DATA:  D1cre-Opto (N=' num2str(length(dataLaserOff)) ', '...
    num2str(nbins) ' bins): ' char(trialType(1))])
legend('Empirical', '\bf Gamma fit',...
    'Location','northeast','Interpreter','Latex')
legend boxoff
box off
box off
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20])
ylim([0 ymax])
%xlabel(xxlabel)
ylabel('PDF')
text(15.5, 0.07, ['\bf shape $\alpha$=' num2str(pdL1_Off.a,'%1.2f')],'Interpreter','Latex')
text(15.5,0.05,['\bf rate $\beta$=' num2str(1/pdL1_Off.b,'%1.2f')],'Interpreter','Latex')
clear mu sigma h yy zz aa bb 

%-------------%
ax2=subplot(3,2,3); 
h=histogram(ax2,TrialDur_Off,nbins); 
h.FaceColor = [0.7 0.7 0.7] ;
h.Normalization = 'pdf';
hold on
%
% try fitting with Gamma distrib
pdL2_Off = fitdist(TrialDur_Off,'Gamma');
% ---- in shape-rate form -----------------------------------------------
aa=pdL2_Off.a; bb = 1/pdL2_Off.b;
zz =  (bb^aa/ gamma(aa)) .* xx.^(aa-1) .* exp(-bb*xx) ; 
%-------------------------------------------------------------------------
plot(xx,zz,'--k','Linewidth',2)
%
title(['MODEL SIMULATIONS (N=' num2str(successfulTrials_Off) ', '...
    num2str(nbins) ' bins): ' char(trialType(1))])
legend('Empirical','\bf Gamma fit','Location','northeast','Interpreter','Latex')
legend boxoff
box off
box off
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20])
ylim([0 ymax])
%xlabel(xxlabel)
ylabel('PDF')
text(15.5, 0.07, ['\bf shape $\alpha$=' num2str(pdL2_Off.a,'%1.2f')],'Interpreter','Latex')
text(15.5,0.05,['\bf rate $\beta$=' num2str(1/pdL2_Off.b,'%1.2f')],'Interpreter','Latex')
clear mu sigma h yy zz  aa bb

%----------------
ax1=subplot(3,2,2); 
h=histogram(ax1,dataLaserOn,nbins); 
h.FaceColor = colorA{2};
h.Normalization = 'pdf';
hold on
%
% try fitting with Gamma distrib
pdL1_On = fitdist(dataLaserOn,'Gamma');
% ---- in shape-rate form -----------------------------------------------
aa=pdL1_On.a; bb = 1/pdL1_On.b;
zz =  (bb^aa/ gamma(aa)) .* xx.^(aa-1) .* exp(-bb*xx) ; 
%-------------------------------------------------------------------------
plot(xx,zz,':','color',colorB{2},'Linewidth',2.5)
%
title(['DATA:  D1cre-Opto (N=' num2str(length(dataLaserOn)) ', '...
    num2str(nbins) ' bins): ' char(trialType(2))])
legend('Empirical', '\bf Gamma fit',...
    'Location','northeast','Interpreter','Latex')
legend boxoff
box off
box off
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20])
ylim([0 ymax])
%xlabel(xxlabel)
ylabel('PDF')
text(1, 0.13, ['\bf shape $\alpha$=' num2str(pdL1_On.a,'%1.2f')],'Interpreter','Latex')
text(1, 0.11,['\bf rate $\beta$=' num2str(1/pdL1_On.b,'%1.2f')],'Interpreter','Latex')
clear mu sigma h yy zz aa bb 

%-------------%
ax2=subplot(3,2,4); 
h=histogram(ax2,TrialDur_On,nbins); 
h.FaceColor = [0.7 0.7 0.7] ;
h.Normalization = 'pdf';
hold on
%
% try fitting with Gamma distrib
pdL2_On = fitdist(TrialDur_On,'Gamma');
% ---- in shape-rate form -----------------------------------------------
aa=pdL2_On.a; bb = 1/pdL2_On.b;
zz =  (bb^aa/ gamma(aa)) .* xx.^(aa-1) .* exp(-bb*xx) ; 
%-------------------------------------------------------------------------
plot(xx,zz,'--k','Linewidth',2)
%
title(['MODEL SIMULATIONS (N=' num2str(successfulTrials_On) ', '...
    num2str(nbins) ' bins): ' char(trialType(2))])
legend('Empirical','\bf Gamma fit','Location','northeast','Interpreter','Latex')
legend boxoff
box off
box off
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20])
ylim([0 ymax])
%xlabel(xxlabel)
ylabel('PDF')
text(1,0.13, ['\bf shape $\alpha$=' num2str(pdL2_On.a,'%1.2f')],'Interpreter','Latex')
text(1,0.11,['\bf rate $\beta$=' num2str(1/pdL2_On.b,'%1.2f')],'Interpreter','Latex')
clear mu sigma h yy zz  aa bb

%------------
subplot(3,2,5)
h=cdfplot(dataLaserOff);
h.Color= colorA{1};
h.LineWidth=2;
h.LineStyle = '-';
hold on
plot(xx,cdf(pdL1_Off,xx),':','color',colorB{1},'Linewidth',2) % for gamma data
%
plot(xx,cdf(pdL2_Off,xx),'--','color', 'k','Linewidth',1) % for gamma model
title(['DATA/MODEL Comparison - D1cre-Opto: ' char(trialType(1))])
legend('Empirical','\bf Gamma fit', '\bf Model Fit',...
    'Location','northwest','Interpreter','Latex')
box off
legend boxoff
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20])
ylim([-0.05 1.05])
xlabel(xxlabel)
ylabel('CDF')
text(12, 0.60, ['Data: $\alpha$=' num2str(pdL1_Off.a,'%1.2f') ...
    ', $\beta$=' num2str(1/pdL1_Off.b,'%1.2f')],'Color',colorB{1},'Interpreter','Latex')
text(12, 0.45, ['Model: $\alpha$=' num2str(pdL2_Off.a,'%1.2f') ...
    ', $\beta$=' num2str(1/pdL2_Off.b,'%1.2f')],'Interpreter','Latex')
% ----------
text(12, 0.25, ['Data mean/CV: ' num2str(m_LaserOff ,'%1.2f') '; '  ...
    num2str(cv_LaserOff,'%1.2f')],'Color',colorB{1},'Interpreter','Latex')
text(12, 0.10, ['Model mean/CV: ' num2str(pdL2_Off.a*pdL2_Off.b ,'%1.2f') '; '  ...
    num2str(1/sqrt(pdL2_Off.a),'%1.2f')],'Color','k','Interpreter','Latex')

%------------
subplot(3,2,6)
h=cdfplot(dataLaserOn);
h.Color= colorA{2};
h.LineWidth=2;
h.LineStyle = '-';
hold on
plot(xx,cdf(pdL1_On,xx),':','color',colorB{2},'Linewidth',2) % for gamma data
%
plot(xx,cdf(pdL2_On,xx),'--','color', 'k','Linewidth',1) % for gamma model
title(['DATA/MODEL Comparison - D1cre-Opto: ' char(trialType(2))])
legend('Empirical','\bf Gamma fit', '\bf Model Fit',...
    'Location','northwest','Interpreter','Latex')
box off
legend boxoff
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20])
ylim([-0.05 1.05])
xlabel(xxlabel)
ylabel('CDF')
text(13, 0.60, ['Data: $\alpha$=' num2str(pdL1_On.a,'%1.2f') ...
    ', $\beta$=' num2str(1/pdL1_On.b,'%1.2f')],'Color',colorB{2},'Interpreter','Latex')
text(13, 0.45, ['Model: $\alpha$=' num2str(pdL2_On.a,'%1.2f') ...
    ', $\beta$=' num2str(1/pdL2_On.b,'%1.2f')],'Interpreter','Latex')
% ----------
text(12, 0.25, ['Data mean/CV: ' num2str(m_LaserOn ,'%1.2f') '; '  ...
    num2str(cv_LaserOn,'%1.2f')],'Color',colorB{2},'Interpreter','Latex')
text(12, 0.10, ['Model mean/CV: ' num2str(pdL2_On.a*pdL2_On.b ,'%1.2f') '; '  ...
    num2str(1/sqrt(pdL2_On.a),'%1.2f')],'Color','k','Interpreter','Latex')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------- plot results for model dynamics --------------------- %%

figure(5)
%------------
subplot(2,4,7)
h=cdfplot(dataLaserOff);
h.Color= [0.1 0.1 0.9];
h.LineWidth=2;
h.LineStyle = '-';
hold on
plot(xx,cdf(pdL1_Off,xx),':','color','b','Linewidth',2) % for gamma data
%
plot(xx,cdf(pdL2_Off,xx),'--','color', 'k','Linewidth',1) % for gamma model
title([char(trialType(1))])
box off
legend boxoff
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)

ylim([-0.05 1.05])
xlabel(xxlabel)
ylabel('CDF')
xlim([0 20]); set(gca, 'xtick', [0 6 18]); xlabel("Time (In Seconds)");
grid off; 
disp('D1-Off')
disp(['Data: Alpha' num2str(pdL1_On.a), 'Beta' num2str(1/pdL1_Off.b)])
disp(['Model: Alpha' num2str(pdL2_On.a), 'Beta' num2str(1/pdL2_Off.b)])

[dataCDF] = cdfcalc(dataLaserOff); xx1 = linspace(0,1,size(dataCDF,1));
gammaCDF = cdf(pdL1_Off,xx1); modelCDF = cdf(pdL2_Off,xx1);
dataVsGamma = fitlm(dataCDF, gammaCDF); dataVsGamma.Rsquared.Adjusted
dataVsModel = fitlm(dataCDF, modelCDF); 
disp(['D1 MSN  data vs Model  ' num2str(dataVsModel.Rsquared.Adjusted)])
%GammaVsModel = fitlm(gammaCDF, modelCDF)
%GammaVsModel.Rsquared.Adjusted

%------------
subplot(2,4,8)
h=cdfplot(dataLaserOn);
h.Color= [0.1 0.1 0.9];
h.LineWidth=2;
h.LineStyle = '-';
hold on
plot(xx,cdf(pdL1_On,xx),':','color','b','Linewidth',2) % for gamma data
%
plot(xx,cdf(pdL2_On,xx),'--','color', 'k','Linewidth',1) % for gamma model
title([char(trialType(2))])
box off
legend boxoff
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20])
ylim([-0.05 1.05])
xlabel(xxlabel)
ylabel('CDF')
xlim([0 20]); set(gca, 'xtick', [0 6 18]); xlabel("Time (In Seconds)");
grid off; 
disp('D1 On')
disp(['Data: Alpha' num2str(pdL1_On.a), 'Beta' num2str(1/pdL1_On.b)])
disp(['Model: Alpha' num2str(pdL2_On.a), 'Beta' num2str(1/pdL2_On.b)])

[dataCDF] = cdfcalc(dataLaserOn); xx1 = linspace(0,1,size(dataCDF,1));
gammaCDF = cdf(pdL1_On,xx1); modelCDF = cdf(pdL2_On,xx1);
% dataVsGamma = fitlm(dataCDF, gammaCDF); dataVsGamma.Rsquared.Adjusted
dataVsModel = fitlm(dataCDF, modelCDF); 
disp(['D1 MSN Disrupted data vs Model  ' num2str(dataVsModel.Rsquared.Adjusted)])
GammaVsModel = fitlm(gammaCDF, modelCDF); GammaVsModel.Rsquared.Adjusted

%------------
subplot(2,4,5)
TT = Ymat_Off(samp_path_Off,:);
TT(isnan(TT)) = FR_target_Off;
y_max = max(max(TT));

y_min = FR_target_Off;

yLIMS = [y_min y_max] + 0.1*[-1 1] .* (y_max - y_min);
yLIMS = [0 yLIMS(2)];
%
t= t0_Off:h_Off:tf_Off;
%
plot(t, Ymat_Off(samp_path_Off,:), 'LineWidth', 2, 'Color', [.1 .1 .9]);
hold on
line([t0_Off tf_Off], [baseline_Off baseline_Off], 'LineWidth', 2, 'Color', 'k', ...
    'LineStyle', ':')
text(14.5, 0.9*baseline_Off, ['baseline $b$=' num2str(baseline_Off,'%1.2f')],...
    'Interpreter','LaTex','FontSize', 12)
line([t0_Off tf_Off], [FR_target_Off FR_target_Off], 'LineWidth', 2, 'Color', 'k', ...
    'LineStyle', '--')
text(15, 0.65*FR_target_Off,['target $T$=' num2str(FR_target_Off,'%1.2f') ],...
    'Interpreter','LaTex','FontSize', 12)
text(6, 0.38,['Drift rate $D$=' num2str(D_Off,'%1.4f')],'Interpreter','LaTex','FontSize', 12)
text(6, 0.34,['Noise $\sigma$=' num2str(sigma_Off,'%1.4f')],'Interpreter','LaTex','FontSize', 12)
ylim(yLIMS)
xlim([0 20])
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
box on
xlim([0 20]); set(gca, 'xtick', [0 6 18]); xlabel("Time (In Seconds)");
ylabel('$r$, Firing Rate',  'Interpreter', 'LaTex')
% title([trialType ': $dr = (1-r)\,D\, dt + \sigma \, d \, \xi(t)$. If $r=T$ then reset to $r=b$. (Here $\sigma$=' num2str(sigma_Off,'%1.4f') ', $D$=' num2str(D_Off,'%1.4f') ')'], ...
%     'Interpreter', 'latex', 'FontSize', 18)
title(['MODEL Trajectories for D1cre-Opto: ' char(trialType(1))])
box off
clear TT y_min y_max yLIMS t

       
%------------
subplot(2,4,6)
TT = Ymat_On(samp_path_On,:);
TT(isnan(TT)) = FR_target_On;
y_max = max(max(TT));

y_min = FR_target_On;

yLIMS = [y_min y_max] + 0.1*[-1 1] .* (y_max - y_min);
yLIMS = [0 yLIMS(2)];
%
t= t0_On:h_On:tf_On;
%
plot(t, Ymat_On(samp_path_On,:), 'LineWidth', 2, 'Color', [.1 .1 .9]);
hold on
line([t0_On tf_On], [baseline_On baseline_On], 'LineWidth', 2, 'Color', 'k', ...
    'LineStyle', ':')
text(14.5, 0.9*baseline_On, ['baseline $b$=' num2str(baseline_On,'%1.2f')],...
    'Interpreter','LaTex','FontSize', 12)
line([t0_On tf_On], [FR_target_On FR_target_On], 'LineWidth', 2, 'Color', 'k', ...
    'LineStyle', '--')
text(15, 0.65*FR_target_On,['target $T$=' num2str(FR_target_On,'%1.2f') ],...
    'Interpreter','LaTex','FontSize', 12)
text(5, 0.44,['Drift rate $D$=' num2str(D_On,'%1.4f')],'Interpreter','LaTex','FontSize', 12)
text(5, 0.4,['Noise $\sigma$=' num2str(sigma_On,'%1.4f')],'Interpreter','LaTex','FontSize', 12)
ylim(yLIMS)
xlim([0 20])
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
box on
xlim([0 20]); set(gca, 'xtick', [0 6 18]); xlabel("Time (In Seconds)");ylabel('$r$, Firing Rate',  'Interpreter', 'LaTex')
% title([trialType ': $dr = (1-r)\,D\, dt + \sigma \, d \, \xi(t)$. If $r=T$ then reset to $r=b$. (Here $\sigma$=' num2str(sigma_Off,'%1.4f') ', $D$=' num2str(D_Off,'%1.4f') ')'], ...
%     'Interpreter', 'latex', 'FontSize', 18)
title(['MODEL Trajectories for D1cre-Opto: ' char(trialType(2))])
box off
clear TT y_min y_max yLIMS t



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------------------------------------------------------------%%
function [x,mx,mdx,cvx] = AnimalData(T,LaserStatus)
% 
x = nan(size(T,1),1);

%%%%% animal = unique(T.ANIMALID(strcmp(T.LASER,LaserStatus))) ;
selection = (strcmp(T.LASER,LaserStatus)) ;

if ~isempty(selection)
   x(1:sum(selection)) = T(selection,1).RSPTIME; 
end
x = x(~isnan(x));

mx=mean(x);
mdx=median(x);
cvx = std(x)./mean(x); 
%
end
%-------------------------------------------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%