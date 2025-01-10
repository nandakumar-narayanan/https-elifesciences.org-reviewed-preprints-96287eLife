%% Written by Rodica Curtu / checked by Narayanna
% Curtu / Narayanna Lab

%folderModelFiles = './DataFiles_Model/';

%% ---------------  read data from file -----------------------------------
T = readtable('OptoBehavioralTableV3.csv');

% this generates a table with 4 columns, in this order:
% T.RSPTIME=Response Time for animal, seconds
% T.LASER = Laser (Whether the laser is inhibiting D1/D2 MSNs or not)
% T.ANIMALID = animal ID
% T.GENOTYPE = either D1-cre or D2-cre or control , i.e.
% here only: OptoD1Halo = D1-cre,  OptoD2Halo = D2-cre, 
%            OptoD2CtrlHalo = Control mice but perturbations to D2 MSNs (?)

%% ---------------- keep only the D2-cre mice ---------------------------
D2creMice = (strcmp(T.GENOTYPE,'OptoD2Halo'));
T = T(D2creMice,:);

% if you want to exclude any animal from the analysis, add their ID here
% % % %animalsToExclude = [];   % e.g.,   animalsToExclude = {'IOR3','SL10'}; 
%
trialType = {'LaserOff', 'LaserOn'};
colorB={'b','r'};

%% calculations
for j=1:length(trialType)
    if strcmp(trialType(j),'LaserOff')
        LaserStatus = 'ctrl';  
        [dataLaserOff,m_LaserOff,md_LaserOff,cv_LaserOff] = AnimalData(T,LaserStatus);
        load('D2_MSNs_LaserOff_ModelExample1_updated.mat');
        %
        TrialDur_Off =TrialDur; Ymat_Off =Ymat;  samp_path_Off = samp_path;
        baseline_Off = baseline; FR_target_Off = FR_target;  D_Off = D;
        sigma_Off = sigma; t0_Off = t0; tf_Off = tf; h_Off = h;
        clear Ymat samp_path baseline FR_target D sigma t0 tf h
       %-----------------------------------------------
    elseif strcmp(trialType(j),'LaserOn')
        LaserStatus = 'inhibition';    
        [dataLaserOn,m_LaserOn,md_LaserOn, cv_LaserOn] = AnimalData(T,LaserStatus);
        load('D2_MSNs_LaserOn_ModelExample1_updated.mat');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------------plot results for behavior ----------------------- %%
xx=linspace(0,20);
xxlabel = 'Time (seconds)';
ymax=0.15;
nbins=30;

set(0,'DefaultFigureWindowStyle','docked');
figure(41)
%----------------
ax1=subplot(3,4,1); 
h=histogram(ax1,dataLaserOff,nbins); 
h.FaceColor = [0.7 0.7 0.7];
h.Normalization = 'pdf';
hold on
%
% try fitting with Gamma distrib
pdL1_Off = fitdist(dataLaserOff,'Gamma');
% ---- in shape-rate form -----------------------------------------------
aa=pdL1_Off.a; bb = 1/pdL1_Off.b;
zz =  (bb^aa/ gamma(aa)) .* xx.^(aa-1) .* exp(-bb*xx) ; 
%-------------------------------------------------------------------------
plot(xx,zz,':','color',[0.7 0.7 0.7],'Linewidth',2.5)
%
title('D2-Cre mice: Laser Off')

box off; legend off
box off; legend off
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");
ylim([0 ymax])
ylabel('PDF')
text(15.5, 0.07, ['\bf shape $\alpha$=' num2str(pdL1_Off.a,'%1.2f')],'Interpreter','Latex')
text(15.5,0.05,['\bf rate $\beta$=' num2str(1/pdL1_Off.b,'%1.2f')],'Interpreter','Latex')
clear mu sigma h yy zz aa bb 

%-------------%
ax2=subplot(3,4,5); cla; 
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
title('MODEL')
box off; legend off
box off; legend off
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");
ylim([0 ymax])
ylabel('PDF')
text(15.5, 0.07, ['\bf shape $\alpha$=' num2str(pdL2_Off.a,'%1.2f')],'Interpreter','Latex')
text(15.5,0.05,['\bf rate $\beta$=' num2str(1/pdL2_Off.b,'%1.2f')],'Interpreter','Latex')
clear mu sigma h yy zz  aa bb

%----------------
ax1=subplot(3,4,2); 
h=histogram(ax1,dataLaserOn,nbins); 
h.FaceColor = [1 0 0];
h.Normalization = 'pdf';
hold on
%
% try fitting with Gamma distrib
pdL1_On = fitdist(dataLaserOn,'Gamma');
% ---- in shape-rate form -----------------------------------------------
aa=pdL1_On.a; bb = 1/pdL1_On.b;
zz =  (bb^aa/ gamma(aa)) .* xx.^(aa-1) .* exp(-bb*xx) ; 
%-------------------------------------------------------------------------
plot(xx,zz,':','color',[1 0 0],'Linewidth',2.5)
%
title('DATA: D2Cre Laser ON')

box off; legend off
box off; legend off
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");
ylim([0 ymax])
%xlabel(xxlabel)
ylabel('PDF')
text(15.5, 0.07, ['\bf shape $\alpha$=' num2str(pdL1_On.a,'%1.2f')],'Interpreter','Latex')
text(15.5, 0.05,['\bf rate $\beta$=' num2str(1/pdL1_On.b,'%1.2f')],'Interpreter','Latex')
clear mu sigma h yy zz aa bb 

%-------------%
ax2=subplot(3,4,6); cla; 
h=histogram(ax2,TrialDur_On,nbins); 
h.FaceColor = [1 0 0 ] ;
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
title('Model Laser ON')
legend off 
box off; legend off
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");
ylim([0 ymax])
ylabel('PDF')
text(15.5, 0.07, ['\bf shape $\alpha$=' num2str(pdL2_On.a,'%1.2f')],'Interpreter','Latex')
text(15.5,0.05,['\bf rate $\beta$=' num2str(1/pdL2_On.b,'%1.2f')],'Interpreter','Latex')
clear mu sigma h yy zz  aa bb
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)"); set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");

%------------
subplot(3,4,9); cla; 
h=cdfplot(dataLaserOff);
h.Color= [0.7 0.7 0.7];
h.LineWidth=2;
h.LineStyle = '-';
hold on
plot(xx,cdf(pdL1_Off,xx),':','color',[0.7 0.7 0.7],'Linewidth',2) % for gamma data
%
plot(xx,cdf(pdL2_Off,xx),'--','color', 'k','Linewidth',2) % for gamma model
title('DATA/MODEL  D2cre-Opto: ')
box off; legend off
grid off; 

set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");
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
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");


%------------
subplot(3,4,10)
h=cdfplot(dataLaserOn);
h.Color= [1 0 0];
h.LineWidth=2;
h.LineStyle = '-';
hold on
plot(xx,cdf(pdL1_On,xx),':','color',[1 0 0],'Linewidth',2) % for gamma data
%
plot(xx,cdf(pdL2_On,xx),'--','color', 'k','Linewidth',2) % for gamma model
title('DATA/MODEL Laser ON')
box off; legend off

set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");
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
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");
grid off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------- plot results for model dynamics --------------------- %%
 set(0,'DefaultFigureWindowStyle','docked');
figure(4); 
%------------
subplot(2,3,3); cla; 

plot(xx,cdf(pdL2_Off,xx),'--','color', 'k','Linewidth',1) % for gamma model
%title([char(trialType(1))])
box off; legend off
grid off
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");
ylim([-0.05 1.05])
ylabel('CDF')
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)"); set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------- Write CSV Files --------------------- %%

fprintf('\n\n>>>>>D2-Off')
[dataCDF] = cdfcalc(dataLaserOff);
xx1 = linspace(0,20,size(dataCDF,1));  % work on 20 s time interval (not 1s)
modelCDF = cdf(pdL2_Off,xx1)';

ThresholdOFF = (1-0.25*baseline_Off) * FR_target_Off + 0.25 * baseline_Off * (1-FR_target_Off);
writematrix([FR_target_Off,baseline_Off,ThresholdOFF, D_Off, ...
    sigma_Off ThresholdOFF], 'ModelD2Off.csv');

% Gamma distribution fit
gammaCDF = cdf(pdL1_Off,xx1)'; 
dataVsGamma = fitlm(dataCDF, gammaCDF); 
writematrix([pdL2_Off.a, 1/pdL2_Off.b dataVsGamma.Rsquared.Adjusted], ...
    'D2OFFGammaParameters.csv');

ModelVsGamma = fitlm(modelCDF, gammaCDF); 
dataVsModel = fitlm(dataCDF, modelCDF); 
writematrix([pdL2_Off.a  1/pdL2_Off.b ModelVsGamma.Rsquared.Adjusted], ...
    'D2Gamma.csv');
D2DataVsModelOFF = dataVsModel.Rsquared.Adjusted;

%
%------------
subplot(2,3,3)
hold on
% plot(xx,cdf(pdL1_On,xx),':','color','r','Linewidth',2) % for gamma data
%
plot(xx,cdf(pdL2_On,xx),'--','color', 'r','Linewidth',1) % for gamma model
%title([ char(trialType(2))])
box off; legend off
grid off

set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)"); set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");
ylim([-0.05 1.05])
ylabel('CDF')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------- Write CSV Files --------------------- %%
fprintf('\n\n>>>>>D2-On')
[dataCDF] = cdfcalc(dataLaserOn); 
xx1 = linspace(0,20,size(dataCDF,1));  % again,  work on 20 s time interval (not 1s)
modelCDF = cdf(pdL2_On,xx1)';
dataVsModel = fitlm(dataCDF, modelCDF); 
D2DataVsModelOn = dataVsModel.Rsquared.Adjusted;
writematrix([D2DataVsModelOFF D2DataVsModelOn],'D2ModelFit.csv'); 
ThresholdON = (1-0.25*baseline_On) * FR_target_On + 0.25 * baseline_On * (1-FR_target_On);
writematrix([FR_target_On,baseline_On, D_On, sigma_On], 'ModelD2On.csv');

% Gamma distribution fit
gammaCDF = cdf(pdL1_On,xx1)'; 
dataVsGamma = fitlm(dataCDF, gammaCDF); 
writematrix([pdL2_On.a, 1/pdL2_On.b dataVsGamma.Rsquared.Adjusted], ...
    'D2ONGammaParameters.csv');
% ModelVsGamma = fitlm(modelCDF, gammaCDF); 



%% ------------
subplot(2,3,1)
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
line([t0_Off tf_Off], [ThresholdOFF ThresholdOFF], 'LineWidth', 2, 'Color', 'k', ...
    'LineStyle', '--')
text(15, 0.97*FR_target_Off,['input-target $F$=' num2str(FR_target_Off,'%1.2f') ],...
    'Interpreter','LaTex','FontSize', 12)
text(15, 0.97*ThresholdOFF,['Threshold $T$=' num2str(ThresholdOFF,'%1.2f') ],...
    'Interpreter','LaTex','FontSize', 12)
text(9, 0.7,['Drift rate $D$=' num2str(D_Off,'%1.4f')],'Interpreter','LaTex','FontSize', 12)
text(9, 0.66,['Noise $\sigma$=' num2str(sigma_Off,'%1.4f')],'Interpreter','LaTex','FontSize', 12)
ylim(yLIMS)
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)"); set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
box on
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)"); set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)"); ylabel('$x$, Firing Rate',  'Interpreter', 'LaTex')
% % % % title([trialType ': $dx = (F-x)\,D\, dt + \sigma \, d \, \xi(t)$. If $r=thresh$ then reset to $x=b$. (Here $\sigma$=' num2str(sigma_Off,'%1.4f') ', $D$=' num2str(D_Off,'%1.4f') ')'], ...
% % % %     'Interpreter', 'latex', 'FontSize', 18)
title([ char(trialType(1))])
box off; legend off
clear TT y_min y_max yLIMS t

       
% ------------
subplot(2,3,2)
TT = Ymat_On(samp_path_On,:);
TT(isnan(TT)) = FR_target_On;
y_min = min(min(TT));
y_max = FR_target_On;

%
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
line([t0_On tf_On], [ThresholdON ThresholdON], 'LineWidth', 2, 'Color', 'k', ...
    'LineStyle', '--')
text(15, 0.97*FR_target_On,['input-target $F$=' num2str(FR_target_On,'%1.2f') ],...
    'Interpreter','LaTex','FontSize', 12)
text(15, 0.97*ThresholdON,['Threshold $T$=' num2str(ThresholdON,'%1.2f') ],...
    'Interpreter','LaTex','FontSize', 12)
text(10, 0.7,['Drift rate $D$=' num2str(D_On,'%1.4f')],'Interpreter','LaTex','FontSize', 12)
text(10, 0.66,['Noise $\sigma$=' num2str(sigma_On,'%1.4f')],'Interpreter','LaTex','FontSize', 12)
ylim(yLIMS)
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)"); set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
box on
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)"); set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)"); ylabel('$r$, Firing Rate',  'Interpreter', 'LaTex')
title([char(trialType(2))])
box off; legend off
clear TT y_min y_max yLIMS t


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------- D1 mice ---------------------------
clear
T = readtable('OptoBehavioralTableV3.csv');
D1creMice = (strcmp(T.GENOTYPE,'OptoD1Halo'));
T = T(D1creMice,:);

% if you want to exclude any animal from the analysis, add their ID here
% % % animalsToExclude = [];   % e.g.,   animalsToExclude = {'IOR3','SL10'}; 
%
trialType = {'LaserOff', 'LaserOn'};
% % % % colorA={'c','m'};
colorB={'b','r'};

%% calculations
for j=1:length(trialType)
    if strcmp(trialType(j),'LaserOff')
        LaserStatus = 'ctrl';  
        [dataLaserOff,m_LaserOff,md_LaserOff,cv_LaserOff] = AnimalData(T,LaserStatus);
        load('D1_MSNs_LaserOff_ModelExample1_updated.mat');
        %
        TrialDur_Off =TrialDur; Ymat_Off =Ymat;  samp_path_Off = samp_path;
        baseline_Off = baseline; FR_target_Off = FR_target;  D_Off = D;
        sigma_Off = sigma; t0_Off = t0; tf_Off = tf; h_Off = h;
        clear Ymat samp_path baseline FR_target D sigma t0 tf h
       %-----------------------------------------------
    elseif strcmp(trialType(j),'LaserOn')
        LaserStatus = 'inhibition';    
        [dataLaserOn,m_LaserOn,md_LaserOn, cv_LaserOn] = AnimalData(T,LaserStatus);
        load('D1_MSNs_LaserOn_ModelExample1_updated.mat');
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------------plot results for behavior ----------------------- %%
xx=linspace(0,20);
xxlabel = 'Time (seconds)';
ymax=0.15;
nbins=30;

figure(41); 
%----------------
ax1=subplot(3,4,3); 
h=histogram(ax1,dataLaserOff,nbins); 
h.FaceColor = [0.7 0.7 0.7] ;
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
title('DATA:  D1cre-mice Laser OFF')
box off; legend off
box off; legend off
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");
ylim([0 ymax])
ylabel('PDF')
text(15.5, 0.07, ['\bf shape $\alpha$=' num2str(pdL1_Off.a,'%1.2f')],'Interpreter','Latex')
text(15.5,0.05,['\bf rate $\beta$=' num2str(1/pdL1_Off.b,'%1.2f')],'Interpreter','Latex')
clear mu sigma h yy zz aa bb 

%-------------%
ax2=subplot(3,4,7); 
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
title('MODEL')

box off; legend off
box off; legend off
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");
ylim([0 ymax])
ylabel('PDF')
text(15.5, 0.07, ['\bf shape $\alpha$=' num2str(pdL2_Off.a,'%1.2f')],'Interpreter','Latex')
text(15.5,0.05,['\bf rate $\beta$=' num2str(1/pdL2_Off.b,'%1.2f')],'Interpreter','Latex')
clear mu sigma h yy zz  aa bb

%----------------
ax1=subplot(3,4,4); cla; 
h=histogram(ax1,dataLaserOn,nbins); 
h.FaceColor = [0 0 1] ;
h.Normalization = 'pdf';
hold on
%
% try fitting with Gamma distrib
pdL1_On = fitdist(dataLaserOn,'Gamma');
% ---- in shape-rate form -----------------------------------------------
aa=pdL1_On.a; bb = 1/pdL1_On.b;
zz =  (bb^aa/ gamma(aa)) .* xx.^(aa-1) .* exp(-bb*xx) ; 
%-------------------------------------------------------------------------
plot(xx,zz,':','color',[0 0 1],'Linewidth',2.5)
%
title('DATA:  D1cre Laser ON')

box off; legend off
box off; legend off
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");
ylim([0 ymax])
ylabel('PDF')
text(1, 0.13, ['\bf shape $\alpha$=' num2str(pdL1_On.a,'%1.2f')],'Interpreter','Latex')
text(1, 0.11,['\bf rate $\beta$=' num2str(1/pdL1_On.b,'%1.2f')],'Interpreter','Latex')
clear mu sigma h yy zz aa bb 

%-------------%
ax2=subplot(3,4,8); 
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
title('MODEL')
box off; legend off
box off; legend off
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");
ylim([0 ymax])
ylabel('PDF')
text(1,0.13, ['\bf shape $\alpha$=' num2str(pdL2_On.a,'%1.2f')],'Interpreter','Latex')
text(1,0.11,['\bf rate $\beta$=' num2str(1/pdL2_On.b,'%1.2f')],'Interpreter','Latex')
clear mu sigma h yy zz  aa bb

%------------
subplot(3,4,11)
h=cdfplot(dataLaserOff);
h.Color= [0.7 0.7 0.7];
h.LineWidth=2;
h.LineStyle = '-';
hold on
plot(xx,cdf(pdL1_Off,xx),':','color',[0.7 0.7 0.7],'Linewidth',2) % for gamma data
%
plot(xx,cdf(pdL2_Off,xx),'--','color', 'k','Linewidth',1) % for gamma model
title('DATA/MODEL Comparison')

box off; legend off

set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");
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
grid off; 
%------------
subplot(3,4,12)
h=cdfplot(dataLaserOn);
h.Color= [0 0 1];
h.LineWidth=2;
h.LineStyle = '-';
hold on
plot(xx,cdf(pdL1_On,xx),':','color',[0 0 1],'Linewidth',2) % for gamma data
%
plot(xx,cdf(pdL2_On,xx),'--','color', 'k','Linewidth',1) % for gamma model
title('DATA/MODEL Comparison')
box off; legend off

set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");
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
grid off; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------- plot results for model dynamics --------------------- %%

figure(4)
%------------
subplot(2,3,6)
hold on
%plot(xx,cdf(pdL1_Off,xx),':','color','k','Linewidth',2) % for gamma data
%
plot(xx,cdf(pdL2_Off,xx),'--','color', 'k','Linewidth',1) % for gamma model
%title([char(trialType(1))])
box off; legend off

set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)

ylim([-0.05 1.05])
xlabel(xxlabel)
ylabel('CDF')
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)"); set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");
grid off; 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------- Write CSV Files --------------------- %%

fprintf('\n\n>>>>>D1-Off')

[dataCDF] = cdfcalc(dataLaserOff);
xx1 = linspace(0,20,size(dataCDF,1));
modelCDF = cdf(pdL2_Off,xx1)';
% dataVsModel = fitlm(dataCDF, modelCDF); 
ThresholdOFF = (1-0.25*(baseline_Off)) * FR_target_Off + 0.25 * baseline_Off * (1-FR_target_Off);
writematrix([FR_target_Off,baseline_Off,ThresholdOFF, D_Off, sigma_Off], ...
    'ModelD1Off.csv');

gammaCDF = cdf(pdL2_Off,xx1)'; 
dataVsGamma = fitlm(dataCDF, gammaCDF); 
ModelVsGamma = fitlm(modelCDF, gammaCDF); 
writematrix([pdL2_Off.a, 1/pdL2_Off.b dataVsGamma.Rsquared.Adjusted], ...
    'D1OFFGammaParameters.csv');

dataVsModel = fitlm(dataCDF, modelCDF); 
writematrix([pdL2_Off.a  1/pdL2_Off.b ModelVsGamma.Rsquared.Adjusted],...
    'D1Gamma.csv');

D1DataVsModelOff = dataVsModel.Rsquared.Adjusted;



%------------
subplot(2,3,6)
hold on
%plot(xx,cdf(pdL1_On,xx),':','color','b','Linewidth',2) % for gamma data
%
plot(xx,cdf(pdL2_On,xx),'--','color', 'b','Linewidth',1) % for gamma model
%title([char(trialType(2))])
box off; legend off

set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");
ylim([-0.05 1.05])
xlabel(xxlabel)
ylabel('CDF')
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)"); set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");
grid off; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------- Write CSV Files --------------------- %%
fprintf('\n\n>>>>>D1-On')

[dataCDF] = cdfcalc(dataLaserOn); 

xx1 = linspace(0,20,size(dataCDF,1));
% gammaCDF = cdf(pdL1_On,xx1)'; 
modelCDF = cdf(pdL2_On,xx1)';

dataVsModel = fitlm(dataCDF, modelCDF); 
D1DataVsModelOn = dataVsModel.Rsquared.Adjusted;
writematrix([D1DataVsModelOff D1DataVsModelOn], 'D1ModelFit.csv'); 

ThresholdON = (1-0.25*baseline_On) * FR_target_On + 0.25 * baseline_On * (1-FR_target_On);
[dataCDF] = cdfcalc(dataLaserOn); 
gammaCDF = cdf(pdL1_On,xx1)'; 
modelCDF = cdf(pdL2_On,xx1)';
dataVsGamma = fitlm(dataCDF, gammaCDF); 
writematrix([pdL2_On.a, 1/pdL2_On.b dataVsGamma.Rsquared.Adjusted], ...
    'D1ONGammaParameters.csv');

ModelVsGamma = fitlm(modelCDF, gammaCDF); 


%------------
subplot(2,3,4)
TT = Ymat_Off(samp_path_Off,:);
TT(isnan(TT)) = FR_target_Off;
y_max = max(max(TT));

y_min = FR_target_Off;

yLIMS = [y_min y_max] + 0.1*[-1 1] .* (y_max - y_min);
yLIMS = [0 yLIMS(2)];
%
t= t0_Off:h_Off:tf_Off;
%
plot(t, Ymat_Off(samp_path_Off,:), 'LineWidth', 2, 'Color', [.1 .1 .7]);
hold on
line([t0_Off tf_Off], [baseline_Off baseline_Off], 'LineWidth', 2, 'Color', 'k', ...
    'LineStyle', ':')
text(14.5, 0.9*baseline_Off, ['baseline $b$=' num2str(baseline_Off,'%1.2f')],...
    'Interpreter','LaTex','FontSize', 12)
line([t0_Off tf_Off], [ThresholdOFF ThresholdOFF], 'LineWidth', 2, 'Color', 'k', ...
    'LineStyle', '--')
text(15, 0.65*FR_target_Off,['input-target $F$=' num2str(FR_target_Off,'%1.2f') ],...
    'Interpreter','LaTex','FontSize', 12)
text(15, 0.65*ThresholdOFF,['Threshold $T$=' num2str(ThresholdOFF,'%1.2f') ],...
    'Interpreter','LaTex','FontSize', 12)
text(6, 0.38,['Drift rate $D$=' num2str(D_Off,'%1.4f')],'Interpreter','LaTex','FontSize', 12)
text(6, 0.34,['Noise $\sigma$=' num2str(sigma_Off,'%1.4f')],'Interpreter','LaTex','FontSize', 12)
ylim(yLIMS)
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
box on
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)"); set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");
ylabel('$r$, Firing Rate',  'Interpreter', 'LaTex')
title(['MODEL Trajectories for D1cre-Opto: ' char(trialType(1))])
box off; legend off
clear TT y_min y_max yLIMS t

       
%------------
subplot(2,3,5)
TT = Ymat_On(samp_path_On,:);
TT(isnan(TT)) = FR_target_On;
y_max = max(max(TT));

y_min = FR_target_On;

yLIMS = [y_min y_max] + 0.1*[-1 1] .* (y_max - y_min);
yLIMS = [0 yLIMS(2)];
%
t= t0_On:h_On:tf_On;
%
%
plot(t, Ymat_On(samp_path_On,:), 'LineWidth', 2, 'Color', [.1 .1 .9]);
hold on
line([t0_On tf_On], [baseline_On baseline_On], 'LineWidth', 2, 'Color', 'k', ...
    'LineStyle', ':')
text(14.5, 0.9*baseline_On, ['baseline $b$=' num2str(baseline_On,'%1.2f')],...
    'Interpreter','LaTex','FontSize', 12)
line([t0_On tf_On], [ThresholdON ThresholdON], 'LineWidth', 2, 'Color', 'k', ...
    'LineStyle', '--')
text(15, 0.65*FR_target_On,['input-target $F$=' num2str(FR_target_On,'%1.2f') ],...
    'Interpreter','LaTex','FontSize', 12)
text(15, 0.65*ThresholdON,['Threshold $T$=' num2str(ThresholdON,'%1.2f') ],...
    'Interpreter','LaTex','FontSize', 12)
text(5, 0.44,['Drift rate $D$=' num2str(D_On,'%1.4f')],'Interpreter','LaTex','FontSize', 12)
text(5, 0.4,['Noise $\sigma$=' num2str(sigma_On,'%1.4f')],'Interpreter','LaTex','FontSize', 12)
ylim(yLIMS)
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1)
box on
xlim([0 20]);  set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)"); set(gca, 'xtick', [0 6 18]); xlabel("Time (seconds)");ylabel('$r$, Firing Rate',  'Interpreter', 'LaTex')

title(['MODEL Trajectories for D1cre-Opto: ' char(trialType(2))])
box off; legend off
clear TT y_min y_max yLIMS t

writematrix([baseline_On, FR_target_On, D_On, sigma_On],'ModelD1On.csv');


%% Print Results - the data are out of order from the results section because
%% variables are cleared.  This section is intended to match the results
fprintf('\n\n>>>>>>Model Results')
ModelD2Off = readmatrix('ModelD2Off.csv');
fprintf('\nD2 OFF Model parameters: FR Target  %0.3f baseline %0.3f Threshold %0.2f  D %0.3f sigma %0.3f', ModelD2Off(1), ModelD2Off(2), ModelD2Off(3),ModelD2Off(4), ModelD2Off(5))

ModelD1Off = readmatrix('ModelD1Off.csv');
fprintf('\nD1 OFF Model parameters: FR Target  %0.3f baseline %0.3f Threshold %0.2f  D %0.3f sigma %0.3f', ModelD1Off(1), ModelD1Off(2), ModelD1Off(3),ModelD1Off(4), ModelD1Off(5))

D2Gamma = readmatrix('D2Gamma.csv');
fprintf('\n D2 Gamma Alpha %2.2f Beta %2.2f %2.2f', D2Gamma(1), D2Gamma(2), D2Gamma(3));

D1Gamma = readmatrix('D1Gamma.csv');
fprintf('\n D1 Gamma Alpha %2.2f Beta %2.2f %2.2f', D1Gamma(1), D1Gamma(2), D1Gamma(3));

ModelD2On = readmatrix('ModelD2On.csv');
fprintf('\nD2 On D %0.3f sigma %0.3f', ModelD2On(3), ModelD2On(4));

ModelD1On = readmatrix('ModelD1On.csv');
fprintf('\nD1 On D %0.3f sigma %0.3f', ModelD1On(3), ModelD1On(4));

D1GammaParameters = readmatrix('D1OFFGammaParameters.csv');
D2GammaParameters = readmatrix('D2OFFGammaParameters.csv');
fprintf('\n Model fit to Gamma for D2 Laser OFF: %0.2f and D1 Laser OFF %0.2f', D2GammaParameters(3), D1GammaParameters(3));

D2ModelFit = readmatrix('D2ModelFit.csv');
fprintf('\n Model fit Data from D2 Laser ON: %0.2f and D2 Laser OFF %0.2f', D2ModelFit(2), D2ModelFit(1));

D1ModelFit = readmatrix('D1ModelFit.csv');
fprintf('\n Model fit Data from D1 Laser ON: %0.2f and D1 Laser OFF %0.2f', D1ModelFit(2), D1ModelFit(1));
fprintf('\n\n')



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