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
plot(xx,cdf(pdL2_Off,xx),'--','color', 'k','Linewidth',1) % for gamma model
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
plot(xx,cdf(pdL1_On,xx),':','color',[0 0 0],'Linewidth',2) % for gamma data
%
plot(xx,cdf(pdL2_On,xx),'--','color', [1 0 0],'Linewidth',1) % for gamma model
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
plot(xx,cdf(pdL1_On,xx),':','color',[0 0 0],'Linewidth',2) % for gamma data
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