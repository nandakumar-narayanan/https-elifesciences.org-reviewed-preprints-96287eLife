

rootdir='./'; 

fprintf('\n For Table S2:')

D2ONGammaParameters = readmatrix('D2ONGammaParameters.csv');
fprintf('\n>> D2 MSN inhibited model: Gamma distribution parameters: Alpha %2.2f; Beta %2.2f,', D2ONGammaParameters(1), D2ONGammaParameters(2));

D1ONGammaParameters = readmatrix('D1ONGammaParameters.csv');
fprintf('\n>> D1 MSN inhibited model: Gamma distribution parameters: Alpha %2.2f; Beta %2.2f,', D1ONGammaParameters(1), D1ONGammaParameters(2));

D2OFFGammaParameters = readmatrix('D2OFFGammaParameters.csv');
fprintf('\n>> D2 MSN Laser Off model: Gamma distribution parameters: Alpha %2.2f; Beta %2.2f,', D2OFFGammaParameters(1), D2OFFGammaParameters(2));

D1OFFGammaParameters = readmatrix('D1OFFGammaParameters.csv');
fprintf('\n>> D1 MSN Laser Off model: Gamma distribution parameters: Alpha %2.2f; Beta %2.2f,', D1OFFGammaParameters(1), D1OFFGammaParameters(2));


for i = 1:2
    if i==1
        opto='D1'; 
    else 
        opto='D2';
    end
clearvars -except T i opto rootdir
%%  Load the experimental data and  the modeling data
T = readtable('OptoBehavioralTableV3.csv');

if strcmp(opto,'D1')
    % keep only D1-mice
    DcreMice = (strcmp(T.GENOTYPE,'OptoD1Halo'));
    colorOFF = [0.3 0.75 0.93]; % light blue
    colorON = 'b' ;
    %
    FR_target = 0; baseline = 0.48;
    D_off = 0.1405; sigmaNoise_off= 0.0525;
    D_on = 0.1224;  sigmaNoise_on = 0.0435;
    %
elseif strcmp(opto,'D2')    
    % keep only D2-mice
    DcreMice = (strcmp(T.GENOTYPE,'OptoD2Halo'));
    colorOFF = [1 0.7 0.7]; % light red
    colorON = 'r' ;
    %
    FR_target = 1; baseline = 0.52;
    D_off = 0.135; sigmaNoise_off= 0.0525;
    D_on = 0.1287;  sigmaNoise_on = 0.0435;
    %
end
T = T(DcreMice,:);

%% Model
FileName= ['FR' num2str(FR_target) '-b0' num2str(baseline*100) '.mat'];  
load([rootdir FileName],...
            'estimated_muHat','estimated_stdHat', 'D', 'sigmaNoise'); 
estimated_muHat = squeeze(median(estimated_muHat,3,"omitnan")); 
           % matrix of same dimension as D x sigmaNoise
estimated_cvHat = ...
 squeeze(median(estimated_stdHat,3,"omitnan")) ./ estimated_muHat ;
fprintf('\n\n>> Model parameters: FR Target  %0.1f baseline %0.2f', FR_target,baseline)



%%  Calculate mean and coeff of variation for switch times

[dataLaserOff,m_LaserOff,md_LaserOff,cv_LaserOff]= AnimalData(T,'ctrl');
[dataLaserOn,m_LaserOn,md_LaserOn, cv_LaserOn] = AnimalData(T,'inhibition');

%---------- CRTL -------------------------
pdL_Off = fitdist(dataLaserOff,'Gamma');
alpha_off = pdL_Off.a;   % shape 
beta_off = 1/pdL_Off.b;  % form
% -- fit-mean and fit-CV calculated from the Gamma distrib aprox
mu_off = alpha_off / beta_off ;
cv_off = 1/sqrt(alpha_off);

fprintf(['\n\n>> CTRL - ' opto ' (empirical data)']);
fprintf('\n>> Gamma distribution parameters: Alpha %2.2f; Beta %2.2f,', alpha_off, beta_off);
fprintf('\n>> Gamma distribution stats: Mean %2.2f; CV %2.2f,', mu_off, cv_off);
fprintf('\n>> Empirical values from data: Mean %2.2f; CV %2.2f,', m_LaserOff, cv_LaserOff);

%---------- INHIBITION -------------------------
pdL_On = fitdist(dataLaserOn,'Gamma');
alpha_on = pdL_On.a;   % shape 
beta_on = 1/pdL_On.b;  % form
% -- fit-mean and fit-CV calculated from the Gamma distrib aprox
mu_on = alpha_on / beta_on ;
cv_on = 1/sqrt(alpha_on);

fprintf(['\n\n>> Inhibition - ' opto ' (empirical data)']);
fprintf('\n>> Gamma distribution parameters: Alpha %2.2f; Beta %2.2f,', alpha_on, beta_on);
fprintf('\n>> Gamma distribution stats: Mean %2.2f; CV %2.2f,', mu_on, cv_on);
fprintf('\n>> Empirical values from data: Mean %2.2f; CV %2.2f,', m_LaserOn, cv_LaserOn) 
fprintf('\n');

%% Draw preliminary diagram to see which values of Ta, D may work
% create matrix coordinates
Dcoord = repmat(D',1,length(sigmaNoise));
SNcoord = repmat(sigmaNoise,length(D),1);

% tranform into column vector
Dcoord = Dcoord(:);
SNcoord = SNcoord(:); 


%% Compare with experimental data from both LaserOff and LaserOn

%figure('units','normalized','outerposition',[0 0 0.25 0.5])
figure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(2,2,1)
    % CTRL
    % ------ plot the values that match the mean/ for LaserOff
    clims = [-0.2 0.2];
    imagesc(D,sigmaNoise, (estimated_muHat'- mu_off)/mu_off, clims);
    set(gca,'YDir','normal');
    nJet = 1000; % needs to be even, high number for better interpolation
    mapJet = jet(nJet);
    newmap= jet(nJet);
    map=[[0.5 0.5 0.5];newmap(1:100,:); mapJet(100:end,:)];
    colormap(map);
    colorbar('northoutside')
    hold on
    % show the parameter values used in Figure 4/eLife paper      
    scatter(D_off,sigmaNoise_off,50,...
        'MarkerEdgeColor','k','MarkerFaceColor','k','Marker','square'); 
    % add details to figure
    ylabel('$\sigma$','Interpreter','Latex','FontSize', 12)
    xlabel('$D$','Interpreter','Latex','FontSize', 12)
    title([opto ' - ctrl - Mean'],'FontSize',14)
    
    set(findall(gcf,'-property','FontSize','-property','LineWidth'),...
            'FontSize',12,'LineWidth',2, 'FontName','Arial'); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CTRL
    % ------ plot the values that match the mean/ for LaserOff
    subplot(2,2,2)
    clims = [-0.14 0.14];
    imagesc(D,sigmaNoise, estimated_cvHat'- cv_off , clims);
    set(gca,'YDir','normal');
    nJet = 1000; % needs to be even, high number for better interpolation
    mapJet = jet(nJet);
    newmap= jet(nJet);
    map=[[0.5 0.5 0.5];newmap(1:100,:); mapJet(100:end,:)];
    colormap(map);
    colorbar('northoutside')
    hold on
    scatter(D_off,sigmaNoise_off,50,...
        'MarkerEdgeColor','k','MarkerFaceColor','k','Marker','square'); 
    ylabel('$\sigma$','Interpreter','Latex','FontSize', 12)
    xlabel('$D$','Interpreter','Latex','FontSize', 12)
    title([opto ' - ctrl - CV'],'FontSize',14)
    
    set(findall(gcf,'-property','FontSize','-property','LineWidth'),...
            'FontSize',12,'LineWidth',2, 'FontName','Arial'); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(2,2,3)
    % INHIBITION
    % ------ plot the values that match the mean/ for LaserOn
    clims = [-0.2 0.2];
    imagesc(D,sigmaNoise, (estimated_muHat'- mu_on)/mu_on, clims);
    set(gca,'YDir','normal');
    nJet = 1000; % needs to be even, high number for better interpolation
    mapJet = jet(nJet);
    newmap= jet(nJet);
    map=[[0.5 0.5 0.5];newmap(1:100,:); mapJet(100:end,:)];
    colormap(map);
    colorbar('northoutside')
    hold on
    % show the parameter values used in Figure 4/eLife paper      
    scatter(D_on,sigmaNoise_on,50,...
        'MarkerEdgeColor','k','MarkerFaceColor','k','Marker','^'); 
   % add details to figure
    ylabel('$\sigma$','Interpreter','Latex','FontSize', 12)
    xlabel('$D$','Interpreter','Latex','FontSize', 12)
    title([opto ' - inhibition - Mean'],'FontSize',14)
  
   set(findall(gcf,'-property','FontSize','-property','LineWidth'),...
            'FontSize',12,'LineWidth',2, 'FontName','Arial'); 

    
   subplot(2,2,4)
    % INHIBITION
    % ------ plot the values that match the mean/ for LaserOn
    clims = [-0.14 0.14];
    imagesc(D,sigmaNoise, estimated_cvHat'- cv_on , clims);
    set(gca,'YDir','normal');
    nJet = 1000; % needs to be even, high number for better interpolation
    mapJet = jet(nJet);
    newmap= jet(nJet);
    map=[[0.5 0.5 0.5];newmap(1:100,:); mapJet(100:end,:)];
    colormap(map);
    colorbar('northoutside')
    hold on
    % show the parameter values used in Figure 4/eLife paper      
    scatter(D_on,sigmaNoise_on,50,...
        'MarkerEdgeColor','k','MarkerFaceColor','k','Marker','^'); 

   % add details to figure
    ylabel('$\sigma$','Interpreter','Latex','FontSize', 12)
    xlabel('$D$','Interpreter','Latex','FontSize', 12)
    title([opto ' - inhibition - CV'],'FontSize',14)
  
   set(findall(gcf,'-property','FontSize','-property','LineWidth'),...
            'FontSize',12,'LineWidth',2, 'FontName','Arial'); 
end

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
