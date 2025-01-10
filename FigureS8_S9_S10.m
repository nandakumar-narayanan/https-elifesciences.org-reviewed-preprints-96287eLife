
%-----------------------------------------------------------------------------------------------------------
% Written by Austin Bruce / Rewritten by Youngcho / checked by Kumar
%% Narayanan Lab

addpath(genpath('./analytical_tools_striatalcircuits/'))
%addpath(genpath('./analysis_tools'))
addpath('./behaviorCode/');

clear all;
cmap = [0 0 0; .9 .1 .1; .1 .1 .75; .5 .5 .5; .4 .4 .4];

medpc_files = 'pharm_data_WT/';
addpath(medpc_files)

centers = 0:0.1:20;
tinterval = centers(2:end);


% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loads OPTO Data
local_dataOpto = 'MPChalo';
addpath(local_dataOpto)

[MSN, group] = getAnProfile_RAB_opto_D1D2_v5();
mpcParsed = getDataIntr(local_dataOpto, MSN, group);
TrialAnSt = getTrialData_v3(mpcParsed);
GroupAnSt = plotGroup(TrialAnSt, group);

figure(56); set(gcf,'Color', 'White');  
subplot(2,2,1); 
hold on;     
data = squeeze(GroupAnSt.OptoD2CtrlHalo.cdf_switch(:,2,1,:))'; % time  x type x (depart,arrive) x animal,  cdf_switch(:,i_type,i_var,uniNameIdx)
plotband(tinterval,mean(data,1),std(data,[],1)/sqrt(size(data,1)),cmap(4,:))
data = squeeze(GroupAnSt.OptoD2CtrlHalo.cdf_switch(:,1,1,:))'; % time x (depart,arrive) x type x animal
plotband(tinterval,mean(data,1),std(data,[],1)/sqrt(size(data,1)),cmap(2,:))
ylim([0,1]); ylabel('Cumulative Density'); xlabel('Time from Start (in Seconds)'); set(gca, 'xtick', [0 6 18], 'ytick', [0 1])
hold off

subplot(2,2,2); 
data = GroupAnSt.OptoD2CtrlHalo.totSwitchTrDe;
data = cell2mat(cellfun(@(x) mean(x,1,'omitnan')',data,'UniformOutput',false))';
LadderPlot(num2cell(data(:,2)), num2cell(data(:,1)),1, 'D2CTLSwitchTime'); ylim([6.5 13]); 


subplot(2,2,3); 
hold on;     
data = squeeze(GroupAnSt.OptoD1CtrlHalo.cdf_switch(:,2,1,:))'; % time  x type x (depart,arrive) x animal,  cdf_switch(:,i_type,i_var,uniNameIdx)
plotband(tinterval,mean(data,1),std(data,[],1)/sqrt(size(data,1)),cmap(4,:))
data = squeeze(GroupAnSt.OptoD1CtrlHalo.cdf_switch(:,1,1,:))'; % time x (depart,arrive) x type x animal
plotband(tinterval,mean(data,1),std(data,[],1)/sqrt(size(data,1)),cmap(3,:))
ylim([0,1]); ylabel('Cumulative Density'); xlabel('Time from Start (in Seconds)'); set(gca, 'xtick', [0 6 18], 'ytick', [0 1])
hold off

subplot(2,2,4);
data = GroupAnSt.OptoD1CtrlHalo.totSwitchTrDe;
data = cell2mat(cellfun(@(x) mean(x,1,'omitnan')',data,'UniformOutput',false))';
LadderPlot(num2cell(data(:,2)), num2cell(data(:,1)),2, 'D1CTLSwitchTime'); ylim([6.5 13]); 
%%



% nose poke duration
% traversal time between the short and long response port
figure(57); set(gcf,'Color', 'White');  
subplot(2,3,2); 
data = GroupAnSt.OptoD2Halo.totRspDuration;
data = cell2mat(cellfun(@(x) mean(x,1,'omitnan')',data,'UniformOutput',false))';
LadderPlot(num2cell(data(:,2)), num2cell(data(:,1)),1, 'D2NosepokeDuration(seconds)'); % ylim([6.5 13]); 
subplot(2,3,3); 
data = GroupAnSt.OptoD2Halo.totTraversalTime;
data = cell2mat(cellfun(@(x) mean(x,1,'omitnan')',data,'UniformOutput',false))';
LadderPlot(num2cell(data(:,2)), num2cell(data(:,1)),1, 'D2TraversalTime(seconds)'); % ylim([6.5 13]); 
subplot(2,3,5); 
data = GroupAnSt.OptoD1Halo.totRspDuration;
data = cell2mat(cellfun(@(x) mean(x,1,'omitnan')',data,'UniformOutput',false))';
LadderPlot(num2cell(data(:,2)), num2cell(data(:,1)),2, 'D1NosepokeDuration(seconds)'); % ylim([6.5 13]); 
subplot(2,3,6); 
data = GroupAnSt.OptoD1Halo.totTraversalTime;
data = cell2mat(cellfun(@(x) mean(x,1,'omitnan')',data,'UniformOutput',false))';
LadderPlot(num2cell(data(:,2)), num2cell(data(:,1)),2, 'D1TraversalTime(seconds)'); % ylim([6.5 13]); 



subplot(2,3,1); 
hold on; 
npinterval = [0:0.01:0.8];
data = GroupAnSt.OptoD2Halo.totRspDuration;
data = cell2mat(cellfun(@(x) mean(x,1,'omitnan')',data,'UniformOutput',false))';
[x,fi] = ksdensity(data(:,1),npinterval);   
pH = patch([0 npinterval max(npinterval)], [0 x 0], 'k');
set(pH, 'FaceColor', [0.8 0 0], 'FaceAlpha', [0.5]); xlim([0 0.8]); 
[x,fi] = ksdensity(data(:,2),npinterval); 
pH = patch([0 npinterval max(npinterval)], [0 x 0], 'k');
set(pH, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', [0.5]); xlim([0 0.8]); 
hold off
ylabel('NP duration (Probability');  
xlabel('Time (Seconds'); 

subplot(2,3,4); 
hold on; 
npinterval = [0:0.01:0.8];
data = GroupAnSt.OptoD1Halo.totRspDuration;
data = cell2mat(cellfun(@(x) mean(x,1,'omitnan')',data,'UniformOutput',false))';
[x,fi] = ksdensity(data(:,1),npinterval);   
pH = patch([0 npinterval max(npinterval)], [0 x 0], 'k');
set(pH, 'FaceColor', [0 0 0.8], 'FaceAlpha', [0.5]); xlim([0 0.8]); 
[x,fi] = ksdensity(data(:,2),npinterval); 
pH = patch([0 npinterval max(npinterval)], [0 x 0], 'k');
set(pH, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', [0.5]); xlim([0 0.8]); 
hold off
ylabel('NP duration (Probability');  
xlabel('Time (Seconds'); 


figure(59) ; set(gcf,'Color', 'White');  

subplot(2,2,1); 
data = GroupAnSt.OptoD2Halo.totSwitchTrDe;
data = cell2mat(cellfun(@(x) std(x,1,'omitnan')',data,'UniformOutput',false))';
LadderPlot(num2cell(data(:,2)), num2cell(data(:,1)),1, sprintf('D2SwitchTimeSTD'));  ylim([2 5]); 
set(gca, 'ytick', [2 5])

subplot(2,2,2); 
data = GroupAnSt.OptoD2Halo.totReward;
data = cell2mat(cellfun(@(x) sum(~isnan(x))',data,'UniformOutput',false))';
LadderPlot(num2cell(data(:,2)), num2cell(data(:,1)),1, 'D2Reward#'); % ylim([6.5 13]); 
set(gca, 'ytick', [0 140])

subplot(2,2,3); 
data = GroupAnSt.OptoD1Halo.totSwitchTrDe;
data = cell2mat(cellfun(@(x) std(x,1,'omitnan')',data,'UniformOutput',false))';
LadderPlot(num2cell(data(:,2)), num2cell(data(:,1)),2, sprintf('D1SwitchTimeSTD'));  ylim([2 5]); 
set(gca, 'ytick', [2 5])

subplot(2,2,4); 
data = GroupAnSt.OptoD1Halo.totReward;
data = cell2mat(cellfun(@(x) sum(~isnan(x))',data,'UniformOutput',false))';
LadderPlot(num2cell(data(:,2)), num2cell(data(:,1)),2, 'D1Reward#'); % ylim([6.5 13]); 
set(gca, 'ytick', [0 140])





