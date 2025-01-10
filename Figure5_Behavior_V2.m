
%-----------------------------------------------------------------------------------------------------------
% Written by Austin Bruce / Rewritten by Youngcho / checked by Kumar
%% Narayanan Lab

addpath(genpath('./analytical_tools_striatalcircuits/'))
%addpath(genpath('./analysis_tools'))
addpath('./behaviorCode/');

clear all;
cmap = [0 0 0; .9 .1 .1; .1 .1 .75; .5 .5 .5; .4 .4 .4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PHARMACOLOGY DATA
medpc_files = 'pharm_data_WT/';
addpath(medpc_files)

[MSN, group] = getAnProfile_D1D2_antagonist_V2();
mpcParsed = getDataIntr(medpc_files, MSN, group);
TrialAnSt = getTrialData_v3(mpcParsed); 
GroupAnStPharm = plotGroup(TrialAnSt, group);

centers = 0:0.1:20;
tinterval = centers(2:end);

% Reviewer wanted us to compare saline and laser off; 
dataD2Sal = GroupAnStPharm.Saline_SUL.totSwitchTrDe;
dataD2Sal = cell2mat(cellfun(@(x) mean(x,1,'omitnan')',dataD2Sal,'UniformOutput',false))';

dataD1Sal = GroupAnStPharm.Saline_SCH.totSwitchTrDe;
dataD1Sal = cell2mat(cellfun(@(x) mean(x,1,'omitnan')',dataD1Sal,'UniformOutput',false))';


% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  OPTO Data
local_dataOpto = 'MPChalo';
addpath(local_dataOpto)

[MSN, group] = getAnProfile_RAB_opto_D1D2_v5();
mpcParsed = getDataIntr(local_dataOpto, MSN, group);
TrialAnSt = getTrialData_v3(mpcParsed);
GroupAnSt = plotGroup(TrialAnSt, group);


figure(5);  
subplot(2,4,1); 
hold on;     
% Note that on Line 87 of GetTrialData_V3, opto==1 is specified as laser
% OFF, which is why the 2nd column is laser off, and the 1st column is
% laser on
data = squeeze(GroupAnSt.OptoD2Halo.pdf_switch(:,2,1,:))'; % time  x type x (depart,arrive) x animal,  cdf_switch(:,i_type,i_var,uniNameIdx)
plotband(tinterval,mean(data,1),std(data,[],1)/sqrt(size(data,1)),cmap(4,:))
data = squeeze(GroupAnSt.OptoD2Halo.pdf_switch(:,1,1,:))'; % time x (depart,arrive) x type x animal
plotband(tinterval,mean(data,1),std(data,[],1)/sqrt(size(data,1)),cmap(2,:))
ylim([0,0.18]); ylabel('Probability Density'); xlabel('Time from Start (in Seconds)'); set(gca, 'xtick', [0 6 18]);
hold off


fprintf('\n\n**Opto***')


hold off

subplot(2,4,2); 
data = GroupAnSt.OptoD2Halo.totSwitchTrDe;
data = cell2mat(cellfun(@(x) mean(x,1,'omitnan')',data,'UniformOutput',false))';
LadderPlot(num2cell(data(:,2)), num2cell(data(:,1)),1, 'D2OptoSwitchTime(seconds)'); ylim([6.5 13]); 

subplot(2,4,3); 
hold on;     
data = squeeze(GroupAnSt.OptoD1Halo.pdf_switch(:,2,1,:))'; % time  x type x (depart,arrive) x animal,  cdf_switch(:,i_type,i_var,uniNameIdx)
plotband(tinterval,mean(data,1),std(data,[],1)/sqrt(size(data,1)),cmap(4,:))
data = squeeze(GroupAnSt.OptoD1Halo.pdf_switch(:,1,1,:))'; % time x (depart,arrive) x type x animal
plotband(tinterval,mean(data,1),std(data,[],1)/sqrt(size(data,1)),cmap(3,:))
ylim([0,0.18]); ylabel('Probability Density'); xlabel('Time from Start (in Seconds)'); set(gca, 'xtick', [0 6 18]);
hold off


subplot(2,4,4); 
data = GroupAnSt.OptoD1Halo.totSwitchTrDe;
data = cell2mat(cellfun(@(x) mean(x,1,'omitnan')',data,'UniformOutput',false))';
LadderPlot(num2cell(data(:,2)), num2cell(data(:,1)),2, 'D1OptoSwitchTime(seconds)'); ylim([6.5 13]); 

fprintf('\n\n**Pharm')
subplot(2,4,5); hold on; 
data = squeeze(GroupAnStPharm.Saline_SUL.pdf_switch(:,1,1,:))'; 
plotband(tinterval,mean(data,1),std(data,[],1)/sqrt(size(data,1)),cmap(4,:))
data = squeeze(GroupAnStPharm.Sulpiride.pdf_switch(:,1,1,:))'; 
plotband(tinterval,mean(data,1),std(data,[],1)/sqrt(size(data,1)),cmap(2,:))
ylim([0,0.18]); ylabel('Probability Density'); xlabel('Time from Start (in Seconds)'); set(gca, 'xtick', [0 6 18]);   set(gca, 'xtick', [0 6 18]); 

subplot(2,4,6); 
data = GroupAnStPharm.Saline_SUL.totSwitchTrDe;
data1 = cell2mat(cellfun(@(x) mean(x,1,'omitnan')',data,'UniformOutput',false))';
trials1 = cell2mat(cellfun(@(x) length(x)',data,'UniformOutput',false))';

data = GroupAnStPharm.Sulpiride.totSwitchTrDe;
data2 = cell2mat(cellfun(@(x) mean(x,1,'omitnan')',data,'UniformOutput',false))';
trials2 = cell2mat(cellfun(@(x) length(x)',data,'UniformOutput',false))';

LadderPlot(num2cell(data1), num2cell(data2),1, 'D2PharmSwitchTime(seconds)'); ylim([6.5 13]); 


fprintf('\n Saline # Trials %2.0f (%2.0f-%2.0f)',  median(trials1), quantile(trials1,0.25),quantile(trials1,0.75))
fprintf('\n SUL # Trials %2.0f (%2.0f-%2.0f)',  median(trials2), quantile(trials2,0.25),quantile(trials2,0.75))

data1 = GroupAnSt.OptoD2Halo.totSwitchTrDe;
data1 = cell2mat(cellfun(@(x) mean(x,1,'omitnan')',data1,'UniformOutput',false))';
data1 = data1(:,2);
fprintf('\n D2 Laser Off vs Saline %0.2f ', ranksum(data1,dataD2Sal));


subplot(2,4,7); hold on; 
data = squeeze(GroupAnStPharm.Saline_SCH.pdf_switch(:,1,1,:))'; 
plotband(tinterval,mean(data,1),std(data,[],1)/sqrt(size(data,1)),cmap(1,:))
data = squeeze(GroupAnStPharm.SCH23390.pdf_switch(:,1,1,:))'; 
plotband(tinterval,mean(data,1),std(data,[],1)/sqrt(size(data,1)),cmap(3,:))
ylim([0,0.18]); ylabel('Probability Density'); xlabel('Time from Start (in Seconds)'); set(gca, 'xtick', [0 6 18]);

subplot(2,4,8); 
data = GroupAnStPharm.Saline_SCH.totSwitchTrDe;
data1 = cell2mat(cellfun(@(x) mean(x,1,'omitnan')',data,'UniformOutput',false))';
trials1 = cell2mat(cellfun(@(x) length(x)',data,'UniformOutput',false))';

data = GroupAnStPharm.SCH23390.totSwitchTrDe;
data2 = cell2mat(cellfun(@(x) mean(x,1,'omitnan')',data,'UniformOutput',false))';
trials2 = cell2mat(cellfun(@(x) length(x)',data,'UniformOutput',false))';

LadderPlot(num2cell(data1), num2cell(data2),2, 'D1PharmSwitchTime(seconds)'); ylim([6.5 13]); ylim([6.5 13]); 
fprintf('\n Saline # Trials %2.0f (%2.0f-%2.0f)',  median(trials1), quantile(trials1,0.25),quantile(trials1,0.75))
fprintf('\n SCH # Trials %2.0f (%2.0f-%2.0f)',  median(trials2), quantile(trials2,0.25),quantile(trials2,0.75))

data1 = GroupAnSt.OptoD1Halo.totSwitchTrDe;
data1 = cell2mat(cellfun(@(x) mean(x,1,'omitnan')',data1,'UniformOutput',false))';
data1 = data1(:,2);
fprintf('\n D1 Laser Off vs Saline %0.2f ', ranksum(data1,dataD1Sal));


%%
behavioralTable = [];
groupFields = fieldnames(GroupAnSt);
for i_group = 1:2
    data = GroupAnSt.(groupFields{i_group}).totSwitchTrDe;
    ids = GroupAnSt.(groupFields{i_group}).ids;
    for i_ani = 1:size(data,2)
        RSPTIME = reshape(data{i_ani},[],1);
        LASER = reshape(repmat({'inhibition' 'ctrl'}, length(data{i_ani}), 1),[],1);
        ANIMALID = reshape(repmat(ids(i_ani), length(data{i_ani}),2),[],1);
        GENOTYPE = reshape(repmat(groupFields(i_group), length(data{i_ani}),2),[],1);
        m_table = table(RSPTIME,LASER,ANIMALID, GENOTYPE);    
        behavioralTable = [behavioralTable; m_table];
    end
end
behavioralTable(isnan(behavioralTable.RSPTIME),:) = [];
writetable(behavioralTable, 'OptoBehavioralTable_V3.csv');


% Reviewer wanted trial for Table 1
for i= 1:length((GroupAnSt.OptoD2Halo.totSwitchTrDe))
    trialNumber(i) = [sum(isfinite(GroupAnSt.OptoD2Halo.totSwitchTrDe{i}(:,1)))] + [sum(isfinite(GroupAnSt.OptoD2Halo.totSwitchTrDe{i}(:,2)))];
end
fprintf('\n\n**')
fprintf('\n D2Opto Trials %2.0f (%2.0f-%2.0f); ON %2.0f (%2.0f-%2.0f) ', median(trialNumber), quantile(trialNumber,0.25),quantile(trialNumber,0.75));

for i= 1:length((GroupAnSt.OptoD2CtrlHalo.totSwitchTrDe))
    trialNumber(i) = [sum(isfinite(GroupAnSt.OptoD2CtrlHalo.totSwitchTrDe{i}(:,1)))] + [sum(isfinite(GroupAnSt.OptoD2CtrlHalo.totSwitchTrDe{i}(:,2)))];
end
fprintf('\n D2Ctrl Trials %2.0f (%2.0f-%2.0f); ON %2.0f (%2.0f-%2.0f) ', median(trialNumber), quantile(trialNumber,0.25),quantile(trialNumber,0.75));


for i= 1:length((GroupAnSt.OptoD1Halo.totSwitchTrDe))
    trialNumber(i) = [sum(isfinite(GroupAnSt.OptoD1Halo.totSwitchTrDe{i}(:,1)))] + [sum(isfinite(GroupAnSt.OptoD1Halo.totSwitchTrDe{i}(:,2)))];
end
fprintf('\n D1Opto Trials %2.0f (%2.0f-%2.0f); ON %2.0f (%2.0f-%2.0f) ', median(trialNumber), quantile(trialNumber,0.25),quantile(trialNumber,0.75));

for i= 1:length((GroupAnSt.OptoD1CtrlHalo.totSwitchTrDe))
    trialNumber(i) = [sum(isfinite(GroupAnSt.OptoD1CtrlHalo.totSwitchTrDe{i}(:,1)))] + [sum(isfinite(GroupAnSt.OptoD1CtrlHalo.totSwitchTrDe{i}(:,2)))];
end
fprintf('\n D1Ctrl Trials %2.0f (%2.0f-%2.0f); ON %2.0f (%2.0f-%2.0f) ', median(trialNumber), quantile(trialNumber,0.25),quantile(trialNumber,0.75));


