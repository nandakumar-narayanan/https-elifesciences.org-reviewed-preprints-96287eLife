% Figure 1; Hannah Stutt, Matt Weber, and Kumar Narayanan; plots basic
% behavioral data from 30 animals

% Sets up directories
addpath(genpath('./analytical_tools_striatalcircuits/'))
addpath('./behaviorCode/');
medpc_files = './BehavioralData/';


% Loads behavioral data from 30 animals
[MSN, group] = getAnProfile_Learning_paper();
mpcParsed = getDataIntr(medpc_files, MSN, group);
TrialAnSt = getTrialData(mpcParsed);
GroupAnSt = plotGroup(TrialAnSt, group);
cmap = [0 0 0; .9 .1 .1; .1 .1 .75; .5 .5 .5; .4 .4 .4];

centers = 0:0.1:20;
tinterval = centers(2:end);

% Gets individual data for short and long NPs
clear shorts longs shortscdf; 
for i_animal=1:size(GroupAnSt.SITI_RI_HS.ids,2) 
[s1,~] = ksdensity(GroupAnSt.SITI_RI_HS.pooled_short_rsp{i_animal},centers(2:end),'Bandwidth',1.0);
[l1,~] = ksdensity(GroupAnSt.SITI_RI_HS.pooled_long_rsp{i_animal},centers(2:end),'Bandwidth',1.0);
[s1cdf,~] = histcounts(GroupAnSt.SITI_RI_HS.pooled_short_rsp{i_animal},centers,'Normalization','cdf');
shorts(i_animal,:) = s1; 
shortscdf(i_animal,:) = s1cdf; 
longs(i_animal,:) = l1; 
end

shorts = shorts'; longs = longs'; 

%Plots short, switch, and long responses
figure(1); set(gcf, 'Color', 'White'); hold on; 
data=squeeze(GroupAnSt.SITI_RI_HS.pdf_switch(:,1,1,:));
plotband(tinterval,mean(data,2),std(data,[],2)/sqrt(size(data,1)),cmap(1,:))
plotband(tinterval,mean(shorts,2),std(shorts,[],2)/sqrt(size(shorts,1)),cmap(2,:))
plotband(tinterval,mean(longs,2),std(longs,[],2)/sqrt(size(longs,1)),cmap(3,:))
ylim([0,0.15]); ylabel(' Density'); xlabel('Time from Start (in Seconds)'); set(gca, 'xtick', [0 6 18]);   set(gca, 'xtick', [0 6 18]); box off; 

% Plots animal-by-animal CDFs; 
figure(11); set(gcf, 'Color', 'White');
data=squeeze(GroupAnSt.SITI_RI_HS.cdf_switch(:,1,1,:));
plot(tinterval, data, 'k')
ylim([0,1]); ylabel('Cumulative Density'); xlabel('Time from Start (in Seconds)'); set(gca, 'xtick', [0 6 18]);   set(gca, 'xtick', [0 6 18]); box off; 

% Plots  CDFS
figure(12); set(gcf, 'Color', 'White');
data=squeeze(GroupAnSt.SITI_RI_HS.cdf_switch(:,1,1,:));
plotband(tinterval,mean(data,2),std(data,[],2)/sqrt(size(data,1)),cmap(1,:))
ylim([0,1]); ylabel('Cumulative Density'); xlabel('Time from Start (in Seconds)'); set(gca, 'xtick', [0 6 18]);   set(gca, 'xtick', [0 6 18]); box off; 

data = GroupAnSt.SITI_RI_HS.totSwitchTrDe;  % grab the number of switch departures
data = cell2mat(cellfun(@(x) mean(x,1,'omitnan')',data,'UniformOutput',false))'; data = mean(data'); % mean across trials; 
fprintf('\n# Mean Switch Time %0.1f (%0.1f-%0.1f)', median(data), quantile(data,0.25),quantile(data,0.75))


for i= 1:length((GroupAnSt.SITI_RI_HS.totSwitchTrDe))
    trialNumber(i) = sum(isfinite([GroupAnSt.SITI_RI_HS.totSwitchTrDe{i}(:,1); GroupAnSt.SITI_RI_HS.totSwitchTrDe{i}(:,2)]));
end
fprintf('\n# Trials %0.0f (%0.0f-%0.0f)', median(trialNumber), quantile(trialNumber,0.25),quantile(trialNumber,0.75))

fprintf('\n D2/D1 cre mouse behavioral stats in Figure 3')
fprintf('\n Trials are scattered throughout the code for Table S1')



