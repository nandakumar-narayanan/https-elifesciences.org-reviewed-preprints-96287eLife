%% Need to run Figure 5 first; 
load neurons_pharm.mat; % Loads the Full Dataset
avgFR = [neurons.avgFR];
num_spike = [neurons.numTS]; 
% valid_idx = avgFR>0.5 & avgFR<20 & strcmp({neurons.type},'MSN') & num_spike > 5000;
valid_idx = avgFR>0.5 & avgFR<20 & strcmp({neurons.type},'MSN');
fprintf('\n\nTotal neurons for analysis in the full dataset %d\n',sum(valid_idx));
binSize = 0.2;  
bw = 1; 

salineNeurons = neurons([neurons.condition]==1 & valid_idx); 
D2Block = neurons([neurons.condition]==2 & valid_idx);
D1Block = neurons([neurons.condition]==3 & valid_idx);

salineStats=EnsembleAnalysisSwitch(salineNeurons, bw, binSize);  % Currently, 20 trial minimum; does better with more trials, obviously
D2Stats=EnsembleAnalysisSwitch(D2Block, bw, binSize);  % Currently, 20 trial minimum; does better with more trials, obviously
D1Stats=EnsembleAnalysisSwitch(D1Block, bw, binSize);  % Currently, 20 trial minimum; does better with more trials, obviously


figure(3); 
times = D1Stats.times(D1Stats.goodtimes==1); b = [times times times];

subplot(2,4,5); cla; 
imagesc(times, times, salineStats.y(salineStats.goodtimes==1,salineStats.goodtimes==1)); axis xy; 
set(gca, 'xtick', [0 6 18], 'xticklabel', [0 6 18]);  set(gca, 'ytick', [0 6 18], 'yticklabel', [0 6 18]); 
xlabel('Objective Time'); ylabel('Predicted Time');  colorbar; 

subplot(2,4,6); cla; 
imagesc(times, times, D2Stats.y(D2Stats.goodtimes==1,D2Stats.goodtimes==1)); axis xy;  colorbar; 
set(gca, 'xtick', [0 6 18], 'xticklabel', [0 6 18]);  set(gca, 'ytick', [0 6 18], 'yticklabel', [0 6 18]);  
xlabel('Objective Time'); ylabel('Predicted Time');  colorbar; 

subplot(2,4,7); cla;   colorbar; 
imagesc(times, times, D1Stats.y(D1Stats.goodtimes==1,D1Stats.goodtimes==1)); axis xy;
set(gca, 'xtick', [0 6 18], 'xticklabel', [0 6 18]);  set(gca, 'ytick', [0 6 18], 'yticklabel', [0 6 18]); 
xlabel('Objective Time'); ylabel('Predicted Time');  colorbar; 

% set(gca, 'xtick', [1 7 19 20 26 38 39 45 57], 'xticklabel', b([1 7 19 20 26 38 39 45 57])); 
% set(gca, 'ytick', [1 7 19], 'yticklabel', times([1 7 19])); ylabel('Predicted Time');
% set(gca, 'xtick', [7 26 45], 'xticklabel', {'Saline'; 'D2 Blockade'; 'D1 Blockade'});

subplot(2,4,8); cla; clear data; 
data{1} =  salineStats.r2_0_6; 
data{2} =  D2Stats.r2_0_6; 
data{3} =  D1Stats.r2_0_6; 

data{4} =  salineStats.r2_6_12; 
data{5} =  D2Stats.r2_6_12; 
data{6} =  D1Stats.r2_6_12; 

data{7} =  salineStats.r2_12_18; 
data{8} =  D2Stats.r2_12_18; 
data{9} =  D1Stats.r2_12_18; 


jitterPlot(data, [0.7 0.7 0.7; 1 0.7 0.7; 0.7 0.7 1]); 
set(gca, 'xtick', [2:3:8], 'xticklabel', {'0-6 s'; '6-12 s'; '6-12s';}); 

ylabel('R2')
fprintf('\n Early: Saline %0.2f (%0.2f-%0.2f)', median(salineStats.r2_0_6), quantile(salineStats.r2_0_6,0.25), quantile(salineStats.r2_0_6, 0.75) )
fprintf('\n Mid: Saline %0.2f (%0.2f-%0.2f)', median(salineStats.r2_6_12), quantile(salineStats.r2_6_12,0.25), quantile(salineStats.r2_6_12, 0.75) )
fprintf('\n Late: Saline %0.2f (%0.2f-%0.2f)', median(salineStats.r2_12_18), quantile(salineStats.r2_12_18,0.25), quantile(salineStats.r2_12_18, 0.75) )
fprintf('\n Decoding Mid vs Early: ranksum 0-6 V. 6-12   p = %0.6f, cohen d %2.1f', ranksum(salineStats.r2_6_12, salineStats.r2_0_6), cohend(salineStats.r2_6_12, salineStats.r2_0_6))
fprintf('\n Decoding Late vs Early: ranksum 0-6 V. 12-18 p = %0.6f, cohen d %2.1f', ranksum(salineStats.r2_12_18, salineStats.r2_0_6), cohend(salineStats.r2_12_18, salineStats.r2_0_6))
csvwrite('./stats/EarlyDecodingSaline.csv', salineStats.r2_0_6');
csvwrite('./stats/MidDecodingSaline.csv', salineStats.r2_6_12');
csvwrite('./stats/LateDecodingSaline.csv', salineStats.r2_12_18');


fprintf('\n D2Stats %0.2f (%0.2f-%0.2f)', median(D2Stats.r2_0_6), quantile(D2Stats.r2_0_6,0.25), quantile(D2Stats.r2_0_6, 0.75) )
fprintf('\n D1Stats %0.2f (%0.2f-%0.2f)', median(D1Stats.r2_0_6), quantile(D1Stats.r2_0_6,0.25), quantile(D1Stats.r2_0_6, 0.75) )
fprintf('\n Early: Decoding: ranksum Saline V. D2 Block p = %0.5f, cohen d %2.1f', ranksum(salineStats.r2_0_6, D2Stats.r2_0_6), cohend(salineStats.r2_0_6, D2Stats.r2_0_6))
fprintf('\n Decoding: ranksum Saline V. D1 Block p = %0.3f, cohen d %2.1f', ranksum(salineStats.r2_0_6, D1Stats.r2_0_6), cohend(salineStats.r2_0_6, D1Stats.r2_0_6))
fprintf('\n Decoding: ranksum D2 V. D1 Block p = %0.3f, cohen d %2.1f', ranksum(D2Stats.r2_0_6, D1Stats.r2_0_6), cohend(D2Stats.r2_0_6, D1Stats.r2_0_6))
csvwrite('./stats/EarlyDecodingD2.csv', D2Stats.r2_0_6');
csvwrite('./stats/EarlyDecodingD1.csv', D1Stats.r2_0_6');

fprintf('\n Mid: Decoding: ranksum Saline V. D2 Block p = %0.3f, cohen d %2.1f', ranksum(salineStats.r2_6_12, D2Stats.r2_6_12), cohend(salineStats.r2_6_12, D2Stats.r2_6_12))
fprintf('\n Mid: Decoding: ranksum Saline V. D1 Block p = %0.3f, cohen d %2.1f', ranksum(salineStats.r2_6_12, D1Stats.r2_6_12), cohend(salineStats.r2_6_12, D1Stats.r2_6_12))
csvwrite('./stats/MidDecodingD2.csv', D2Stats.r2_6_12');
csvwrite('./stats/MidDecodingD1.csv', D1Stats.r2_6_12');


fprintf('\n Late: Decoding: ranksum Saline V. D2 Block p = %0.3f, cohen d %2.1f', ranksum(salineStats.r2_12_18, D2Stats.r2_12_18), cohend(salineStats.r2_12_18, D2Stats.r2_12_18))
fprintf('\n Late: Decoding: ranksum Saline V. D1 Block p = %0.3f, cohen d %2.1f', ranksum(salineStats.r2_12_18, D1Stats.r2_12_18), cohend(salineStats.r2_12_18, D1Stats.r2_12_18))
csvwrite('./stats/LateDecodingD2.csv', D2Stats.r2_6_12');
csvwrite('./stats/LateDecodingD1.csv', D1Stats.r2_6_12');


se = [salineStats.r2_0_6' 0*ones(1,length(salineStats.r2_0_6))' ones(1,length(salineStats.r2_0_6))' [1:length(salineStats.r2_0_6)]'];
sm = [salineStats.r2_6_12' 0*ones(1,length(salineStats.r2_6_12))' 2*ones(1,length(salineStats.r2_6_12))' [1:length(salineStats.r2_6_12)]'];
sl = [salineStats.r2_12_18' 0*ones(1,length(salineStats.r2_12_18))' 3*ones(1,length(salineStats.r2_12_18))' [1:length(salineStats.r2_12_18)]'];

d1e = [D1Stats.r2_0_6' 1*ones(1,length(D1Stats.r2_0_6))' ones(1,length(D1Stats.r2_0_6))' [1:length(D1Stats.r2_0_6)]'];
d1m = [D1Stats.r2_6_12' 1*ones(1,length(D1Stats.r2_6_12))' 2*ones(1,length(D1Stats.r2_6_12))' [1:length(D1Stats.r2_6_12)]'];
d1l = [D1Stats.r2_12_18' 1*ones(1,length(D1Stats.r2_12_18))' 3*ones(1,length(D1Stats.r2_12_18))' [1:length(D1Stats.r2_12_18)]'];


d2e = [D2Stats.r2_0_6' 2*ones(1,length(D2Stats.r2_0_6))' ones(1,length(D2Stats.r2_0_6))' [1:length(D2Stats.r2_0_6)]'];
d2m = [D2Stats.r2_6_12' 2*ones(1,length(D2Stats.r2_6_12))' 2*ones(1,length(D2Stats.r2_6_12))' [1:length(D2Stats.r2_6_12)]'];
d2l = [D2Stats.r2_12_18' 2*ones(1,length(D2Stats.r2_12_18))' 3*ones(1,length(D2Stats.r2_12_18))' [1:length(D2Stats.r2_12_18)]'];

t = [se; sm; sl; d1e; d1m; d1l; d2e; d2m; d2l];

t = table(t(:,1), t(:, 2), t(:,3), t(:,4), 'VariableNames', {'R2', 'Drug', 'Epoch', 'Trial'}); 
writetable(t, 'DecodingData.csv');
% 
% % Takes 10 minutes
% neuronsEnsembleAnalysesPharm
