%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supplement for all neurons, unmatched dataset
load neurons_pharm.mat; % Loads the Full Dataset
avgFR = [neurons.avgFR];
num_spike = [neurons.numTS]; 
% valid_idx = avgFR>0.5 & avgFR<20 & strcmp({neurons.type},'MSN') & num_spike > 5000;
valid_idx = avgFR>0.5 & avgFR<20 & strcmp({neurons.type},'MSN');
fprintf('\n\nTotal neurons for analysis in the full dataset %d\n',sum(valid_idx));
neurons=neurons(valid_idx);

% Interval parameters
intStart =  -4;
trialStart = 0; 
trialEnd = 6; 
intEnd = 22;
binSize = 0.2;  
bw = 1;  
interval = intStart:binSize:intEnd;
interval_bins = trialStart:binSize:trialEnd;


counter=1; clear *PETH* sPETH condition;
for i_neuron = 1:length(neurons) % 1
    periEventSpike = peSpike(neurons(i_neuron).spikeTS, neurons(i_neuron).events.switchTrialInit, interval);
    PETH(counter, :) = gksmooth(periEventSpike, interval, bw);  % Smooth for PCA
    condition(counter) = neurons(i_neuron).condition;
    counter=counter+1;
end

condition = [neurons.condition];

TimeInterval = find(interval>[trialStart-binSize]&interval<[trialEnd+binSize]);
cPETH = PETH(:, TimeInterval);  
cInterval=interval(TimeInterval); 
zPETH = zscore(cPETH')';

warning off;
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(zPETH,'NumComponents',3);
warning on; 

figure(62); set(gcf, 'Color', 'White')
subplot(1,3,1); cla; hold on; 
plot(interval(TimeInterval), COEFF(:,1)+0.3, 'k', 'LineWidth', 7); 
set(gca, 'ytick', []); xlabel('Time from  Start (seconds)')
axis tight;  set(gca, 'YDir', 'reverse'); % to match Fig 6

subplot(1,3,2); cla; hold on; 
plot(EXPLAINED./sum(EXPLAINED), 'ko-', 'MarkerFaceColor','k','MarkerSize', 10)
xlim([0 10]); ylabel('Fraction Variance Explained'); xlabel('Component #')

subplot(1,3,3); cla; hold on; 
vdata{1} = (SCORE(condition==1,1)); 
vdata{2} = (SCORE(condition==2,1)); 
vdata{3} = (SCORE(condition==3,1)); 
data = [vdata{1}; vdata{2}; vdata{3}]';
groups = [ones(1,length(vdata{1})) 2*ones(1,length(vdata{2})) 3*ones(1,length(vdata{3})) ];
set(gca, 'YDir', 'reverse'); % to match Fig 6

scatter(ones(size(vdata{1})).*(1+(rand(size(vdata{1}))-0.5)/10),vdata{1},'k','filled')
scatter(2*ones(size(vdata{2})).*(1+(rand(size(vdata{2}))-0.5)/10),vdata{2},'r','filled')
scatter(3*ones(size(vdata{3})).*(1+(rand(size(vdata{3}))-0.5)/10),vdata{3},'b','filled')
boxplot(data, groups); box off; 

% have to rebuild an ensmemble for each animal 
animalCounter = 1;  clear animalArray;
counter=1; clear *PETH* sPETH;
for i_neuron = 1:length(neurons) % 1
    if i_neuron>2&strcmp(neurons(i_neuron).fn(1:4), neurons(i_neuron-1).fn(1:4))~=1; animalCounter=animalCounter+1; end
    animalArray(i_neuron)=animalCounter; 
end

t = table(SCORE(:,1), condition',  animalArray', 'VariableNames', {'PC1', 'Type', 'Animal'});
writetable(t, 'UnmatchedPharmLMM.csv');

fprintf('\n Unmatched Data \n PC1 Variance Explained p=%0.2f', EXPLAINED(1)./sum(EXPLAINED));


SUL = find(t.Type==1|t.Type==2); tSUL = t(SUL, :);
PC1ModelSUL = fitglme(tSUL, 'PC1~Type+(1|Animal)'); PC1ModelSUL = anova(PC1ModelSUL);
fprintf('\n PC1 Sal vs SUL Model F=%2.1f, p=%0.2f', PC1ModelSUL.FStat(2), PC1ModelSUL.pValue(2));

SCH = find(t.Type==1|t.Type==3); tSCH = t(SCH, :);
PC1ModelSCH = fitglme(tSCH, 'PC1~Type+(Animal)'); PC1ModelSCH = anova(PC1ModelSCH);
fprintf('\n PC1 Sal vs SCH Model F=%2.1f, p=%0.2f', PC1ModelSCH.FStat(2), PC1ModelSCH.pValue(2));

