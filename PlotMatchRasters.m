
%-----------------------------------------------------------------------------------------------------------
% Written by Austin Bruce / Rewritten by Youngcho / checked by Kumar
%% Narayanan Lab% 
% % %%%%%% PHARM
% Loading code
% tic; load_neurons; toc; % Takes about 12 minutes
% for i = 1:length(neuronDB)     
%     trialLength(i) = length(neuronDB(i).events.switchTrialInit); 
%     numTimeStamps(i) = [neuronDB(i).numTS];
%     FR(i) = length(neuronDB(i).spikeTS)./max(neuronDB(i).spikeTS);
% end
% neuronDB = neuronDB(find(FR>0.5&FR<20)); FR=FR(find(FR>0.5&FR<20)); % drops neurons < 0.5 and >20 Hz;
% 
% save PharmNeuronDB neuronDB; % 456 neurons
% clear; 
% 
% load PharmNeuronDB.mat

addpath(genpath('./analytical_tools_striatalcircuits/'))
load neurons_pharm.mat 
%%

avgFR = [neurons.avgFR];
num_spike = [neurons.numTS]; 
% valid_idx = avgFR>0.5 & avgFR<20 & strcmp({neurons.type},'MSN') & num_spike > 5000;
valid_idx = avgFR>0.5 & avgFR<20 & strcmp({neurons.type},'MSN');
fprintf('total neurons for analysis %d\n',sum(valid_idx))

%%

intStart =  -4;
trialStart = 0; 
trialEnd = 6; 
intEnd = 22;
binSize = 0.25;  % Holds 0.01 - 1 s bins
bw = 1.0;  % Holds 0.1 - 1.5;
interval = intStart:binSize:intEnd;
interval_bins = trialStart:binSize:18;

mBinSize =0.1; 
motorInterval = -1.5:mBinSize:1.5; mStart=-1; mEnd=1; 

counter=1; clear *PETH* sPETH condition;
for i_neuron = 1:length(neurons) % 1
    periEventSpike = peSpike(neurons(i_neuron).spikeTS, neurons(i_neuron).events.switchTrialInit, interval);

    periEventSpikeSwitch = peSpike(neurons(i_neuron).spikeTS, neurons(i_neuron).events.switchDeparts, motorInterval);

    allPokes = [neurons(i_neuron).events.rightPoke; neurons(i_neuron).events.leftPoke];
    periEventSpikePoke = peSpike(neurons(i_neuron).spikeTS, allPokes, motorInterval);
 
    
    PETH(counter, :) = gksmooth(periEventSpike, interval, bw);  % Smooth for PCA
 
    %FR(i_neuron) = length(neuronDB(i_neuron).spikeTS)./max(neuronDB(i_neuron).spikeTS);
    try
      a=fixedEffects(neurons(i_neuron).LMtime);
        slope(i_neuron) = a(2);  
    catch 
        slope(i_neuron)=NaN; 
    end
      
    behFR(i_neuron) =    length(find(isnan(periEventSpike)))./[size(periEventSpike,1)*[intEnd-intStart]];
    
    condition(counter) = neurons(i_neuron).condition;
    counter=counter+1;
end


condition = [neurons.condition];
condition = condition(valid_idx);


TimeInterval = find(interval>[trialStart-binSize]&interval<[trialEnd+binSize]);


cPETH = PETH(:, TimeInterval);  
cInterval=interval(TimeInterval); 
zPETH = zscore(cPETH')';


warning off;
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(zPETH,'NumComponents',3);

warning on;


figure(6); clf;  set(gcf, 'Color', 'White');
subplot(3,3,1);
PCtoSort = 1; 
PETHtoPlot = cPETH(condition==1,:); 
[~,sortKey] = sort(SCORE(condition==1, PCtoSort));  % For sorting by PCA
imagesc(cInterval, [], zscore(PETHtoPlot(sortKey,:)')', [-3 3]);  ylabel('Early Ramping Neurons (#)');
xlabel('Time from Start (Seconds)'); title('Saline');
set(gca, 'xtick', [0 6 18], 'ytick', [1 size(PETHtoPlot,1)]);
colorbar; colormap('jet');

subplot(3,3,2);
PETHtoPlot = cPETH(condition==2,:);
[~,sortKey] = sort(SCORE(condition==2,PCtoSort));  % For sorting by PCA
imagesc(cInterval, [], zscore(PETHtoPlot(sortKey,:)')', [-3 3]);  ylabel('Early Ramping Neurons (#)');
xlabel('Time from Start (Seconds)'); title('D2Block');
set(gca, 'xtick', [0 6 18], 'ytick', [1 size(PETHtoPlot,1)]);
colorbar; colormap('jet');

subplot(3,3,3);
PETHtoPlot = cPETH(condition==3,:);
[~,sortKey] = sort(SCORE(condition==3, PCtoSort));  % For sorting by PCA
imagesc(cInterval, [], zscore(PETHtoPlot(sortKey,:)')', [-3 3]);  ylabel('Early Ramping Neurons (#)');
xlabel('Time from Start (Seconds)'); title('D1Block');
set(gca, 'xtick', [0 6 18], 'ytick', [1 size(PETHtoPlot,1)]);
colorbar; colormap('jet');




%%% Refocus on first 6 seconds
neuronsFull = neurons(valid_idx); 
intStart =  -4;
trialStart = 0; 
trialEnd = 6; 
intEnd = 22;
binSize = 0.25;  % Holds 0.01 - 1 s bins
bw = 1.0;  % Holds 0.1 - 1.5;
interval = intStart:binSize:intEnd;

mBinSize =0.1; 

counter=1; clear *PETH* sPETH condition;
for i_neuron = 1:length(neuronsFull) % 1
    periEventSpike = peSpike(neuronsFull(i_neuron).spikeTS, neuronsFull(i_neuron).events.switchTrialInit, interval);
    PETH(counter, :) = gksmooth(periEventSpike, interval, bw);  % Smooth for PCA
    condition(counter) = neuronsFull(i_neuron).condition;
    counter=counter+1;

    behFR(i_neuron) =    length(find(isnan(periEventSpike)))./[size(periEventSpike,1)*[intEnd-intStart]];

end


TimeInterval = find(interval>[trialStart-binSize]&interval<[trialEnd+binSize]);

cPETH = PETH(:, TimeInterval);  
cInterval=interval(TimeInterval); 
zPETH = zscore(cPETH')';


warning off;
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(zPETH,'NumComponents',3);
warning on;


subplot(3,3,4); cla; hold on; 
plot(interval(TimeInterval), COEFF(:,1)+0.3, 'k', 'LineWidth', 7); 
plot(interval(TimeInterval), COEFF(:,2)-0.3, 'k', 'LineWidth', 7); 
set(gca, 'ytick', [])
xlabel('Time from Trial Start (seconds)')
axis tight; 

subplot(3,3,5); cla; hold on; 
plot(EXPLAINED./sum(EXPLAINED), 'ko-', 'MarkerFaceColor','k','MarkerSize', 10)
xlim([0 10]);
ylabel('Fraction Variance Explained')
xlabel('Component #')



subplot(3,3,7); cla; hold on; 
vdata{1} = (SCORE(condition==1,1)); 
vdata{2} = (SCORE(condition==2,1)); 
vdata{3} = (SCORE(condition==3,1)); 
data = [vdata{1}; vdata{2}; vdata{3}]';
groups = [ones(1,length(vdata{1})) 2*ones(1,length(vdata{2})) 3*ones(1,length(vdata{3})) ];

scatter(ones(size(vdata{1})).*(1+(rand(size(vdata{1}))-0.5)/10),vdata{1},'k','filled')
scatter(2*ones(size(vdata{2})).*(1+(rand(size(vdata{2}))-0.5)/10),vdata{2},'r','filled')
scatter(3*ones(size(vdata{3})).*(1+(rand(size(vdata{3}))-0.5)/10),vdata{3},'b','filled')
boxplot(data, groups); box off; 


fprintf('\n *****')
fprintf('\n PC1 data saline %0.2f (%0.2f-%0.2f)', median(vdata{1}), quantile(vdata{1},0.25),quantile(vdata{1},0.75))
fprintf('\n PC1 data D2 %0.2f (%0.2f-%0.2f)', median(vdata{2}), quantile(vdata{2},0.25),quantile(vdata{2},0.75))
fprintf('\n PC1 data D1 %0.2f (%0.2f-%0.2f)', median(vdata{3}), quantile(vdata{3},0.25),quantile(vdata{3},0.75))

fprintf('\n PC1: signrank Saline V. D2 Block p = %0.3f, cohen d %2.1f', ranksum(SCORE(condition==1,1), SCORE(condition==3,1)), cohend(SCORE(condition==1,1), SCORE(condition==3,1)))
fprintf('\n PC1: signrank Saline V. D1 Block p = %0.3f, cohen d %2.1f', ranksum(SCORE(condition==1,1), SCORE(condition==2,1)), cohend(SCORE(condition==1,1), SCORE(condition==2,1)))
fprintf('\n PC1: signrank  D2 Block V. D1 Block p = %0.3f, cohen d %2.1f', ranksum(SCORE(condition==3,1), SCORE(condition==2,1)), cohend(SCORE(condition==3,1), SCORE(condition==2,1)))
csvwrite('./stats/PC1Saline.csv', (SCORE(condition==1,1)));
csvwrite('./stats/PC1D2.csv', (SCORE(condition==2,1)));
csvwrite('./stats/PC1D1.csv', (SCORE(condition==3,1)));



ylabel('PC1 Score'); set(gca, 'xtick', [1:3], 'xticklabel', {'Saline'; 'D2Block'; 'D1Block'});
box off; 

subplot(3,3,8); cla; hold on; 
vdata{1} = (SCORE(condition==1,2)); 
vdata{2} = (SCORE(condition==2,2)); 
vdata{3} = (SCORE(condition==3,2)); 
data = [vdata{1}; vdata{2}; vdata{3}]';
groups = [ones(1,length(vdata{1})) 2*ones(1,length(vdata{2})) 3*ones(1,length(vdata{3})) ];
boxplot(data, groups, 'k'); 

scatter(ones(size(vdata{1})).*(1+(rand(size(vdata{1}))-0.5)/10),vdata{1},'k','filled')
scatter(2*ones(size(vdata{2})).*(1+(rand(size(vdata{2}))-0.5)/10),vdata{2},'r','filled')
scatter(3*ones(size(vdata{3})).*(1+(rand(size(vdata{3}))-0.5)/10),vdata{3},'b','filled')
ylabel('PC2 Score'); set(gca, 'xtick', [1:3], 'xticklabel', {'Saline'; 'D2Block'; 'D1Block'});
box off; 


fprintf('\n *****')
fprintf('\n PC2: signrank Saline V. D2 Block p = %0.3f, cohen d %2.1f', ranksum(SCORE(condition==1,2), SCORE(condition==3,2)), cohend(SCORE(condition==1,2), SCORE(condition==3,2)))
fprintf('\n PC2: signrank Saline V. D1 Block p = %0.3f, cohen d %2.1f', ranksum(SCORE(condition==1,2), SCORE(condition==2,2)), cohend(SCORE(condition==1,2), SCORE(condition==2,2)))
fprintf('\n PC2: signrank  D2 Block V. D1 Block p = %0.3f, cohen d %2.1f', ranksum(SCORE(condition==3,2), SCORE(condition==2,2)), cohend(SCORE(condition==3,2), SCORE(condition==2,2)))
csvwrite('./stats/PC2Saline.csv', (SCORE(condition==1,2)));
csvwrite('./stats/PC2D2.csv', (SCORE(condition==2,2)));
csvwrite('./stats/PC2D1.csv', (SCORE(condition==3,2)));


neuronsPharm=neurons(valid_idx);


% have to build an ensmemble for each animal - as animals predict trials!
animalCounter = 1;  clear animalArray;
counter=1; clear *PETH* sPETH;
for i_neuron = 1:length(neuronsPharm) % 1
    if i_neuron>2&strcmp(neuronsPharm(i_neuron).fn(1:4), neuronsPharm(i_neuron-1).fn(1:4))~=1; animalCounter=animalCounter+1; end
    animalArray(i_neuron)=animalCounter; 

    FR(i_neuron) = neuronsPharm(i_neuron).avgFR;
    % EL is switching L; ER is switching Right
    if contains(neuronsPharm(i_neuron).nameDB, 'EL'); direction(i_neuron)=0; else; direction(i_neuron)=1; end;

    % Gets the Sex of each animal
    Sex(i_neuron) =  GetAnimalSex(neuronsPharm(i_neuron).nameDB);

end


t = table(SCORE(:,1), condition',  animalArray', direction', Sex', FR', 'VariableNames', {'PC1', 'Type', 'Animal', 'direction', 'Sex', 'FR'});
a = fitglme(t, 'PC1~Type+(1|Animal)');
writetable(t, 'SlopePharmLMM.csv');


SFR  = FR(condition==1);  D2FR = FR(condition==2); D1FR = FR(condition==3);
fprintf('\n S FR %0.2f (%0.2f-%0.2f)', median(SFR), quantile(SFR,0.25),quantile(SFR,0.75));
fprintf('\n D2 FR %0.2f (%0.2f-%0.2f)', median(D2FR), quantile(D2FR,0.25),quantile(D2FR,0.75));
fprintf('\n D1 FR %0.2f (%0.2f-%0.2f)', median(D1FR), quantile(D1FR,0.25),quantile(D1FR,0.75));



