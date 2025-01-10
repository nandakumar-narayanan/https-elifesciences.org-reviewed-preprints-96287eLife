
fprintf('\n**load_neurons.m will create neurons_pharm.mat')

load neuronDB_matched; 

% key variable is neurons - a DB of all neurons
% neurons.condition gives you the session; saline = 1; D
% %D2block with sulpiride = 2; 
% D1 block with SCH = 3;
% 99 neurons - each neurons should be aligned 

condition = [neurons.condition]; 
salIndex = find(condition==1);
sulIndex = find(condition==2);
schIndex = find(condition==3);

% Interval parameters
intStart =  -4;
trialStart = 0; 
trialEnd = 6; 
intEnd = 22;
binSize = 0.2;  
bw = 1;  
interval = intStart:binSize:intEnd;
interval_bins = trialStart:binSize:trialEnd;

% Makes PETHS for PCA
counter=1; clear *PETH* sPETH condition;
for i_neuron = 1:length(neurons) % 1
    periEventSpike = peSpike(neurons(i_neuron).spikeTS, neurons(i_neuron).events.switchTrialInit, interval);
    PETH(counter, :) = gksmooth(periEventSpike, interval, bw);  % Smooth for PCA;     % Makes PETHs
    
    % Makes a table for trial by trial glme
    a = histc(periEventSpike',interval_bins); % time x trial
    spikes = reshape(a, [], 1); 
    c = repmat(interval_bins, 1, size(a,2))';
    T_neuron = table(spikes, c,'VariableNames', {'FiringRate', 'Time'});
    lmNeuron = fitglme(T_neuron, 'FiringRate~Time');  % In the case of poor fit, where there are few nosepokes, we revert to time. 
    lmNeuronAnova = anova(lmNeuron);
    neurons(i_neuron).pTime2 = [lmNeuronAnova{2,5}];
    neurons(i_neuron).LMtime2 = lmNeuron;
    a = fixedEffects(lmNeuron); 
    slope(counter) = a(2);

    counter=counter+1; 
end


%% Find the time interval for each PETH --- first 6 seconds
% get the firing rate 1 second before start and 1 second after 18 s trial
TimeIntervalLong = ...
    find((interval> trialStart-binSize) & (interval<trialEnd+binSize));
cPETH_all = PETH(:, TimeIntervalLong);  
cInterval_all=interval(TimeIntervalLong); 
zPETH_all= zscore(cPETH_all')';

warning off; [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(zPETH_all,'NumComponents',3); warning on; 

[~,sortIndex] = sort(SCORE(salIndex,1)); % for PCA sorting; 

figure(6); CLIM=[-3 3]; set(gcf, 'Color', 'white'); colormap('jet')
subplot(231); imagesc(cInterval_all,[], zPETH_all(salIndex(sortIndex),:),CLIM); title('Saline'); set(gca,'xtick', [0 6 18]);  
ylabel('MSNs'); xlabel('Time from Trial  (seconds)'); colorbar;  
subplot(232); imagesc(cInterval_all,[], zPETH_all(sulIndex(sortIndex),:),CLIM); title('Sul'); set(gca,'xtick', [0 6 18]);  
ylabel('Matched MSNs'); xlabel('Time from  Start (seconds)'); colorbar;  
subplot(233); imagesc(cInterval_all,[], zPETH_all(schIndex(sortIndex),:),CLIM); title('SCH'); set(gca,'xtick', [0 6 18]);  
ylabel('Matched MSNs'); xlabel('Time from  Start (seconds)'); colorbar;  

subplot(234);
plot(cInterval_all, COEFF(:,1), 'k', 'LineWidth', 7); 
set(gca, 'ytick', []); xlabel('Time from  Start (seconds)')
axis tight; box off; colormap('jet')

subplot(235); cla; hold on; 
plot(EXPLAINED./sum(EXPLAINED), 'ko-', 'MarkerFaceColor','k','MarkerSize', 10)
xlim([0 10]);
ylabel('Fraction Variance Explained')
xlabel('Component #')



condition = [neurons.condition]; 

% Quantify how similar the waveforms are for matched sorting
counter = 1; 
for i=1:99;
    matchedSulWF(i) = corr(neurons(i).wf.mean', neurons(i+99).wf.mean'); 
    matchedSchWF(i) = corr(neurons(i).wf.mean', neurons(i+198).wf.mean'); 
    for j = 1:99; 
           if i~=j; 
            unmatchedWF(counter) = corr(neurons(i).wf.mean', neurons(j).wf.mean'); 
            counter = counter+1; 
           end
     end

end

fprintf('\n Correlation: Waveforms vs. Sulpiride %0.2f (%0.2f-%0.2f), vs unmatched WF %e ',median(matchedSulWF), quantile(matchedSulWF,0.25),quantile(matchedSulWF,0.75),ranksum(matchedSulWF,unmatchedWF ));
fprintf('\n Correlation: Waveforms vs. Sch 23390 %0.2f (%0.2f-%0.2f) vs unmatched WF %e ',median(matchedSchWF), quantile(matchedSchWF,0.25),quantile(matchedSchWF,0.75),  ranksum(matchedSchWF,unmatchedWF));

% have to build an ensmemble for each animal
animalCounter = 1; 
sessionCounter = 1;  clear sessionArray Sex FR direction animalArray;
counter=1; clear *PETH* sPETH ; 
for i_neuron = 1:length(neurons) % 1
    if i_neuron>2&strcmp(neurons(i_neuron).fn(1:4), neurons(i_neuron-1).fn(1:4))~=1; sessionCounter=sessionCounter+1; end
    sessionArray(i_neuron)=sessionCounter;
    if sessionCounter<12; animalArray(i_neuron)=sessionCounter;
    elseif sessionCounter>11&sessionCounter<23 animalArray(i_neuron)=sessionCounter-11;
    elseif sessionCounter>22; animalArray(i_neuron)=sessionCounter-22;
    end
        
    FR(i_neuron) = neurons(i_neuron).avgFR;
    % EL is switching L; ER is switching Right
    if contains(neurons(i_neuron).nameDB, 'EL'); direction(i_neuron)=0; else; direction(i_neuron)=1; end;
    % Gets the Sex of each session
    Sex(i_neuron) =  GetAnimalSex(neurons(i_neuron).nameDB);
end
n = [1:99 1:99 1:99]'; % Index for neurons

% Linear models
t = table(SCORE(:,1), condition',  animalArray', direction', Sex', FR', n,slope', 'VariableNames', {'PC1', 'Type', 'Animal', 'direction', 'Sex', 'FR', 'Neuron', 'Slope'});
writetable(t, 'PCAPharm.csv');
fprintf('\n\n Linear models - all account for between neuron variance');
SUL = find(t.Type==1|t.Type==2); tSUL = t(SUL, :);
SULSCH = find(t.Type==3|t.Type==2); tSULSCH = t(SULSCH, :);
SCH = find(t.Type==1|t.Type==3); tSCH = t(SCH, :);

SFR  = FR(condition==1);  D2FR = FR(condition==2); D1FR = FR(condition==3);
fprintf('\n SAL FR %0.1f (%0.1f-%0.1f)', median(SFR), quantile(SFR,0.25),quantile(SFR,0.75));
fprintf('\n D2 FR %0.1f (%0.1f-%0.1f)', median(D2FR), quantile(D2FR,0.25),quantile(D2FR,0.75));

FRModel = fitglme(tSUL, 'FR~Type+(1|Neuron)'); FRModel = anova(FRModel);
fprintf('\n SUL FR Model F=%2.1f, p=%0.2f', FRModel.FStat(2), FRModel.pValue(2));

fprintf('\n\n D1 FR %0.1f (%0.1f-%0.1f)', median(D1FR), quantile(D1FR,0.25),quantile(D1FR,0.75));
FRModel = fitglme(tSCH, 'FR~Type+(1|Neuron)'); FRModel = anova(FRModel);
fprintf('\n SCH FR Model F=%2.1f, p=%0.2f', FRModel.FStat(2), FRModel.pValue(2));


subplot(236); cla; hold on; 
vdata{1} = (SCORE(salIndex,1)); 
vdata{3} = (SCORE(schIndex,1)); 
vdata{2} = (SCORE(sulIndex,1)); 
data = [vdata{1}; vdata{2}; vdata{3}]';
groups = [ones(1,length(vdata{1})) 2*ones(1,length(vdata{2})) 3*ones(1,length(vdata{3})) ];

scatter(ones(size(vdata{1})).*(1+(rand(size(vdata{1}))-0.5)/10),vdata{1},'k','filled')
scatter(2*ones(size(vdata{2})).*(1+(rand(size(vdata{2}))-0.5)/10),vdata{2},'r','filled')
scatter(3*ones(size(vdata{3})).*(1+(rand(size(vdata{3}))-0.5)/10),vdata{3},'b','filled')
boxplot(data, groups); box off; 

fprintf('\n PCA*****')
fprintf('\n PC1 Variance Explained: %2.0f', 100*EXPLAINED(1)./sum(EXPLAINED))

load pharmRAND.mat; 

fprintf('\n PC1 Variance Explained in Random Data %0.03f', length(find(EXPLAINED(1)<sort(rEXPLAINED_array(:,1))))./size(rEXPLAINED_array,1));
fprintf('\n Variance in Random Data %1.0f (%1.0f-%1.0f)', median(rEXPLAINED_array(:,1)), quantile((rEXPLAINED_array(:,1)), 0.25), quantile((rEXPLAINED_array(:,1)), 0.75))

fprintf('\n\n PC1 data saline %0.1f (%0.1f-%0.1f)', median(vdata{1}), quantile(vdata{1},0.25),quantile(vdata{1},0.75))
fprintf('\n PC1 data D2 %0.1f (%0.1f-%0.1f)', median(vdata{2}), quantile(vdata{2},0.25),quantile(vdata{2},0.75))
fprintf('\n PC1 data D1 %0.1f (%0.1f-%0.1f)', median(vdata{3}), quantile(vdata{3},0.25),quantile(vdata{3},0.75))

PC1ModelSUL = fitglme(tSUL, 'PC1~Type+(1|Neuron)'); PC1ModelSUL = anova(PC1ModelSUL);
fprintf('\n D2 PCA Model F=%2.1f, p=%0.2f', PC1ModelSUL.FStat(2), PC1ModelSUL.pValue(2));


PC1ModelSCH = fitglme(tSCH, 'PC1~Type+(1|Neuron)'); PC1ModelSCH = anova(PC1ModelSCH);
fprintf('\n D1 PCA Model F=%2.1f, p=%0.2f', PC1ModelSCH.FStat(2), PC1ModelSCH.pValue(2));


PC1ModelSULSCH = fitglme(tSULSCH, 'PC1~Type+(1|Neuron)'); PC1ModelSULSCH = anova(PC1ModelSULSCH);
fprintf('\n D2 vs D1 PCA Model F=%2.1f, p=%0.2f', PC1ModelSULSCH.FStat(2), PC1ModelSULSCH.pValue(2));



SexDirectionModel = fitglme(tSUL, 'PC1~Sex+(1|Neuron)'); SexDirectionModel = anova(SexDirectionModel);
fprintf('\n SUL: Effects of Sex F=%2.1f, p=%0.2f', SexDirectionModel.FStat(2), SexDirectionModel.pValue(2))

SexDirectionModel = fitglme(tSUL, 'PC1~direction+(1|Neuron)'); SexDirectionModel = anova(SexDirectionModel);
fprintf('\n SUL: Effects of Direction F=%2.1f, p=%0.2f', SexDirectionModel.FStat(2), SexDirectionModel.pValue(2))

SexDirectionModel = fitglme(tSCH, 'PC1~Sex+(1|Neuron)'); SexDirectionModel = anova(SexDirectionModel);
fprintf('\n SCH: Effects of Sex F=%2.1f, p=%0.2f', SexDirectionModel.FStat(2), SexDirectionModel.pValue(2))

SexDirectionModel = fitglme(tSCH, 'PC1~direction+(1|Neuron)'); SexDirectionModel = anova(SexDirectionModel);
fprintf('\n SCH: Effects of Direction F=%2.1f, p=%0.2f', SexDirectionModel.FStat(2), SexDirectionModel.pValue(2))


uniquesessions = [find(diff(sessionArray)==1) length(sessionArray)];
for i = 1:length(uniquesessions)
  % This the variable being predicted
            switchTimes = neurons(uniquesessions(i)).events.switchDeparts'-neurons(uniquesessions(i)).events.switchTrialInit;
            numTrials(i,:) = length(switchTimes); 
end
fprintf('\n# Saline: Trials %0.0f (%0.0f-%0.0f)', median(numTrials(1:11)), quantile(numTrials(1:11),0.25),quantile(numTrials(1:11),0.75))
fprintf('\n# D2Blockade: Trials %0.0f (%0.0f-%0.0f)', median(numTrials(12:22)), quantile(numTrials(12:22),0.25),quantile(numTrials(12:22),0.75))
fprintf('\n# D1Blockade: Trials %0.0f (%0.0f-%0.0f)', median(numTrials(23:33)), quantile(numTrials(23:33),0.25),quantile(numTrials(23:33),0.75))


[r, p] = corr(t.Slope(1:99), t.PC1(1:99)); % Only for saline session
fprintf('\nCorrelation Between PC1 and Slope %0.2f %e and %2.0f variance', r, p, 100*r.^2)





