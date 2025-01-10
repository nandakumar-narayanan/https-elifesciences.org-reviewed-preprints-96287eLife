%-----------------------------------------------------------------------------------------------------------
% Written by Austin Bruce / Rewritten by Youngcho / checked by Kumar
%% Narayanan Lab
% load_neurons.m % Load

fprintf('\n**load_neurons.m will create neurons_tagging.mat')

load neurons_tagging.mat 

%%

for i_neuron = 1:length(neurons)
    if neurons(i_neuron).tagged == 1 && strcmp(neurons(i_neuron).type,'MSN')
        if strcmp(neurons(i_neuron).genotype, "D1-cre")
            neurons(i_neuron).celltype = 1; %D1-MSNs
            neurons(i_neuron).cellclass = "D1-MSN"; %D1-MSNs
        elseif strcmp(neurons(i_neuron).genotype, "D2-cre")
            neurons(i_neuron).celltype = 2; %D2-MSNs
            neurons(i_neuron).cellclass = "D2-MSN";
        end
    elseif neurons(i_neuron).tagged == 0 || ~strcmp(neurons(i_neuron).type,'MSN')
        neurons(i_neuron).celltype = 3; %Untagged
        neurons(i_neuron).cellclass = "uMSN";
    end

end

% Interval parameters
intStart =  -4;
trialStart = 0; 
trialEnd = 6; 
intEnd = 22;
binSize = 0.2;  % Holds 0.01 - 1 s bins
bw = 1.0;  % Holds 0.1 - 1.5;
interval = intStart:binSize:intEnd;
interval_binsShort = trialStart:binSize:9;
interval_bins = trialStart:binSize:trialEnd;

% For the Reviewer wanted a long interval 
intStartLong =  -11;
trialStartLong = 0; 
trialEndLong = 18; 
intEndLong = 22;
intervalLong = intStartLong:binSize:intEndLong;
interval_bins_Long = trialStartLong:binSize:trialEndLong;

 clear *PETH*
counter=1; clear tagPETH stagPETH tagCondition;
for i_neuron = 1:length(neurons) % 1
 warning off; 
        periEventSpike = peSpike(neurons(i_neuron).spikeTS, neurons(i_neuron).events.switchTrialInit, interval);
        periEventMotor = peSpike(neurons(i_neuron).events.shortPoke, neurons(i_neuron).events.switchTrialInit, interval);  % few switches in the first 6 seconds, and sometimes no nosepokes. 

        behFR(i_neuron) =    length(find(isnan(periEventSpike)))./[size(periEventSpike,1)*[intEnd-intStart]];

        tagPETH(i_neuron, :) = gksmooth(periEventSpike, interval, bw);  % Smooth for PCA


        % Makes a table for trial by trial glme
        a = histc(periEventSpike',interval_bins); % time x trial
        b = histc(periEventMotor',interval_bins); % time x trial
        spikes = reshape(a, [], 1); 
        motor = reshape(b, [], 1); 
        c = repmat(interval_bins, 1, size(a,2))';
        T_neuron = table(spikes, c, motor,'VariableNames', {'FiringRate', 'Time', 'Motor'});
        cr(i_neuron) = neurons(i_neuron).cr; % Correlation between waveforms; 

        try
        lmNeuron = fitglme(T_neuron, 'FiringRate~Time+Motor');  % Our GLM tries to account for both
        lmNeuronAnova = anova(lmNeuron);
        neurons(i_neuron).pTime2 = [lmNeuronAnova{2,5}];
        neurons(i_neuron).LMtime2 = lmNeuron;
        a = fixedEffects(lmNeuron); 
        neurons(i_neuron).slope = a(2);
        neurons(i_neuron).motorSlope = a(3);

        catch
        lmNeuron = fitglme(T_neuron, 'FiringRate~Time');  % In the case of poor fit, where there are few nosepokes, we revert to time. 
        lmNeuronAnova = anova(lmNeuron);
        neurons(i_neuron).pTime2 = [lmNeuronAnova{2,5}];
        neurons(i_neuron).LMtime2 = lmNeuron;
        a = fixedEffects(lmNeuron); 
        neurons(i_neuron).slope = a(2);
        neurons(i_neuron).motorSlope = NaN;
        end    
        % 
        % % Linear models to different epochs; first, the full interval
        periEventSpikeLong = peSpike(neurons(i_neuron).spikeTS, neurons(i_neuron).events.switchTrialInit, intervalLong);
        tagPETHLongInterval(i_neuron, :) = gksmooth(periEventSpikeLong, intervalLong, bw);  % Smooth for PCA


        a = histc(periEventSpikeLong',interval_bins_Long); % time x trial
        spikes = reshape(a, [], 1);
        c = repmat(interval_bins_Long, 1, size(a,2))';
        T_neuron = table(spikes, c, 'VariableNames', {'FiringRate', 'Time'});

        lmNeuron = fitglme(T_neuron, 'FiringRate~Time');
        lmNeuronAnova = anova(lmNeuron);
        neurons(i_neuron).pTime2 = [lmNeuronAnova{2,5}];
        neurons(i_neuron).LMtime2 = lmNeuron;
        a = fixedEffects(lmNeuron); 
        slopeLong(i_neuron) = a(2);

        % %Linear models to the time of the switch trial-by-trial; for a
        % % %reviewer
        spikes = []; times = []; motor = []; interval_bins_switch = []; 
        for j_trial = 1:size(periEventSpike,1)
            switchTime = neurons(i_neuron).events.switchDeparts(j_trial)-neurons(i_neuron).events.switchTrialInit(j_trial);
            interval_bins_switch = [0:binSize:switchTime];
            spikes = [spikes histc(periEventSpike(j_trial,:),interval_bins_switch)]; % time x trial
            motor = [motor histc(periEventMotor(j_trial,:),interval_bins_switch)]; % time x trial
            times =  [times interval_bins_switch]; 
        end
        T_neuron = table(spikes', motor', times','VariableNames', {'FiringRate', 'Motor', 'Time'});
        lmNeuron = fitglme(T_neuron, 'FiringRate~Time');  % In the case of poor fit, where there are few nosepokes, we revert to time.
        lmNeuronAnova = anova(lmNeuron);
        pTagTime(i_neuron) = [lmNeuronAnova{2,5}];
        a = fixedEffects(lmNeuron); 
        slopeSwitch(i_neuron) = a(2);

warning on; 
end

% Valid neurons based on average firing rate over the whole session
valid_idx = [neurons.avgFR]>0.5 & [neurons.avgFR] <20 & [neurons.celltype] < 3;
fprintf('total neurons for analysis %d\n',sum(valid_idx));

tagCondition = [neurons.celltype];
tagCondition = tagCondition(valid_idx);
slope = [neurons.slope];  
slope = slope(valid_idx);  
slopeLong = slopeLong(valid_idx);
slopeSwitch = slopeSwitch(valid_idx);

% For PCA
TimeInterval = find(interval>[trialStart-binSize]&interval<[trialEnd+binSize]);
ctagPETH = tagPETH(valid_idx, TimeInterval); 
cInterval = interval(TimeInterval);
ztagPETH = zscore(ctagPETH')';


earlyInterval = find(cInterval>0&cInterval<1); earlyIntervalZscore = mean(ctagPETH(:, earlyInterval)')'; 
lateInterval = find(cInterval>5&cInterval<6); lateIntervalZScore = mean(ctagPETH(:, lateInterval)')'; 

% Does PCA
warning off; [tagCOEFF, tagSCORE, LATENT, TSQUARED, tagEXPLAINED] = pca(ztagPETH); warning on;


neuronsTagged=neurons(valid_idx); 

% have to build an ensmemble for each animal - as animals predict trials!
animalCounter = 1;  clear animalArray;
counter=1; clear *PETH* sPETH condition FR Sex direction;
animalCounter = 1;  clear animalArray;

%% Loop through the neuron array and find all the animals, and give them a number (now they are strings)
for i_neuron = 1:length(neuronsTagged) % 1
    if i_neuron>2&strcmp(neuronsTagged(i_neuron).fn(1:8), neuronsTagged(i_neuron-1).fn(1:8))~=1; animalCounter=animalCounter+1; end
    animalArray(i_neuron)=animalCounter; 
    FR(i_neuron) = neuronsTagged(i_neuron).avgFR;
    % EL is switching L; ER is switching Right
    if contains(neuronsTagged(i_neuron).nameDB, 'EL'); direction(i_neuron)=0; else; direction(i_neuron)=1; end;

    % Gets the Sex of each animal
    Sex(i_neuron) =  GetAnimalSex(neuronsTagged(i_neuron).nameDB);
end

uniqueAnimals = [find(diff(animalArray)==1) length(animalArray)];
for i = 1:length(uniqueAnimals)
  % This the variable being predicted
            switchTimes = neuronsTagged(uniqueAnimals(i)).events.switchDeparts'-neuronsTagged(uniqueAnimals(i)).events.switchTrialInit;
            sTArray(i,:) = [mean(switchTimes) neuronsTagged(uniqueAnimals(i)).celltype]; % Switch Times per animal
            numTrials(i,:) = length(switchTimes); 
end

% Non-significant stats in the paper
motorSlope=[neurons.motorSlope];
motorSlope = motorSlope(valid_idx);  % Even those where the model could fit motor slopes, there was no difference
t = table(slope'./binSize, slopeSwitch'./binSize, tagSCORE(:,1), tagCondition',  animalArray', direction', Sex', FR',earlyIntervalZscore, lateIntervalZScore, motorSlope', 'VariableNames', {'Slope', 'SlopeSwitch', 'PC1', 'Type', 'Animal', 'direction', 'Sex', 'FR', 'Early', 'Late', 'Motor'});

%% Reviewers wanhted plots by animal
animals = unique(t.Animal);
clear data; 
figure(38); clf; 
subplot(121); 
for i_animal = 1:length(animals)
data{i_animal*3-2} = t.Slope(t.Animal==i_animal&t.Type==2); 
end
jitterPlot(data, [1 0 0]); set(gca, 'xtick', [1:3:27], 'xticklabel', [1:9]); ylim([-2.1 2]); box off; 
ylabel('Slope'); title('D2 MSNs'); xlabel('Animal #');

subplot(122); 
for i_animal = 1:length(animals)
data{i_animal*3-2} = t.Slope(t.Animal==i_animal&t.Type==1); 
end
jitterPlot(data, [0 0 1]); set(gca, 'xtick', [1:3:27], 'xticklabel', [1:9]); ylim([-2.1 2]); box off; 
ylabel('Slope'); title('D1 MSNs'); xlabel('Animal #');


clear data; 
figure(37); clf; 
subplot(121); 
for i_animal = 1:length(animals)
data{i_animal*3-2} = t.PC1(t.Animal==i_animal&t.Type==2); 
end
jitterPlot(data, [1 0 0]); set(gca, 'xtick', [1:3:27], 'xticklabel', [1:9]); ylim([-6 6]); box off; 
ylabel('PC1'); title('D2 MSNs'); xlabel('Animal #');

subplot(122); 
for i_animal = 1:length(animals)
data{i_animal*3-2} = t.PC1  (t.Animal==i_animal&t.Type==1); 
end
jitterPlot(data, [0 0 1]); set(gca, 'xtick', [1:3:27], 'xticklabel', [1:9]); ylim([-6 6]); box off; 
ylabel('PC1'); title('D1 MSNs'); xlabel('Animal #');
