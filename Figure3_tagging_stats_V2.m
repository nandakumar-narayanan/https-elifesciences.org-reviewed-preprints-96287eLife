


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


% Non-significant stats in the paper
motorSlope=[neurons.motorSlope];
motorSlope = motorSlope(valid_idx);  % Even those where the model could fit motor slopes, there was no difference

neuronsTagged=neurons(valid_idx); 
% have to build an ensmemble for each animal - as animals predict trials!
animalCounter = 1;  clear animalArray;
counter=1; clear *PETH* sPETH condition FR Sex direction;
animalCounter = 1;  clear animalArray;
minSwitchTime = []; 
%% Loop through the neuron array and find all the animals, and give them a number (now they are strings)
for i_neuron = 1:length(neuronsTagged) % 1
    if i_neuron>2&strcmp(neuronsTagged(i_neuron).fn(1:8), neuronsTagged(i_neuron-1).fn(1:8))~=1; animalCounter=animalCounter+1; end
    animalArray(i_neuron)=animalCounter; 
    FR(i_neuron) = neuronsTagged(i_neuron).avgFR;
    % EL is switching L; ER is switching Right
    if contains(neuronsTagged(i_neuron).nameDB, 'EL'); direction(i_neuron)=0; else; direction(i_neuron)=1; end;

    % Checks to see if a neuron has a switch before 6 seconds; 
    minSwitchTime(i_neuron) = min(neuronsTagged(i_neuron).events.switchDeparts-neuronsTagged(i_neuron).events.switchTrialInit');

    % Gets the Sex of each animal
    Sex(i_neuron) =  GetAnimalSex(neuronsTagged(i_neuron).nameDB);

    periEventSpike = peSpike(neuronsTagged(i_neuron).spikeTS, neuronsTagged(i_neuron).events.switchTrialInit, interval);
    tagPETH(i_neuron, :) = gksmooth(periEventSpike, interval, bw);  % Smooth for PCA
end

uniqueAnimals = [find(diff(animalArray)==1) length(animalArray)];

allSwitches = []; 
for i = 1:length(uniqueAnimals)
  % This the variable being predicted
            switchTimes = neuronsTagged(uniqueAnimals(i)).events.switchDeparts'-neuronsTagged(uniqueAnimals(i)).events.switchTrialInit;
            sTArray(i,:) = [mean(switchTimes) neuronsTagged(uniqueAnimals(i)).celltype]; % Switch Times per animal
            numTrials(i,:) = length(switchTimes); 
            allSwitches = [allSwitches switchTimes];
end

% For PCA
TimeInterval = find(interval>[trialStart-binSize]&interval<[trialEnd+binSize]);
ctagPETH = tagPETH(:, TimeInterval); 
cInterval = interval(TimeInterval);
ztagPETH = zscore(ctagPETH')';


earlyInterval = find(cInterval>0&cInterval<5); earlyIntervalZscore = mean(ctagPETH(:, earlyInterval)')'; 
lateInterval = find(cInterval>5&cInterval<6); lateIntervalZScore = mean(ctagPETH(:, lateInterval)')'; 

% Does PCA
warning off; [tagCOEFF, tagSCORE, LATENT, TSQUARED, tagEXPLAINED] = pca(ztagPETH); warning on;

t = table(slope'./binSize, slopeSwitch'./binSize, tagSCORE(:,1), tagCondition',  animalArray', direction', Sex', FR',earlyIntervalZscore, lateIntervalZScore, motorSlope', 'VariableNames', {'Slope', 'SlopeSwitch', 'PC1', 'Type', 'Animal', 'direction', 'Sex', 'FR', 'Early', 'Late', 'Motor'});
writetable(t, 'SlopeForLMM.csv');


clear data
data{1} = sTArray(find(sTArray(:,2)==1),1);
data{2} = sTArray(find(sTArray(:,2)==2),1);
fprintf('\n\n Switch Time data D2Cre %0.1f (%0.1f-%0.1f)', median(data{2}), quantile(data{2},0.25),quantile(data{2},0.75))
fprintf('\n Switch Time data D1Cre  %0.1f (%0.1f-%0.1f)', median(data{1}), quantile(data{1},0.25),quantile(data{1},0.75))
fprintf('\n Switch Time in D2Cre vs D1Cre mice p = %0.2f', ranksum(data{1}, data{2}) )

data{1} = numTrials(find(sTArray(:,2)==1),1);
data{2} = numTrials(find(sTArray(:,2)==2),1);
fprintf('\n # Trials data D2Cre %0.0f (%0.0f-%0.0f)', median(data{2}), quantile(data{2},0.25),quantile(data{2},0.75))
fprintf('\n # Trials data D1Cre  %0.0f (%0.0f-%0.0f)', median(data{1}), quantile(data{1},0.25),quantile(data{1},0.75))

D1FR = FR(tagCondition==1); D2FR = FR(tagCondition==2);
fprintf('\n\n D2 FR %0.1f (%0.1f-%0.1f)', median(D2FR), quantile(D2FR,0.25),quantile(D2FR,0.75));
fprintf('\n D1 FR %0.1f (%0.1f-%0.1f)', median(D1FR), quantile(D1FR,0.25),quantile(D1FR,0.75));
FRModel = fitglme(t, 'FR~Type+(1|Animal)'); FRModel = anova(FRModel);
fprintf('\n Differences in Firing Rate F=%2.1f, p=%0.2f\n\n', FRModel.FStat(2), FRModel.pValue(2))


EarlyModel = fitglme(t, 'Early~Type+(1|Animal)'); EarlyModel = anova(EarlyModel);
fprintf('\n\n D2MSN vs D1 MSN activity Early Activity Model F=%2.1f, p=%0.2f', EarlyModel.FStat(2), EarlyModel.pValue(2))

LateModel = fitglme(t, 'Late~Type+(1|Animal)'); LateModel = anova(LateModel);
fprintf('\n D2MSN vs D1 MSN activityLate Activity Model F=%2.1f, p=%0.2f', LateModel.FStat(2), LateModel.pValue(2))

%% PCA Analysis
figure(3); 
subplot(2,2,2);
cla; hold on; clear data; 
data{1} = tagSCORE(tagCondition==2,1);  %Switched so that the plot appears first
data{2} = tagSCORE(tagCondition==1,1);
jitterPlot(data, [1 0 0; 0 0 1]);
fprintf('\n\n\n PC1 Variance Explained %1.0f', tagEXPLAINED(1)./sum(tagEXPLAINED)*100);
load tagRAND.mat; 
fprintf('\n Variance in Random Data %1.0f (%1.0f-%1.0f)', median(rEXPLAINED_array(:,1)), quantile((rEXPLAINED_array(:,1)), 0.25), quantile((rEXPLAINED_array(:,1)), 0.75))
fprintf('\n PC1 Variance Explained vs Random Data %0.3f', length(find(tagEXPLAINED(1)<sort(rEXPLAINED_array(:,1))))./size(rEXPLAINED_array,1));

fprintf('\n\n PC1 data D2s  %0.1f (%0.1f-%0.1f)', median(data{1}), quantile(data{1},0.25),quantile(data{1},0.75))
fprintf('\n PC1 data D1s %0.1f (%0.1f-%0.1f)', median(data{2}), quantile(data{2},0.25),quantile(data{2},0.75))
PC1Model = fitglme(t, 'PC1~Type+(1|Animal)'); PC1Model = anova(PC1Model);
fprintf('\n PC1Model F=%2.1f, p=%0.3f', PC1Model.FStat(2), PC1Model.pValue(2))
pwr = sampsizepwr('t2',[mean(data{1}) std(data{1})],mean(data{2}),[],length(data{1})); 
SexModel = fitglme(t, 'PC1~Sex+(1|Animal)');
SexModel = anova(SexModel);
fprintf('\n Effects of Sex F=%2.1f, p=%0.2f', SexModel.FStat(2), SexModel.pValue(2))

DirectionModel = fitglme(t, 'PC1~direction+(1|Animal)');
DirectionModel = anova(DirectionModel);
fprintf('\n Effects of Direction F=%2.1f, p=%0.2f', DirectionModel.FStat(2), DirectionModel.pValue(2))

fprintf('\n PC1: cohen d %2.1f, D1 vs D2 power = %0.2f', cohend(data{1}, data{2}), pwr )
fprintf('\n PC1, which ramps down, is more than 0 for D1-MSNs: signrank p=%0.3f', signrank(tagSCORE(tagCondition==1,1)', zeros(1,length(find(tagCondition==1)))));
fprintf('\n PC1, which ramps down, is less than 0 for D2-MSNs: signrank p=%0.3f', signrank(tagSCORE(tagCondition==2,1)', zeros(1,length(find(tagCondition==2)))));

groups = [ones(1,length(data{1})) 2*ones(1,length(data{2}))];
data = [data{1}; data{2}]';
boxplot(data, groups);
set(gca, 'xtick', [1:2], 'xticklabel', {'D2'; 'D1'}); 
ylabel('PC1 Score'); 




%% Slope Analysis
subplot(2,2,4); 
cla; hold on; clear data; 
data{1} = slope(tagCondition==2)./binSize; %Switched so that the plot appears first
data{2} = slope(tagCondition==1)./binSize;

jitterPlot(data, [1 0 0; 0 0 1]);
set(gca, 'xtick', [1:2], 'xticklabel', {'D2'; 'D1'}); 
fprintf('\n\n D2 Slope %0.2f (%0.2f-%0.2f)', median(data{1}), quantile(data{1},0.25),quantile(data{1},0.75))
fprintf('\n D1 Slope %0.2f (%0.2f-%0.2f)', median(data{2}), quantile(data{2},0.25),quantile(data{2},0.75))

pwr =  sampsizepwr('t2',[mean(data{1}) std(data{1})], mean(data{2}),[],32); 
fprintf('\ncohen d %2.1f, Slope: D1 vs D2 power = %0.2f', cohend(data{1}, data{2}), pwr )

groups = [ones(1,length(data{1})) 2*ones(1,length(data{2}))];
data = [data{1} data{2}]';
boxplot(data, groups); 
ylabel('Slope'); set(gca, 'xtick', [1:2], 'xticklabel', {'D2'; 'D1'}); 

subplot(2,2,1); cla;  hold on; 
plot(cInterval, tagCOEFF(:,1));
set(gca,'xtick', [0 6 18]); 
box off; 
set(gcf, 'Color', 'White'); 
set(gca, 'ytick', []); 


subplot(2,2,3); cla;  hold on; 
pH = plot(tagEXPLAINED./sum(tagEXPLAINED)*100, 'ko-'); ylabel('Variance Explained'); xlabel('Principal Component'); 
box off; xlim([0.5 10]); set(pH, 'MarkerFaceColor', [0 0 0], 'MarkerSize', 10); 

slopeModel = fitglme(t, 'Slope~Type+(1|Animal)'); slopeModel = anova(slopeModel);
fprintf('\n SlopeModel F=%2.1f, p=%0.3f', slopeModel.FStat(2), slopeModel.pValue(2))

SexModel = fitglme(t, 'Slope~Sex+(1|Animal)');
SexModel = anova(SexModel);
fprintf('\n Effects of Sex F=%2.1f, p=%0.2f', SexModel.FStat(2), SexModel.pValue(2))

DirectionModel = fitglme(t, 'Slope~direction+(1|Animal)');
DirectionModel = anova(DirectionModel);
fprintf('\n Effects of Direction F=%2.1f, p=%0.2f', DirectionModel.FStat(2), DirectionModel.pValue(2))

nonOutlier = find(t.Slope<quantile(t.Slope, 0.975)&t.Slope>quantile(t.Slope, 0.025));
toutliersremoved = t(nonOutlier, :);

slopeModelNoOutliers = fitglme(toutliersremoved, 'Slope~Type+(1|Animal)'); slopeModelNoOutliers = anova(slopeModelNoOutliers);
fprintf('\n SlopeModel No Outliers F=%2.1f, p=%0.3f', slopeModelNoOutliers.FStat(2), slopeModelNoOutliers.pValue(2));

nonOutlier = find(t.Slope>-1.3);
toutliersremoved = t(nonOutlier, :);

slopeModelNoOutliers = fitglme(toutliersremoved, 'Slope~Type+(1|Animal)'); slopeModelNoOutliers = anova(slopeModelNoOutliers);
fprintf('\n SlopeModel Just D1 Outliers removed F=%2.1f, p=%0.2f', slopeModelNoOutliers.FStat(2), slopeModelNoOutliers.pValue(2));

slopeSwitchModel = fitglme(t, 'SlopeSwitch~Type+(1|Animal)'); slopeSwitchModel = anova(slopeSwitchModel);
fprintf('\n SlopeModel F=%2.1f, p=%0.2f', slopeSwitchModel.FStat(2), slopeSwitchModel.pValue(2))

[r, p] = corr(t.Slope, t.PC1); 
fprintf('\nCorrelation Between PC1 and Slope %0.2f %e and %2.0f variance', r, p, 100*r.^2)


MotorModel = fitglme(t, 'Motor~Type+(1|Animal)'); MotorModel = anova(MotorModel);
fprintf('\n\n Motor Model F=%2.1f, p=%0.2f\n\n', MotorModel.FStat(2), MotorModel.pValue(2))
%
fprintf('\n\n Switches Less than 6 seconds %2.0f percent', length(find(allSwitches<6))./length(allSwitches)*100);

fprintf('\n\n')