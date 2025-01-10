%-----------------------------------------------------------------------------------------------------------
% Written by Austin Bruce / Rewritten by Youngcho / checked by Kumar
%% Narayanan Lab
% load_neurons.m % Load

fprintf('\n**load_neurons.m will create neurons_tagging.mat\n')
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


earlyInterval = find(cInterval>0&cInterval<5); earlyIntervalZscore = mean(ctagPETH(:, earlyInterval)')'; 
lateInterval = find(cInterval>5&cInterval<6); lateIntervalZScore = mean(ctagPETH(:, lateInterval)')'; 

% Does PCA
warning off; [tagCOEFF, tagSCORE, LATENT, TSQUARED, tagEXPLAINED] = pca(ztagPETH); warning on;


% Plots Figure 2 - PETHS and averages
figure(2); clf;  set(gcf, 'Color', 'White');
subplot(2,2,1);
tagPETHtoPlot = ctagPETH(tagCondition==2,:); PCtoSort = 1; 

[~,sortKey] = sort(tagSCORE(tagCondition==2, PCtoSort));  % For sorting by PCA
imagesc(cInterval, [], zscore(tagPETHtoPlot(sortKey,:)')', [-3 3]);  ylabel('TAGGED D2 MSNs  (#)');
xlabel('Time from Start (Seconds)'); 
set(gca, 'xtick', [0 6 18], 'ytick', [1 size(tagPETHtoPlot,1)]);
colorbar;  colormap('jet');

subplot(2,2,2);
tagPETHtoPlot = ctagPETH(tagCondition==1,:); PCtoSort = 1; 
[~,sortKey] = sort(tagSCORE(tagCondition==1, PCtoSort));  % For sorting by PCA
imagesc(cInterval, [], zscore(tagPETHtoPlot(sortKey,:)')', [-3 3]);  ylabel('TAGGED D1 MSNs  (#)');
xlabel('Time from Start (Seconds)'); 
set(gca, 'xtick', [0 6 18], 'ytick', [1 size(tagPETHtoPlot,1)]);
colorbar; colormap('jet');

% Plots averages
subplot(2,2,3); hold on; window = 10; 
tagPETHtoPlot = ctagPETH(tagCondition==2,:); 
means  = (mean(zscore(tagPETHtoPlot')')); %D2s
errors  = (std(zscore(tagPETHtoPlot')'))./sqrt(size(tagPETHtoPlot,1)); %D2s
plot(cInterval, means, 'r', 'LineWidth', 5);  plot(cInterval, means+errors, 'r'); plot(cInterval,means-errors, 'r'); axis tight; 
set(gca, 'xtick', [0 6 18]); colorbar; 
colorbar; % to preserve alignment

subplot(2,2,4); hold on; window = 10; 
tagPETHtoPlot = ctagPETH(tagCondition==1,:); % D1s
means  = (mean(zscore(tagPETHtoPlot')')); %D2s
errors  = (std(zscore(tagPETHtoPlot')'))./sqrt(size(tagPETHtoPlot,1)); %D2s
plot(cInterval, means, 'b', 'LineWidth', 5);  plot(cInterval, means+errors, 'b'); plot(cInterval,means-errors, 'b'); axis tight; 
set(gca, 'xtick', [0 6 18]); colorbar; 
colorbar; % to preserve alignment



cr = cr(valid_idx);
fprintf('\n Waveform correlation between laser on / laser off trials:   %0.2f (%0.2f-%0.2f)\n\n', median(cr), quantile(cr,0.25),quantile(cr,0.75))



