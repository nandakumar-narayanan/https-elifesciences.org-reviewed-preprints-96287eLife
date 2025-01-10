
% Interval parameters
intStart =  -4;
trialStart = 0; 
trialEnd = 6; 
intEnd = 22;
binSize = 0.2;  % Holds 0.01 - 1 s bins
bw = 1.0;  % Holds 0.1 - 1.5;
interval = intStart:binSize:intEnd;


 clear *PETH*
counter=1; clear tagPETH stagPETH tagCondition;
for i_neuron = 1:length(neurons) % 1
 
        periEventSpike = peSpike(neurons(i_neuron).spikeTS, neurons(i_neuron).events.switchTrialInit, interval);
        periEventMotor = peSpike(neurons(i_neuron).events.shortPoke, neurons(i_neuron).events.switchTrialInit, interval);  % few switches in the first 6 seconds, and sometimes no nosepokes. 


        spikes = []; times = []; motor = []; 
        for j_trial = 1:size(periEventSpike,1)
            switchTime = neurons(i_neuron).events.switchDeparts(j_trial)-neurons(i_neuron).events.switchTrialInit(j_trial);
            interval_bins = [0:binSize:switchTime];
            spikes = [spikes histc(periEventSpike(j_trial,:),interval_bins)]; % time x trial
            motor = [motor histc(periEventMotor(j_trial,:),interval_bins)]; % time x trial
            times =  [times interval_bins]; 
        end


        T_neuron = table(spikes', motor', times','VariableNames', {'FiringRate', 'Motor', 'Time'});

        lmNeuron = fitglme(T_neuron, 'FiringRate~Time');  % In the case of poor fit, where there are few nosepokes, we revert to time.

        lmNeuronAnova = anova(lmNeuron);
        pTagTime(i_neuron) = [lmNeuronAnova{2,5}];
        a = fixedEffects(lmNeuron); 
        slopeSwitch(i_neuron) = a(2);
    
end

tagCondition = [neurons.celltype];
tagCondition = tagCondition(valid_idx);
slopeSwitch = slopeSwitch(valid_idx);


figure(39)

cla; hold on; clear data; 
data{1} = slopeSwitch(tagCondition==2)./binSize; %Switched so that the plot appears first
data{2} = slopeSwitch(tagCondition==1)./binSize;
jitterPlot(data, [1 0 0; 0 0 1]);
ylabel('Slope'); 
set(gca, 'xtick', [1:2], 'xticklabel', {'D2'; 'D1'}); 

fprintf('\n\n D2 Slope %0.2f (%0.2f-%0.2f)', median(data{1}), quantile(data{1},0.25),quantile(data{1},0.75))
fprintf('\n D1 Slope %0.2f (%0.2f-%0.2f)', median(data{2}), quantile(data{2},0.25),quantile(data{2},0.75))
fprintf('\n Slope: D1 vs D2 p = %0.3f, cohen d %2.1f', ranksum(data{1}, data{2}),cohend(data{1}, data{2}) )
csvwrite('./stats/zTaggedSlope_D2.csv', data{1}');
csvwrite('./stats/zTaggedSlope_D1.csv', data{2}');


% Have to run Figure 3 first; 
t = table(slope'./binSize, slopeSwitch'./binSize, tagSCORE(:,1), tagCondition',  animalArray', direction', Sex', FR', 'VariableNames', {'Slope', 'slopeSwitch', 'PC1', 'Type', 'Animal', 'direction', 'Sex', 'FR'});
writetable(t, 'zSlopeForLMM.csv');




corrcoef(t.Slope, t.PC1);







