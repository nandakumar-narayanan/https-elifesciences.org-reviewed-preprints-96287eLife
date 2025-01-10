

clear; 

load neurons_tagging.mat 

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

valid_idx = [neurons.avgFR]>0.5 & [neurons.avgFR] <20 & [neurons.celltype] < 3;
neurons =neurons(valid_idx);

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
TimeInterval = find(interval>[trialStart-binSize]&interval<[trialEnd+binSize]);

for i_loop = 1:1000;
    
    if mod(i_loop,50)==0; fprintf('%d', i_loop); end

    clear randPETH;
    for i_neuron = 1:length(neurons) % 1
            numFR = length(neurons(i_neuron).spikeTS);
            randTS = sort(rand(1,numFR) * max(neurons(i_neuron).spikeTS)); 
            randSpike = peSpike(randTS, neurons(i_neuron).events.switchTrialInit, interval);
            randPETH(i_neuron, :) = gksmooth(randSpike, interval, bw);  % Smooth for PCA
    
    end
    
    crandPETH = randPETH(:, TimeInterval); 
    zrandPETH = zscore(crandPETH')';
    
    
    warning off; [rCOEFF, rSCORE, LATENT, TSQUARED, r_loop] = pca(zrandPETH); warning on;    
    rCOEFF_array(i_loop, :) = rCOEFF(:,1);
    rEXPLAINED_array(i_loop, :) = r_loop; 
    
end; 


rEXPLAINED = mean(rEXPLAINED_array)';

% figure(41); clf; 
% plot(rEXPLAINED./sum(rEXPLAINED)*100, 'ko-')
% box off; xlim([0.5 10]); set(pH, 'MarkerFaceColor', [0 0 0], 'MarkerSize', 10); 
% fprintf('\n PC1 Variance Explained %1.0f', tagEXPLAINED(1)./sum(tagEXPLAINED)*100);
% fprintf('\n PC2 Variance Explained %1.0f', tagEXPLAINED(2)./sum(tagEXPLAINED)*100);
save tagRAND rEXPLAINED*; 
