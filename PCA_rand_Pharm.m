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

figure(41); clf; 
plot(rEXPLAINED./sum(rEXPLAINED)*100, 'ko-')
box off; xlim([0.5 10]); set(pH, 'MarkerFaceColor', [0 0 0], 'MarkerSize', 10); 
fprintf('\n PC1 Variance Explained %1.0f', EXPLAINED(1)./sum(EXPLAINED)*100);
fprintf('\n PC2 Variance Explained %1.0f', EXPLAINED(2)./sum(EXPLAINED)*100);

length(find(rEXPLAINED_array(:,1)>EXPLAINED(1)))
length(find(rEXPLAINED_array(:,2)>EXPLAINED(2)))

save pharmRAND rEXPLAINED*; 

