load neuronDB_matched.mat; 


intStart =  -11;
trialStart = 0; 
trialEnd = 18; 
intEnd = 22;
binSize = 0.2;  % Holds 0.01 - 1 s bins
bw = 1.0;  % Holds 0.1 - 1.5;
%% Will do a supplement exploring bin size and smoothing bandwidth
interval = intStart:binSize:intEnd;
interval_bins = trialStart:binSize:18;

%params for rasters
nexFiles = struct;
params = struct;
params.Raster = 'on';
params.CI = 'on';
params.Bar = 'off';
params.ZeroLine = 'on';
params.RandPerm = 0;
params.Color = [ 0 0 0; 1 0 0; 0 0 1; 1 0 1];  % don't need a colormap command



for i = 1:99; 
j = i+99; 
k = i+99+99; 
name = [neurons(i).fn(1:4) neurons(i).name]; 

                    periEventSpike{1} = peSpike(neurons(i).spikeTS, neurons(i).events.switchTrialInit, interval);
                    switchTime = neurons(i).events.switchDeparts'-neurons(i).events.switchTrialInit;
                    [~,sortKey] = sort(switchTime); periEventSpike{1} = periEventSpike{1}(sortKey,:); 

                    periEventSpike{2} = peSpike(neurons(j).spikeTS, neurons(j).events.switchTrialInit, interval);
                    switchTime = neurons(j).events.switchDeparts'-neurons(j).events.switchTrialInit;
                    [~,sortKey] = sort(switchTime); periEventSpike{2} = periEventSpike{2}(sortKey,:); 

                    periEventSpike{3} = peSpike(neurons(k).spikeTS, neurons(k).events.switchTrialInit, interval);
                    switchTime = neurons(k).events.switchDeparts'-neurons(k).events.switchTrialInit;
                    [~,sortKey] = sort(switchTime); periEventSpike{3} = periEventSpike{3}(sortKey,:); 
                    
                    figure(34); clf; set(gcf, 'Color', 'White');
                    periEventPlot(periEventSpike, interval, params);
                    title(name); set(gca, 'xtick', [0 6 18]); xlabel('Time from Trial Start (Seconds'); 
                    set(0, 'DefaultFigureRenderer', 'painters');
                    saveas(gcf, ['./PharmRasters/' name '.svg']);
        

 
        

end
