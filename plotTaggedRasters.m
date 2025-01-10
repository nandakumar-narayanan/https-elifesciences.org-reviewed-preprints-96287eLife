load neuronsTagged.mat; 


intStart =  -2;
trialStart = 0; 
trialEnd = 6; 
intEnd = 8;
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


for i = 1:73; 
                     name = [neuronsTagged(i).fn(1:8) neuronsTagged(i).name]; 
                    clear periEvent*
                    periEventSpike{1} = peSpike(neuronsTagged(i).spikeTS, neuronsTagged(i).events.switchTrialInit, interval);
                    switchTime = neuronsTagged(i).events.switchDeparts'-neuronsTagged(i).events.switchTrialInit;
                    [~,sortKey] = sort(switchTime); periEventSpike{1} = periEventSpike{1}(sortKey,:);

                    figure(22); clf; set(gcf, 'Color', 'White');

                    if neuronsTagged(i).celltype==1; params.Color = [0 0 1]; else; params.Color = [1 0 0]; end; 
 
                    periEventPlot(periEventSpike, interval, params);
                    title(name); set(gca, 'xtick', [0 6 18]); xlabel('Time from Trial Start (Seconds'); 
                    set(0, 'DefaultFigureRenderer', 'painters');
                    saveas(gcf, ['./TaggedRasters/' name '_A.Short.svg']);
        
end
