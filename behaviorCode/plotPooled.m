function PooledAnSt = plotPooled(TrialAnSt, savefolder_name, no_plot)
if nargin < 3
    no_plot = false;
end
%% Plots each session for each animal (across days)
animalNames = fieldnames(TrialAnSt);

centers = 0:0.1:18;
cdf_switch_depart = NaN(numel(centers)-1,length(animalNames));
pdf_switch_depart = NaN(numel(centers)-1,length(animalNames));  
for i_animal = 1:length(animalNames)
    animal = char(animalNames(i_animal));
    numOfDay = sum(~cellfun('isempty', {TrialAnSt.(animal)}));
    TotSwitchTrDe = []; TotSwitchTrAr = [];
    trialWithDataIdx = find(~cellfun('isempty', {TrialAnSt.(animal)}));
    for i_day = 1:numOfDay
        tmp_struct = TrialAnSt(trialWithDataIdx(i_day)).(animal);
        % Pulls out the long trials that should contain a switch
        SwitchTrDe = [tmp_struct.SwitchDepart]; SwitchTrAr = [tmp_struct.SwitchArrival];
        TotSwitchTrDe = [TotSwitchTrDe, SwitchTrDe]; TotSwitchTrAr = [TotSwitchTrAr, SwitchTrAr];
    end
    % Compiles all departure and arrival across animals and days
    [N1,~] = histcounts(TotSwitchTrDe,centers,'Normalization','cdf');
    [f1,~] = ksdensity(TotSwitchTrDe,centers(2:end),'Bandwidth',0.6);
    cdf_switch_depart(:,i_animal) = N1;
    pdf_switch_depart(:,i_animal) = f1;

    % Calculates and builds summary stats table for each animal across days
    [Avg, STD, N, SEM, Median, Q1, Q2, IntQRange, Limen, WR] = getStat(TotSwitchTrDe);

    tSummaryStatsTable = array2table([Avg, STD, N, SEM, Median, Q1, Q2, IntQRange, Limen, WR]);
    tSummaryStatsTable.Properties.VariableNames = {'tAvg', 'tSTD', 'tN', 'tSEM', 'tMedian', 'tQ1', 'tQ2', 'tIntQRange', 'tLimen', 'tWR'};     
    PooledAnSt.(animal).Stats = tSummaryStatsTable;
    PooledAnSt.(animal).cdf_switch_depart = cdf_switch_depart(:,i_animal);
    PooledAnSt.(animal).pdf_switch_depart = pdf_switch_depart(:,i_animal);
    PooledAnSt.(animal).TotSwitchTrDe = TotSwitchTrDe';
    PooledAnSt.(animal).TotSwitchTrAr = TotSwitchTrAr';
    
    if no_plot
        continue
    end        
    
    % Plots the total CDF for each animal across days
    txt = sprintf(' tN = %d\n tMedian = %0.2f\n tAvg = %0.2f  +/- %0.2f\n tInt. Quart. Range = %0.2f',N, Median, Avg, SEM, IntQRange);
    figure; % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0.04 1, 0.96]);
    hold on;
    ecdf(TotSwitchTrDe);
    ecdf(TotSwitchTrAr); 
    text(1, 0.9, txt);
    hold off;
    grid on;
    title(animal);
    yticks(0:0.25:1); ylabel({'Cumulative Probability of a Switch'});
    xlim([0 20]); xticks([0 6 12 18]); xlabel({'Trial Duration (s)'});
    set(gca, 'Color', 'none', 'box','off'); 
    tCDF_FileName = sprintf('tCDF_%s.png', animal); tCDF_fullFileName = fullfile(savefolder_name, tCDF_FileName); saveas(gcf, tCDF_fullFileName);
    close(gcf)
end

if no_plot
    return
end     

%% Plots the pooled CDF for all animals across all days
cmap = lines(1);
trial_time = centers(2:end);
figure; % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0.04 1, 0.96]);
hold on
data = cdf_switch_depart;
cdf50 = [];
for i_animal = 1:size(data,2)
    idx1 = find(data(:,i_animal)<0.5,1,'last');
    idx2 = find(data(:,i_animal)>=0.5,1,'first');
    cdf50(i_animal) = interp1(data([idx1 idx2],i_animal),trial_time([idx1 idx2]),0.5);
end
plotband(trial_time,mean(data,2),std(data,[],2)/sqrt(size(data,2)),cmap)
text(0.25, 0.95, sprintf('cdf50: %0.2f +- %0.2f',mean(cdf50),std(cdf50)/sqrt(size(data,2))), 'Color', 'k', 'FontSize', 12);
data_mean = mean(data,2);
idx1 = find(data_mean<0.5,1,'last');
idx2 = find(data_mean>=0.5,1,'first');
cdf50_mean = interp1(data_mean([idx1 idx2]),trial_time([idx1 idx2]),0.5);
plot([0 cdf50_mean cdf50_mean], [0.5 0.5 0],'Color',cmap)

hold off
grid on
xlabel('Trial Duration'); ylabel({'Cumulative Probability of a Switch'});
set(gca, 'Color', 'none', 'box','off'); 
ylim([0 1])

pCDF_FileName = sprintf('pCDF.png', animal); pCDF_fullFileName = fullfile(savefolder_name, pCDF_FileName); saveas(gcf, pCDF_fullFileName);
close(gcf)
%% Plots the average normalized switch rate
figure
hold on
data = pdf_switch_depart;
plotband(centers(2:end),mean(data,2),std(data,[],2)/sqrt(size(data,2)),cmap)
hold off
xlabel('Trial Duration'); ylabel('Probability Distrubtion Function');
set(gca, 'Color', 'none', 'box','off'); 
aANSR_fullFileName = fullfile(savefolder_name, 'ANSR.png');
saveas(gcf, aANSR_fullFileName); hold off;
close(gcf)
%%
subplot_idx = 1;
animalNames = fieldnames(PooledAnSt);
numplot_y = ceil(sqrt(length(animalNames))); %  positions and sizes for subplot
numplot_x = ceil(length(animalNames)/numplot_y);
figure('Units', 'Normalized', 'OuterPosition', [0 0.04 1, 0.96])
for i_animal = 1:length(animalNames)
    animal = char(animalNames(i_animal));
    data = PooledAnSt.(animal).TotSwitchTrDe;
    [h,p,stats] = chi2gof(data);
    pd = fitdist(data,'Normal');
    fwhm = pd.sigma * 2*sqrt(2*log(2));
    x_pdf = [0:0.1:24];
    y = pdf(pd,x_pdf);

    subplot(numplot_y, numplot_x, subplot_idx);
    hold on
    histogram(data,'Normalization','pdf')
    plot(x_pdf,y,'LineWidth',2)
    plot([pd.mu pd.mu],[0 max(y)],'g')
    plot([pd.mu-fwhm/2 pd.mu+fwhm/2],[max(y)/2 max(y)/2],'r')
    text(0.7,0.70,max(y)/2,sprintf(' mean   %0.2f\n median %0.2f\n\n PDF-gof %0.3f\n PDF-FWHM %0.2f\n PDF-Sigma %0.2f\n PDF-IQR   %0.2f\n', ...
        mean(data), median(data), p, fwhm, pd.sigma, iqr(pd)),'Units','normalized')
    hold off
    title(animal)
    xlabel('Trial Duration'); ylabel({'Switch Departure PDF'});
    set(gca, 'Color', 'none', 'box','off');     
    subplot_idx = subplot_idx+1;
end


end

% builds summary stats table
function [tAvg, tSTD, tN, tSEM, tMedian, tQ1, tQ2, tIntQRange, tLimen, tWR] = getStat(x)
    tAvg = mean(x); 
    tSTD = std(x); 
    tN = size(x,2); 
    tSEM = tSTD/(sqrt(tN));
    tMedian = median(x); 
    tQ1 = quantile(x, 0.25); 
    tQ2 = quantile(x, 0.75); 
    tIntQRange = tQ2 - tQ1; tLimen = tIntQRange/2; tWR = tLimen/tMedian;
    % bstats = [tAvg, tSTD, tN, tSEM, tMedian, tQ1, tQ2, tIntQRange, tLimen, tWR];
end

