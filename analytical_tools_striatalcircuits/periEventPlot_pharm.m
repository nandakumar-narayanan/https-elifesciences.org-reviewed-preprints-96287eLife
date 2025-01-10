 function varargout = periEventPlot(periEvtSPK, interval, params, drug)
% plot peri event raster plot and kernel estimated neuronal firing rate
% Inputs
% periEvtSPK
%   single condition - periEvtSPK is matrix (trial x spike timestamps)
%   multiple condition - periEvtSPK is cell with matrix in (trial x spike timestamps)
% interval
%   time vector of peri-event time, required (eg: interval = -1:0.1:+1)
%   should be same interval for periEvtSPK
% variable input arguments
%   - 'Raster', 0 % default 1 'on', trial by trial spike raster
%   - 'Title','title text'
%   - 'Legend', {'a' 'b' 'c'}
%   - 'CI', 1 % default 0 'off', confidence intervals (95% and 5%)
%   - 'Bar', 1 % default 0 'off', histogram in estimated firing rate line plot
%   - 'ZeroLine', 0 % default 1 'on', add vertical line at x=0 in estimated firing rate line plot
%
% Output - optional, figure handle

% default kernel bandwidth and time output limit to remove edge artifact
binsize = interval(2)-interval(1);
kernelWin = binsize*2;
tlimit = [interval(1) + kernelWin*4 interval(end) - kernelWin*4];

% variable input arguments handling
p = getParams(params);

% single condition 
if ~iscell(periEvtSPK)
    periEvtSPK = {periEvtSPK};
end

numCond = length(periEvtSPK);
% randomly select trials for raster plot
% ipermut = cell(numCond,1);
% if isfield(params, 'RandPerm')
%     numTrialSel = params.RandPerm;
%     for condIdx = 1:numCond
%         trialSpkIdx = find(sum(~isnan(periEvtSPK{condIdx}), 2) > 20);
%         numTrialTotal = length(trialSpkIdx);
%         if numTrialSel > 1
%             if numTrialSel > numTrialTotal, numTrialSel = numTrialTotal; end
%             ipermut{condIdx} = sort(trialSpkIdx(randperm(numTrialTotal, numTrialSel)));
%         else
%             ipermut{condIdx} = sort(trialSpkIdx(randperm(numTrialTotal, numTrialTotal*numTrialSel)));
%         end
%     end
% else
%     for condIdx = 1:numCond
%         numTrialTotal = size(periEvtSPK{condIdx},1);
%         ipermut{condIdx} = 1:numTrialTotal;
%     end
% end

%h = figure;
h = gcf;

% Raster plot
if p.addRaster
    numsubplot = 2;
    for condIdx = 1:numCond
        perTrialSpk = periEvtSPK{condIdx}(:,:);  % Removed the permutation stuff
        hold on
        subplot(numsubplot,1,1)
        if condIdx == 1
            periEventRaster(perTrialSpk, tlimit, 'Color', params.Color(drug,:), 'FontSize',4);
%             periEventRaster(perTrialSpk, tlimit, 'Color', p.cmap(condIdx,:), 'FontSize',4);
        elseif condIdx > 1
%             periEventRaster(perTrialSpk, tlimit, 'Color', p.cmap(condIdx,:), 'FontSize', 8, 'LineWidth', 12)
            periEventRaster(perTrialSpk, tlimit, 'Color', params.Color(drug,:), 'FontSize', 8, 'LineWidth', 12)
        end    
    end

    % Trick to remove gaps between raster subplots
%     rasterPadSize = 0.05;
%     subHeight = (0.5 - rasterPadSize*2)/numCond;
%     bottomPostion = 1-(rasterPadSize+subHeight*(1:numCond));
%     if numCond > 1
%         bottomPostion(2:end) = bottomPostion(2:end)-cumsum(ones(numCond-1,1)*rasterPadSize/3)';
%     end
%     hAllAxes = findobj(gcf,'type','axes');
%     for axisIdx = 1:length(hAllAxes)
%         set(hAllAxes(axisIdx),'Position',[hAllAxes(1).Position(1) bottomPostion(axisIdx) hAllAxes(1).Position(3) subHeight]) % [left bottom width height]
%     end
%   
end

hold off

% determine time scale, smaller than 500ms use msec unit
if range(interval) > 0.5
    timeScale = 1;
    timeLabel = 'Time (sec)';
else
    timeScale = 1000;
    timeLabel = 'Time (msec)';
end

% plot estimated firing rate

if p.addRaster
    subplot(numsubplot,1,2)
%     subplot(numsubplot,1,numCond+1:numsubplot)
end

% Trick to legend to be proper.
% for condIdx = 1:numCond
for condIdx = 1
    if ~isempty(periEvtSPK{condIdx})
        data = gksmooth(periEvtSPK{condIdx}, interval, kernelWin);
        subplot(numsubplot,1,2)
        plot(interval.*timeScale, data,'LineWidth',2, 'Color', params.Color(drug,:))
    end      
end
%legend(p.legendText, 'Box','off');
title(p.titleText);
for condIdx = 1
    if ~isempty(periEvtSPK{condIdx})
        if p.addBar
            periEvtSPKCond = periEvtSPK{condIdx};
            cellTime = periEvtSPKCond(~isnan(periEvtSPKCond));
            nelements = histc(cellTime, interval);
            barY = nelements(1:length(interval))./size(periEvtSPKCond, 1)./ binsize;
%             bar(interval.*timeScale, barY, 'FaceColor', p.cmap(condIdx,:), 'FaceAlpha',0.3, 'EdgeColor','none');
            bar(interval.*timeScale, barY, 'FaceColor', params.Color(drug,:), 'FaceAlpha',0.3, 'EdgeColor','none');
        end
        if ~p.addCI
            data = gksmooth(periEvtSPK{condIdx}, interval, kernelWin);
%             plot(interval.*timeScale, data,'LineWidth',2, 'Color', p.cmap(condIdx,:))
            plot(interval.*timeScale, data,'LineWidth',2, 'Color', params.Color(drug,:))
        else
            [data, ci] = gksmooth(periEvtSPK{1}, interval, kernelWin);
%             plotband(interval.*timeScale, data, ci, p.cmap(1,:))
            plotband(interval.*timeScale, data, ci, params.Color(drug,:))
        end
    end
end
xlim(tlimit.*timeScale);
set(gca,'box','off');
set(gca,'color','none');
hold on; 
if p.addZeroline
    plot([0 0],get(gca,'ylim'),'g')  %trial init
    plot([18 18],get(gca,'ylim'),'m')%trial end
end

if p.epochstat
    if ~iscell(params.epochstat)
        params.epochstat = {params.epochstat};
    end
    epochstat = params.epochstat{1};
    xtest = [epochstat.testRange(1) epochstat.testRange(1) epochstat.testRange(2) epochstat.testRange(2)];
    xpre = [epochstat.preRange(1) epochstat.preRange(1) epochstat.preRange(2) epochstat.preRange(2)];
    yrange = get(gca,'ylim');
    fill(xtest,[fliplr(yrange) yrange] ,'r', 'LineStyle', 'none')
    fill(xpre,[fliplr(yrange) yrange] ,'b', 'LineStyle', 'none')
    alpha(0.05)
    %text(max(xpre),yrange(2)-range(yrange)*0.05,sprintf('ttest: %0.3f',epochstat.testp'));
    for condIdx = 1:numCond
        text(max(xpre),yrange(2)-range(yrange)*0.05*condIdx,sprintf('permp: %0.3f',params.epochstat{condIdx}.perm.PValues'));
    end
end
ax = gca;
grayness = 0.3;
ax.XColor = [grayness grayness grayness];
ax.YColor = ax.XColor;
ylabel('Firing Rate (Hz)')
xlabel(timeLabel);
hold off
        
if nargout==1
    varargout{1} = h;
end

end

function varargout = periEventRaster(periEventSpike, tlimit, varargin)

% plot spike raster using text '|', scale better than line plot
alphaThresh = 40000;
dataIdx = ~isnan(periEventSpike);
[row,~] = find(dataIdx);
rasterSpike = periEventSpike(dataIdx);
limitIdx = rasterSpike > tlimit(1) & rasterSpike < tlimit(2);
rasterSpike = rasterSpike(limitIdx);
row = row(limitIdx);
plot(rasterSpike, row, 'LineStyle', 'none', 'Marker', 'none');
t = text(rasterSpike, row, '|', 'HorizontalAlignment','center', varargin{:});
if length(rasterSpike) > alphaThresh
    alphaNum = alphaThresh/length(rasterSpike);
    alpha(t,alphaNum);
end
%axis tight
xlim(tlimit);
axis off
        
if nargout==1
    varargout{1} = t;
end

end

function plotband(x, means, variance, color)

px = [x fliplr(x)];
if size(variance,1) == 1 % single variance like std or sem
    patch(px, [means+variance' fliplr(means-variance')], color, 'EdgeAlpha', 0.1, 'FaceAlpha', 0.3, 'linewidth',0.5);
    hold on
    plot(x, means, 'LineWidth', 4, 'Color', 'k')
%     plot(x, means, 'LineWidth', 4, 'Color', color)
elseif size(variance,1) == 2 % dual variance like confidence intervals
    patch(px, [variance(1,:) fliplr(variance(2,:))], color, 'EdgeAlpha', 0.1, 'FaceAlpha', 0.3, 'linewidth',0.5);
    hold on
    plot(x, means, 'LineWidth', 2, 'Color',color)
end

end

function p = getParams(params)

if isfield(params, 'Raster')
    p.addRaster = strcmpi(params.Raster, 'on');
else
    p.addRaster = 1;
end

if isfield(params, 'CI')
    p.addCI = strcmpi(params.CI, 'on');
else
    p.addCI = 0;
end

if isfield(params, 'Bar')
    p.addBar = strcmpi(params.Bar, 'on');
else
    p.addBar = 0;
end

if isfield(params, 'ZeroLine')
    p.addZeroline = strcmpi(params.ZeroLine, 'on');
else
    p.addZeroline = 1;
end

if isfield(params, 'Title')
    p.titleText = params.Title;
else
    p.titleText = '';
end

if isfield(params, 'Legend')
    p.legendText = params.Legend;
else
    p.legendText = {};
end

if isfield(params, 'Color')
    p.cmap = params.Color;
else
    p.cmap = get(groot,'DefaultAxesColorOrder');
end

if isfield(params, 'epochstat')
    p.epochstat = 1;
else
    p.epochstat = 0;
end

end