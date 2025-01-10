function  [periEventSpikeCondi, epochstat, celldataCondi] = peSpike(spikeTS, eventTS, interval, testEpoch)
% peri-event spike timestamps alignment
% [periEventSpike, celldata] = peSpike(spikeTS, eventTS, interval, testInterval)

% Inputs
%   - spikeTS: raw spike timestamps in seconds, required
%   - eventTS: raw event timestamps in seconds, required
%   - interval: time vector of peri-event time, required (eg: interval = -1:0.1:+1)
%   - testInterval: time vector of peri-event time, optional (eg: interval = -1:0.2:+1)
%     testInterval = [-0.5 0 0.1];  pre-stimulus epoch, stimulus epoch, test bin size
% Outputs
%   - periEventSpike: raw spike timestamps in seconds, required
%       matrix - trial x spike timestamps (aligned to event timestamps)
%   - celldata: raw event timestamps in seconds, optional
%       sparse matrix - trial x binned spike timestamps (aligned to event timestamps)


intervalDOWN = interval(1);
intervalUP = interval(end); 
maxSpike = round((intervalUP-intervalDOWN)*100); % estimate max firing rate < 100Hz

if ~iscell(eventTS)
    eventTS = {eventTS};
end
numCond = length(eventTS);
periEventSpikeCondi = cell(numCond,1);
for condIdx = 1:numCond
    eventTSCond = eventTS{condIdx};
    periEventSpike = NaN(length(eventTSCond),maxSpike); 
    for numTrial = 1:length(eventTSCond)
        trialTs = eventTSCond(numTrial);
        periEvtIdx = spikeTS > (trialTs+intervalDOWN) & spikeTS <= trialTs+intervalUP;
        cellAligned = spikeTS(periEvtIdx)- trialTs;
        periEventSpike(numTrial,1:length(cellAligned)) = cellAligned;
    end
    if size(periEventSpike,2) > maxSpike % if any > 100Hz 
        extraco = [false(size(periEventSpike,1),maxSpike) periEventSpike(:,maxSpike+1:end) == 0];
        periEventSpike(extraco) = NaN;
    end
    lastDataIdx = find(sum(~isnan(periEventSpike)), 1, 'last');
    periEventSpike = periEventSpike(:,1:lastDataIdx);
    periEventSpikeCondi{condIdx} = periEventSpike;
end

if nargout >= 2
    if nargin < 4, error('test ephoch require'); end
    epochstat = cell(numCond,1);
    celldataCondi = cell(numCond,1);
    testBin = testEpoch(3);
    % make sure testEpoch(1) is center +- bin/2
    testEpochLowEdge = testEpoch(1)-testBin/2;
    testEpochHighEdge = testEpoch(1)+testBin/2;
    intervalLowM = fix((interval(1) - testEpochLowEdge)/testBin);
    intervalHighM = fix((interval(end) - testEpochHighEdge)/testBin);
    testInterval = testBin*intervalLowM+testEpochLowEdge:testBin:testBin*intervalHighM+testEpochHighEdge;

    intervalCenter = testInterval(1:end-1) + testBin/2;
    [~, preIdx] = min(abs(intervalCenter - testEpoch(2)));
    [~, testIdx] = min(abs(intervalCenter - testEpoch(1)));
    for condIdx = 1:numCond
        periEventSpike = periEventSpikeCondi{condIdx};
        celldata = histc(periEventSpike', testInterval,1);
        [~, p] = ttest(celldata(preIdx,:), celldata(testIdx,:));
%         [~,ksp,ks2stat] = kstest2(celldata(preIdx,:), celldata(testIdx,:));
%         [perm.PValues, perm.TScores, perm.DFs] = mattestmi(celldata(preIdx,:), celldata(testIdx,:),'permute',5000);
        stat.testp = p;
%         stat.kstestp = ksp;
%         stat.ksteststat = ks2stat;
%         stat.perm = perm;
        stat.testInterval = testInterval;
        stat.preRange = [testInterval(preIdx) testInterval(preIdx+1)];
        stat.testRange = [testInterval(testIdx) testInterval(testIdx+1)];
        epochstat{condIdx} = stat;
        if nargout == 3, celldataCondi{condIdx} = sparse(celldata); end
    end
    if numCond == 1
        epochstat = epochstat{1};
        if nargout == 3, celldataCondi = celldataCondi{1}; end
    end
end

if numCond == 1
    periEventSpikeCondi = periEventSpikeCondi{1};
end


