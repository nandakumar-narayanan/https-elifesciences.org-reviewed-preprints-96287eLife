function SpikeData = calcSpikeLatency(dataSt)

%written 5-2-22 to look at evt 32 in dataSt rather than the pl2 event.
addpath('E:\D1D2_Tagging_Ephys\Re__matlab_optotagging\Matlab Offline Files SDK')
%
PathName = './NeuralData/Plexon_NN/'
FileName = dataSt.pl2_name;

file = fullfile(PathName, FileName);

pl2File = file;
pl2Info = PL2GetFileIndex(pl2File);

numSpikeChannels = numel(pl2Info.SpikeChannels);
sortedNeuronNum =[];
for i=1:numSpikeChannels
    tscounts = pl2Info.SpikeChannels{i}.UnitCounts;
    sortedNeuronNum(i) = sum(tscounts(2:end)>0);
end

% Load spike time stamps
sortedNumCh = find(sortedNeuronNum>0);
neuronNum = 1;
unitCharacters={'a','b','c','d','e','f'};
for channelIdx = 1:length(sortedNumCh)
    for neuronIdx = 1:sortedNeuronNum(sortedNumCh(channelIdx))
        [neurons.num_timestamps, neurons.timestamps] = plx_ts(fullfile(PathName,FileName), sortedNumCh(channelIdx), neuronIdx);
        neuronDB{neuronNum} = neurons;
        neuronFileNames{neuronNum} = pl2File;
        neuronName{neuronNum} = sprintf('%s%s',pl2Info.SpikeChannels{sortedNumCh(channelIdx)}.Name,unitCharacters{neuronIdx});
        neuronCh{neuronNum} = pl2Info.SpikeChannels{sortedNumCh(channelIdx)}.Name;
        neuronUnit{neuronNum} = neuronIdx;
        neuronNum = neuronNum+1;
    end
end


eventNum =[];
numEvnetChannels = numel(pl2Info.EventChannels);
for i=1:numEvnetChannels
    % source 9 is Digital IN
    if pl2Info.EventChannels{i}.Source == 9 && pl2Info.EventChannels{i}.NumEvents > 0
        eventNum(i) = pl2Info.EventChannels{i}.Channel;
        %pl2Info.EventChannels{i}.NumEvents
    end
end



evtIdx = find(~(eventNum == 0));
for i=1:length(evtIdx)
    events.(sprintf('evt%d',eventNum(evtIdx(i)))) = PL2EventTs(pl2File, evtIdx(i));
end


[~,name,~] = fileparts(FileName);
name = [name '-01'];
mkdir(name)

% if it's a kwik it needs to be event 13, else it needs to be event 10
% eventTS = oephys_dir_list(2).openEphysEvents.evt32.ts

if dataSt.format == "bin"
    eventTS = events.evt10.Ts;
    eventTS = dataSt.OEphysEvents.evt32.ts
elseif dataSt.format == "kwik"
    eventTS = events.evt13.Ts;
end
% %
% % if contains(pl2File, 'dE')
% %     eventTS = eventTS - 20;
% % end

%
% eventTS = events.evt10.Ts;% - double(events.startTime(1))/30000;
% % eventTS = events.evt5.Ts;% - double(events.startTime(1))/30000;

intervalUP = 0.025; intervalDOWN = -0.025;
% intervalUP = 0.5; intervalDOWN = -0.5;
samplingRate = 30000;
sampleBin = 1/samplingRate;
wnarray = 0.001:0.001:0.02;
wnbinsbase = intervalDOWN:sampleBin:intervalUP;
baseRangeCnt = abs(intervalDOWN/sampleBin);

% for plotting
cmap = [0 0 1;1 0 0];
binSize = 0.0005; %0.5ms
histBins = intervalDOWN:binSize:intervalUP;


%neuronDB-like struct to incorporate the stuff
neuronData = struct()
neuronDataCounter = 1
%latency most important for now
for i_neuron = 1:length(neuronDB)
    
    neurons.timestamps = neuronDB{i_neuron}.timestamps;
    if length(neurons.timestamps) > 4000
    
    edges = sort([eventTS+intervalDOWN;eventTS+intervalUP]);
    [~,~,bin] = histcounts(neurons.timestamps,edges);
    Lia = ismember(bin,1:2:length(edges)); % select within trial(odd number index)
    
    trialTs = neurons.timestamps(Lia);
    trialbin = floor(bin(Lia)/2)+1;
    trialCol = sum(triu(bsxfun(@eq,trialbin,trialbin.'))); % count index number cummulative
    cellTimeTrial = accumarray([trialbin trialCol'],trialTs,[],[],NaN); % make matrix with TS based on raw Idx and colmn Idx
    trialTsIdx = find(Lia); % find index numbers that are within trials
    cellTimeIdxTrial = accumarray([trialbin trialCol'],trialTsIdx,[],[],NaN); % make matrix with TS Idx based on raw Idx and colmn Idx
    cellTimeTrial = cellTimeTrial - repmat(eventTS(1:length(cellTimeTrial)),[1 max(trialCol)]); % raw time stamps to trial time stamps
    preWfIdx = cellTimeIdxTrial(cellTimeTrial<0); postWfIdx = cellTimeIdxTrial(cellTimeTrial>0);
    cellTimeTrial = [cellTimeTrial; nan(length(eventTS) - length(cellTimeTrial),size(cellTimeTrial,2))]; %add NAN raw to match number of trial #ok<AGROW>
    N = zeros(length(cellTimeTrial),length(wnbinsbase)-1);
    for trialNum = 1:length(cellTimeTrial)
        N(trialNum,:) = histcounts(cellTimeTrial(trialNum,:),wnbinsbase);
    end
    
    NBase = N(:,1:baseRangeCnt);
    NStim = N(:,baseRangeCnt+1:end-1);
    for numwin=1:length(wnarray)
        [p, ~] = saltP(NBase,NStim,sampleBin,wnarray(numwin));
        pwn(numwin) = p;
    end
    % find first spike within trial
    [h,cmin2] = max(cellTimeTrial>0,[],2);
    indx = cmin2.*h;
    indx2 = find(indx>0);
    rasterStX = cellTimeTrial(sub2ind(size(cellTimeTrial),indx2,indx(indx2)));
    rasterStY = indx2;
    % for the raster plotting
    [rasterY,col] = find(~isnan(cellTimeTrial));
    cellTime = cellTimeTrial(sub2ind(size(cellTimeTrial),rasterY,col));
    cellTime = cellTimeTrial(~isnan(cellTimeTrial));
    rasterX = cellTime;
    nelements = histc(cellTime,histBins);
    area = sum(nelements) * (max(histBins) - min(histBins))/binSize / (length(nelements)-1);
    [f,xi] = ksdensity(cellTime,'width',binSize );
    %if ~isempty(cellTime) && (length(cellTime)>500)
    
    % spike waveform check
    wave = internalPL2Waves(neuronFileNames{i_neuron},neuronCh{i_neuron} , neuronUnit{i_neuron});
    % find waveforms index per trials, sperate by pre-/post- stimulus
    
    
    wf1 = mean(wave.Waves(preWfIdx,:),1);
    wf2 = mean(wave.Waves(postWfIdx,:),1);
    cr = xcorr(wf1,wf2,0,'coeff');
    
    % waveform analysis
    
    waveformMean = mean(wave.Waves,1);
    upsFactor = 6;
    upsWf = interp(waveformMean,upsFactor); % up sampling to 40*4 KHz
    [peakAmp,peakTime] = min(upsWf);
    [hmaxAmp1, hmaxTime1] = min(abs(upsWf(1:peakTime) - peakAmp/2));
    [hmaxAmp2, hmaxTime2] = min(abs(upsWf(peakTime:end) - peakAmp/2));
    hmaxTime2 = hmaxTime2 + peakTime - 1;
    halfPeakWidthMs = abs(hmaxTime2 - hmaxTime1)/40000/upsFactor*1000;
    [troughAmp, troughTimeIdx] = max(upsWf(peakTime:end));
    troughTime = troughTimeIdx + peakTime - 1;
    peak2troughDur = (troughTime - peakTime) /40000/upsFactor*1000;
    peak2troughRatio = abs(peakAmp / troughAmp);
    
    neuronWfAn(i_neuron,1) = halfPeakWidthMs;
    neuronWfAn(i_neuron,2) = peak2troughDur;
    neuronWfAn(i_neuron,3) = peak2troughRatio;
    
    figure
    subplot(6,1,1:2)
    hold on
    plot(rasterX*1000, rasterY,'.');
    plot(rasterStX*1000, rasterStY,'.r');
    hold off
    xlimp = get(gca,'xlim');
    
    axis off
    axis tight;
    xlimp = get(gca,'xlim');
    title(neuronName{i_neuron});
    % axis off
    
    %set(gca, 'xtick', [], 'ytick', []);
    
    
    subplot(6,1,3:4)
    hold on
    bar(histBins*1000,nelements(1:length(histBins))/size(cellTimeTrial,1)/binSize,'hist');
    histPatch = findobj(gca,'Type','patch'); set(histPatch,'FaceColor','g','EdgeColor','w');
    plot(xi*1000,f*area/length(eventTS),'LineWidth',2)
    plot([0 0],get(gca,'ylim'),'r')
    ylimp = get(gca,'ylim');
    hold off
    axis tight;
    xlim(xlimp);
    xlabel('Time (mSec)'); ylabel('Firing Rate (Hz)')
    
    
    hs3 = subplot(6,1,5:6);
    hold on
    wfx = -9:9/length(wf1):0;
    wfy = ([wf1;wf2;waveformMean]'-min(wf1))*ylimp(2)/range(wf1)*0.8;
    plot(wfx(1:end-1),wfy)
    text(xlimp(1)*0.9, ylimp(2)*0.8, sprintf('CR %0.4f ms',cr))
    text(xlimp(1)*0.9, ylimp(2)*0.6, sprintf('HMPW %2.2f ms',halfPeakWidthMs))
    text(xlimp(1)*0.9, ylimp(2)*0.4, sprintf('P2TD %2.2f ms',peak2troughDur))
    
    
    plot(1:1:20,pwn*ylimp(2)*0.8+ylimp(2)*0.2)
    latency = find(pwn<0.01,1);
    if ~isempty(latency)
        text(latency, ylimp(2)*0.8,['\uparrow' sprintf(' Latency %2.1f ms\n p=%1.4f',latency,pwn(find(pwn<0.01,1)))])
    end
    hold off
    %axis tight;
    axis off
    xlim(xlimp);
    
    p = get(hs3 , 'pos');
    p(2) = p(2) - 0.1;
    set(hs3, 'pos', p);
    
    figStimFileName = fullfile(name,sprintf('%s-stimAligned',neuronName{i_neuron}));
    saveas(gcf,figStimFileName,'png');
    %pause
    close(gcf)
    
    neuronData(neuronDataCounter).latency = latency;
    neuronData(neuronDataCounter).cr = cr;
    neuronDataCounter = neuronDataCounter + 1
    
    else 
    
     neuronData(neuronDataCounter).latency = 99;
    neuronData(neuronDataCounter).cr = .99;
    neuronDataCounter = neuronDataCounter + 1
    
    end
end

SpikeData = neuronData