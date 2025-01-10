% combine loadPharmNeuronDB / loadRawData_tagging

%Create_dataSt-------------------------------------------------------------
addpath(genpath('./analysis_tools'))
addpath('./behaviorCode')
addpath('./analytical_tools_striatalcircuits')

process_raw_data('pharm');
process_raw_data('tagging');


function process_raw_data(data_type)

if strcmp(data_type,'pharm')
    load_pharm = true;
    load_tagging = false;
elseif strcmp(data_type,'tagging')
    load_pharm = false;
    load_tagging = true;
end

     
if load_pharm
    load('OEphysEvents_Pharma12.mat') %Extracted event data from OEphysfiles
    behavior_data = './MPCpharm'; %the medPC data
    pl2_dirs = './NeuralData/Plexon_NN/';   %MergeSort is for matched rasters
    %behaviour analysis from MPC code
    typecell = {
        'ER5' 1,1,NaN;...
        'ER4' 1,1,NaN;...
        'ER3' 1,1,NaN;...
        'EP8' 1,1,NaN;...
        'iER1' 1,1,NaN;...
        'iER2' 1,1,NaN;...
        'iER4' 1,1,NaN;...
        'iEL3' 0,1,NaN;...
        'iEL5' 0,1,NaN;...
        'dER1' 1,1,NaN;...
        'dER2' 1,1,NaN;...
        'dER6' 1,1,NaN;...
        'dER9' 1,1,NaN;...
        'dEL3' 0,1,NaN;...
        'dEL7' 0,1,NaN;...
        'dEL8' 0,1,NaN;...
        'RCL3' 0,1,NaN;...
        'RCL4' 0,1,NaN;...
        };
    idMatch = {'iER1' 'iER2' 'ER3' 'ER4' 'ER5' 'EP8' 'iER4' 'iEL3' 'iEL4' 'dER1' 'dER2' 'dEL3' 'iEL5' 'dEL7' 'dER6' 'RCL3' 'dER9' 'RCL4' 'dEL8'};
    
    MSN = {'Ben_Switch_18L6R_SHORTviITI' 'Ben_Switch_6L18R_SHORTviITI' 'Switch_6L18R_S_ITI' 'Switch_18L6R_S_ITI' 'Switch_18L6R_SITI_RI_AUSTIN' 'Switch_18L6R_SITI_RI_AUSTIN_OPTO' 'Switch_6L18R_SITI_RI_AUSTIN_OPTO' 'Switch_18L6R_SITI_RI_AUSTIN' 'Switch_6L18R_SITI_RI_AUSTIN'};
    mpcParsed = getDataIntr_pharm(behavior_data, MSN, idMatch); % need to double check on this.

elseif load_tagging
    load ('OEphysEvents_tagging_v4.mat')  %open ephys has the event codes, but these files are too large to host
    behavior_data = './mpcData_tagging'; %the medPC data directories for files
    pl2_dirs = './NeuralDataTagging/Sorted/';

    %behaviour analysis from MPC code
    typecell = {
        'iER1' 1,1,NaN;...
        'iER2' 1,1,NaN;...
        'iEL3' 0,1,NaN;...
        'iER4' 1,1,NaN;...
        'iEL5' 0,1,NaN;...
        'dER1' 1,1,NaN;...
        'dER2' 1,1,NaN;...
        'dEL1' 0,1,NaN;...
        'dEL2' 0,1,NaN;...
        'dEL3' 0,1,NaN;...
        'dEL4' 0,1,NaN;...
        'dEL7' 0,1,NaN;...
        'dEL8' 0,1,NaN;...
        'dER9' 0,1,NaN;...
        };
    MSN = {'Switch_18L6R_SITI_RI_AUSTIN', 'Switch_6L18R_SITI_RI_AUSTIN'};
    %mpcParsed1 = getDataIntr_v1(behavior_data, MSN);
    mpcParsed = getDataIntr_pharm(behavior_data, MSN);
    mpcParsed = mpcParsed(~cellfun(@isempty, regexpi({mpcParsed.Subject}, '(dE|iE)\D+')));
        
end

%%
pl2_dir_list = dir([pl2_dirs '*.pl2']) ;
%creates dataSt, a struct concatenating info from pl2, openephys, and mpc files
dataSt = struct;
for i_file = 1:length(pl2_dir_list)
    fileName = fullfile(pl2_dir_list(i_file).folder, pl2_dir_list(i_file).name);
    [~,name,~] = fileparts(fileName);
    animal_id = regexp(name,'\w{2,3}\d{1}','match','once');
    odate = regexp(name,'\d+\-\d+\-\d+','match','once');
    nex_rec_date = datetime(odate,'InputFormat','yyyy-MM-dd');
    mpc_rec_date = datestr(nex_rec_date,'mm/dd/yy');
    
    dataSt(i_file).animal_id = animal_id;
    dataSt(i_file).pl2 = fileName;
    dataSt(i_file).odate = odate;
    dataSt(i_file).nex_rec_date = nex_rec_date;
    dataSt(i_file).mpc_rec_date = mpc_rec_date;
    dataSt(i_file).pl2_name = strcat(name, '.pl2');
    dataSt(i_file).oephys_name = name(1:regexp(name, 'sort')-5);

    %coding early as 1
    if contains(name, 'SALINE')
        dataSt(i_file).drug = 1;
    elseif contains(name, 'SUL')
        dataSt(i_file).drug = 2;
    elseif contains(name, 'SCH')
        dataSt(i_file).drug = 3;
    else
        dataSt(i_file).drug = NaN; %no drug administered
    end

    
    %may be redundant?
%     if contains(name, 'kwik')
%         dataSt(i_file).format = 'kwik';
%     elseif contains(name, 'bin')
%         dataSt(i_file).format = 'bin';
%     end

    if datetime(odate) < datetime('01/30/2021')
        dataSt(i_file).format = 'kwik';
    else
        dataSt(i_file).format = 'bin';
    end  
    
    if contains(name, 'iE')
        dataSt(i_file).genotype = 'D2-cre';
    elseif contains(name, 'dE')
        dataSt(i_file).genotype = 'D1-cre';
    else
        dataSt(i_file).genotype = 'C57';
    end
end
%%


for i_file = 1:length(dataSt)
    match_idx = find(strcmp(dataSt(i_file).animal_id, {mpcParsed.Subject}) & strcmp(dataSt(i_file).mpc_rec_date, {mpcParsed.StartDate}));
    if numel(match_idx) == 1
        dataSt(i_file).mpc = mpcParsed(match_idx);
        dataSt(i_file).TrialAnSt = getTrialData_Switch(dataSt(i_file).mpc, typecell);
    end
end

for i_file = 1:length(dataSt)
    match_idx = find(strcmp(dataSt(i_file).animal_id, {oephys_dir_list.sid}) & strcmp(dataSt(i_file).odate, {oephys_dir_list.odate}));
    if numel(match_idx) == 1 && dataSt(i_file).format == "bin"
        dataSt(i_file).OEphysEvents = oephys_dir_list(match_idx).openEphysEvents;
    elseif dataSt(i_file).format == "kwik"
        dataSt(i_file).OEphysEvents = get_mpc_bin_event_pl2_v2(dataSt(i_file).pl2);
    end
end

%%

%ttl checker
for i_file = 1:length(dataSt)
    
    animal_idx = dataSt(i_file).animal_id;
    %     for i_step = 1:length({[dataSt(j_file).TrialAnSt.(animal_idx).realTrialStart]})
    
    % % % let step_value be a constant of some type amenable to adding and subtracting with timestamps
    % % % let list1 and list2 be lists of arbitrary length with at least 2 fields; a timestamp of a type you can do math on and data associated with the trial record.
    % % % List2 is longer than list1, but has some number of records (equal, of course, to the number of records in list1) that form a subset
    % % %
    % % % we'll call list1 the "records" and list2 the "candidates". then we will build a third list: the "matches".
    % % %
    records = [dataSt(i_file).TrialAnSt.(animal_idx).realTrialStart];
    candidates =  [dataSt(i_file).OEphysEvents.evt3.ts];
    candidates(find(diff(candidates)<0.0001)+1) = [];
    
    matches = [];
    expected_diff = candidates(1) - records(1);
    
    for i_record = 1:length(records)
        expected_value = records(i_record) + expected_diff;
        for i_candidate = 1:length(candidates)
            if candidates(i_candidate) > expected_value - 1 &  candidates(i_candidate) < expected_value + 1  %can expect 100ms jitter here.
                matches(i_record) = candidates(i_candidate);
            end
        end
    end

    end_records = [dataSt(i_file).TrialAnSt.(animal_idx).trialEnd];
    end_candidates = dataSt(i_file).OEphysEvents.evt17.ts;
    end_candidates(find(diff(candidates)<0.0001)+1) = [];
    
    matches_end = [];
    expected_diff_end = end_candidates(1) - end_records(1);
    
    for i_record = 1:length(end_records)
        expected_value = end_records(i_record) + expected_diff_end;
        for i_candidate = 1:length(end_candidates);
            if end_candidates(i_candidate) > expected_value - 1 &  end_candidates(i_candidate) < expected_value + 1 ; %can expect 100ms jitter here.
                matches_end(i_record) = end_candidates(i_candidate);
            end
        end
    end
    
    if ~length(matches_end) == length(matches)
        fprintf('uh-oh')
    end 
    
    if ~length(matches) == length([dataSt(i_file).TrialAnSt.(animal_idx).realTrialStart])
        fprintf('double uh-oh')
    end 

    dataSt(i_file).OEphysEvents.evt3.ts(:) = [];
    dataSt(i_file).OEphysEvents.evt3.ts = [matches];
    dataSt(i_file).OEphysEvents.evt17.ts(:) = [];
    dataSt(i_file).OEphysEvents.evt17.ts = [matches_end];
    
    dataSt(i_file).pl2_events = dataSt(i_file).OEphysEvents;
end

%%

for i_file = 1:length(dataSt)

    animal_idx = dataSt(i_file).animal_id;
    
    %extract all trial inits
    
    %startTime = [dataSt(i_file).OEphysEvents.evt1.ts(1)];
    trialTypes = [dataSt(i_file).TrialAnSt.(animal_idx).programmedDuration];
    result = [dataSt(i_file).TrialAnSt.(animal_idx).outcome];
    
    % Determine if there's a switch
    
    switches = ~cellfun(@isempty, {dataSt(i_file).TrialAnSt.(animal_idx).SwitchDepart});
    
    shortTrials = [trialTypes == 6000 & (result == 4 | result == 2)];
    longTrials =  [trialTypes == 18000 & (result == 4 | result == 2)];
    shorterrorTrials = [trialTypes == 6000 & (result == 3 | result == 1)];
    longerrorTrials = [trialTypes == 18000 & (result == 3 | result == 1)];
    switchTrials = [switches == 1];
    
    %logical loop to exlcude cases in which MedPC terminated prior to final
    %trial completion. 2 possible scenarios--animal never completed trial,
    %OR medPC protocol timed out (90 minutes) during Trial.
    
    % calculate SwitchDepart
    animals = dataSt(i_file).animal_id;
    trialStart = dataSt(i_file).pl2_events.evt3.ts; % CUE ON
    trialEnd = dataSt(i_file).pl2_events.evt17.ts; % TRIAL END
    rpInLeft = dataSt(i_file).pl2_events.evt7.ts; % LEFT RESPONSE
    rpOutLeft = dataSt(i_file).pl2_events.evt11.ts; % LEFT RELEASE
    rpInRight = dataSt(i_file).pl2_events.evt19.ts; % RIGHT RESPONSE
    rpOutRight = dataSt(i_file).pl2_events.evt15.ts;% RIGHT RELEASE
    rewards = dataSt(i_file).pl2_events.evt9.ts';% REWARD DISPENSE
    
    shortTrials(find(diff(trialStart)<0.0001)+1) = [];
    longTrials(find(diff(trialStart)<0.0001)+1) = [];
    shorterrorTrials(find(diff(trialStart)<0.0001)+1) = [];
    longerrorTrials(find(diff(trialStart)<0.0001)+1) = [];
    switchTrials(find(diff(trialStart)<0.0001)+1) = [];
    trialEnd(find(diff(trialStart)<0.0001)+1) = [];
    trialStart(find(diff(trialStart)<0.0001)+1) = [];
    
    trial = struct;
    trialNum = min(length(trialStart), length(trialEnd)); % reevaluate this

    type = typecell{strcmp(animals,typecell(:,1)),2};

    for j = 1:trialNum
        % trial start/end time, duration, and ITI
        curTS = trialStart(j); % H
        curTE = trialEnd(j);
        trial(j).trialStart = curTS;
        % response time within trial
        trial(j).leftRspTimeTrial = rpInLeft(rpInLeft>curTS & rpInLeft<=curTE) - curTS;
        trial(j).leftRelTimeTrial = rpOutLeft(rpOutLeft>curTS & rpOutLeft<=curTE) - curTS;
        trial(j).rightRspTimeTrial = rpInRight(rpInRight>curTS & rpInRight<=curTE) - curTS;
        trial(j).rightRelTimeTrial = rpOutRight(rpOutRight>curTS & rpOutRight<=curTE) - curTS;
        
        if type == 0    % 0 indicates that a short latency trial is rewarded at the left nose poke
            if ~isempty(trial(j).leftRelTimeTrial) && ~isempty(trial(j).rightRspTimeTrial)
                trial(j).SwitchArrival = min(trial(j).rightRspTimeTrial(trial(j).rightRspTimeTrial > min(trial(j).leftRelTimeTrial)));
                if ~isempty(trial(j).SwitchArrival)
                    trial(j).SwitchDepart = max(trial(j).leftRelTimeTrial(trial(j).leftRelTimeTrial < trial(j).SwitchArrival));
                else
                    trial(j).SwitchArrival=[];
                end
            end
            
        elseif type == 1    % 1 indicates that a short latency trial is rewarded at the right nose poke
            if ~isempty(trial(j).rightRelTimeTrial) && ~isempty(trial(j).leftRspTimeTrial)
                trial(j).SwitchArrival = min(trial(j).leftRspTimeTrial(trial(j).leftRspTimeTrial > min(trial(j).rightRelTimeTrial)));
                if ~isempty(trial(j).SwitchArrival)
                    trial(j).SwitchDepart = max(trial(j).rightRelTimeTrial(trial(j).rightRelTimeTrial < trial(j).SwitchArrival));
                else
                    trial(j).SwitchArrival=[];
                end
            end
        end
        
    end
    
      
    starts = trialStart;
    
    shortStarts = starts(shortTrials);
    longStarts = starts(longTrials);
    %switchStarts = starts(switchTrials & longTrials);
    shortErrorStarts = starts(shorterrorTrials);
    longErrorStarts = starts(longerrorTrials);

    sw_trial = ~cellfun(@(x) isempty(x), {trial.SwitchDepart});

    
    if ~isequal(sw_trial, switchTrials)
        fprintf("uh-oh %d\n", i_file)
    end
    
    switchStarts = starts(switchTrials);
    switchDeparts = [trial(sw_trial).trialStart]' + [trial(sw_trial).SwitchDepart]';
    switchArrivals = [trial(sw_trial).trialStart]' + [trial(sw_trial).SwitchArrival]';
    
    events = struct;
    events.shortTrialInit =  shortStarts;
    events.longTrialInit =  longStarts;
    events.switchTrialInit = switchStarts;
    events.switchDeparts = switchDeparts;
    events.switchArrives = switchArrivals;
    events.Rewards = rewards;
    events.shortErrorTrialInit = shortErrorStarts;
    events.longErrorTrialInit = longErrorStarts;

    events.rightPoke = dataSt(i_file).pl2_events.evt19.ts;
    events.leftPoke = dataSt(i_file).pl2_events.evt7.ts;
    if type == 0
        events.shortPoke = dataSt(i_file).pl2_events.evt19.ts;
       events.longPoke = dataSt(i_file).pl2_events.evt7.ts;

    elseif type == 1
        events.shortPoke = dataSt(i_file).pl2_events.evt7.ts;
        events.longPoke = dataSt(i_file).pl2_events.evt19.ts;

    end
    if dataSt(i_file).format == "bin"
        %eventTS = events.evt10.Ts;
        events.laser_ttl = dataSt(i_file).OEphysEvents.evt32.ts;
    elseif dataSt(i_file).format == "kwik"
        evt13 = PL2EventTs(dataSt(i_file).pl2, 13);
        events.laser_ttl = evt13.Ts;
    end    
   dataSt(i_file).events = events;        
end

%%


neurons = struct;
neuronNum = 1;
for i_file = 1:length(dataSt)
    pl2File = dataSt(i_file).pl2;
    pl2Info = PL2GetFileIndex(pl2File);
    numSpikeChannels = numel(pl2Info.SpikeChannels);
    sortedNeuronNum =[];
    for i=1:numSpikeChannels
        tscounts = pl2Info.SpikeChannels{i}.UnitCounts;
        sortedNeuronNum(i) = sum(tscounts(2:end)>0);
    end

    % Load spike time stamps
    sortedNumCh = find(sortedNeuronNum>0);
    unitCharacters={'a','b','c','d','e','f'};
    for channelIdx = 1:length(sortedNumCh)
        for neuronIdx = 1:sortedNeuronNum(sortedNumCh(channelIdx))
            wave = internalPL2Waves(pl2File,pl2Info.SpikeChannels{sortedNumCh(channelIdx)}.Name , neuronIdx);
            waveformMean = mean(wave.Waves,1);
            neuronName = sprintf('%s%s',pl2Info.SpikeChannels{sortedNumCh(channelIdx)}.Name,unitCharacters{neuronIdx});
            [~, timestamps] = plx_ts(pl2File,pl2Info.SpikeChannels{sortedNumCh(channelIdx)}.Name , neuronIdx);
            neurons(neuronNum).name = neuronName;
            neurons(neuronNum).spikeTS = timestamps;
            neurons(neuronNum).events = dataSt(i_file).events;
            neurons(neuronNum).fn = dataSt(i_file).pl2_name;
            neurons(neuronNum).genotype = dataSt(i_file).genotype;
            neurons(neuronNum).condition = dataSt(i_file).drug;
            neurons(neuronNum).wf.mean = waveformMean;
            neurons(neuronNum).pl2File = pl2File;
            neurons(neuronNum).ch_name = pl2Info.SpikeChannels{sortedNumCh(channelIdx)}.Name;
            neurons(neuronNum).ch_idx = neuronIdx;
            neurons(neuronNum).nameDB = strcat(dataSt(i_file).animal_id,'_',neuronName, '_', dataSt(i_file).odate);
            neurons(neuronNum).numTS = length(neurons(neuronNum).spikeTS); 
            neurons(neuronNum).avgFR = length(neurons(neuronNum).spikeTS)./max(neurons(neuronNum).spikeTS);   
            neuronNum = neuronNum+1;
        end
    end
end

% up-scaling and extract features

upsFactor = 6;
sampleRate = 60000;
timeUnitConv = 1000; % sec to milli second    
for i_neuron = 1:length(neurons)
    % waveform analysis
    waveformMean = neurons(i_neuron).wf.mean ;

    upsWf = interp(waveformMean,upsFactor); % up sampling to 40*4 KHz
    [peakAmp,peakTime] = min(upsWf);
    [hmaxAmp1, hmaxTime1] = min(abs(upsWf(1:peakTime) - peakAmp/2));
    [hmaxAmp2, hmaxTime2] = min(abs(upsWf(peakTime:end) - peakAmp/2));
    hmaxTime2 = hmaxTime2 + peakTime - 1;
    halfPeakWidthMs = abs(hmaxTime2 - hmaxTime1)/sampleRate/upsFactor*timeUnitConv;

    [troughAmp, troughTimeIdx] = max(upsWf(peakTime:end));
    troughTime = troughTimeIdx + peakTime - 1;
    peak2troughDur = (troughTime - peakTime) /sampleRate/upsFactor*timeUnitConv;
    peak2troughRatio = abs(peakAmp / troughAmp);

    wf.mean = waveformMean;
    wf.halfPeakWidthMs = halfPeakWidthMs;
    wf.peak2troughDur = peak2troughDur;
    wf.peak2troughRatio = peak2troughRatio;
    neurons(i_neuron).wf = wf;
end

%% Clustering code

for i_neuron = 1:length(neurons)
    if ~isempty(neurons(i_neuron).wf)
        neuronWfAn(i_neuron,1) = neurons(i_neuron).wf.halfPeakWidthMs;
        neuronWfAn(i_neuron,2) = neurons(i_neuron).wf.peak2troughDur;
    end
end

% manually inspect for outlier
% figure
% plot(neuronWfAn(:,1),neuronWfAn(:,2),'*')


% GMM clustering
selUnits = neuronWfAn(:,1)< 0.4 & neuronWfAn(:,1) ~= 0; % exclude a two outliers

wfFeatures = neuronWfAn(selUnits,1:2);

options = statset('Display','final');
gm = fitgmdist(wfFeatures,2,'Options',options);
idx = cluster(gm,wfFeatures);
P = posterior(gm,wfFeatures);

% 99% confidence
idxP = nan(size(idx));
idxP(P(:,1)>0.99) = 1;
idxP(P(:,2)>0.99) = 2;

% GMM uses a random seem
% peak2troughDur should be smaller for interneuron
if mean(wfFeatures(idxP==1,2)) < mean(wfFeatures(idxP==2,2)), idxInt = idxP == 2; idxNonInt = idxP == 1; idxP(idxInt) = 1; idxP(idxNonInt) = 2; end

cluster1 = (idxP == 1);
cluster2 = (idxP == 2);

figure;
scatter(wfFeatures(cluster1,1),wfFeatures(cluster1,2),10,P(cluster1,1),'.')
hold on
scatter(wfFeatures(cluster2,1),wfFeatures(cluster2,2),10,P(cluster2,1),'.')
hold off
legend('Cluster 1','Cluster 2','Location','NW')
xlim([0 0.35]); ylim([0 0.6]);
xlabel('halfPeakWidtH(msec)'); ylabel('peak2trough(msec)'); 
clrmap = jet(80); colormap(clrmap(9:72,:))
title('STR'); 



selIdx = find(selUnits);
for i_selIdx = 1:length(selIdx)
    if idxP(i_selIdx) == 1
        neurons(selIdx(i_selIdx)).type = 'MSN';
    elseif idxP(i_selIdx) == 2
        neurons(selIdx(i_selIdx)).type = 'INT';
    else
        neurons(selIdx(i_selIdx)).type = 'NA';
    end
end


fprintf('\nTotal MSN neurons: %d', sum(strcmp({neurons.type},'MSN')));
fprintf('\nTotal INT neurons: %d', sum(strcmp({neurons.type},'INT')));
fprintf('\nTotal Excluded neurons: %d\n', sum(strcmp({neurons.type},'NA')));




%%
intStart = 0;
intEnd = 18; %latest switch time in dataset is at 18.5 seconds, therefore the interval has been expanded +/- .5 secs.
binSize = 0.2500; % in seconds
timep = (intStart+binSize:binSize:intEnd)';
interval_bins = intStart:binSize:intEnd;
for i_neuron = 1:length(neurons)
    if length(neurons(i_neuron).spikeTS) < 5000
        neurons(i_neuron).pTime = NaN;
        continue
    end    
    spike_trial = peSpike(neurons(i_neuron).spikeTS, neurons(i_neuron).events.switchTrialInit, timep');
    time_x = timep;
    a = histc(spike_trial',interval_bins); % time x trial
    b = a(1:end-1,:);
    t_firing_rate = reshape(b,[],1);
    t_times = repmat(time_x,size(b,2),1);

    %T_neuron = T_SPK(T_SPK.NeuronNum==i_neuron,:);
    T_neuron = table( t_firing_rate,t_times,   'VariableNames',{'FiringRate','Time'});

    lmNeuron = fitglme(T_neuron, 'FiringRate~Time');
    lmNeuronAnova = anova(lmNeuron);
    neurons(i_neuron).pTime = [lmNeuronAnova{2,5}];
    neurons(i_neuron).LMtime = lmNeuron;

end


%%

if load_tagging

intervalUP = 0.025; intervalDOWN = -0.025;

samplingRate = 30000;
sampleBin = 1/samplingRate;
wnarray = 0.001:0.001:0.02;
wnbinsbase = intervalDOWN:sampleBin:intervalUP;
baseRangeCnt = abs(intervalDOWN/sampleBin);

% for plotting
cmap = [0 0 1;1 0 0];
binSize = 0.0005; %0.5ms
histBins = intervalDOWN:binSize:intervalUP;

%
is_gen_fig = false;

%latency most important for now
for i_neuron = 1:length(neurons)
    
    spike_ts = neurons(i_neuron).spikeTS;
    eventTS = neurons(i_neuron).events.laser_ttl;
    if length(spike_ts) < 4000
        neurons(i_neuron).latency = 99;
        neurons(i_neuron).cr = .99;    
        neurons(i_neuron).tagged = 0;
        continue
    end
    
    edges = sort([eventTS+intervalDOWN;eventTS+intervalUP]);
    [~,~,bin] = histcounts(spike_ts,edges);
    Lia = ismember(bin,1:2:length(edges)); % select within trial(odd number index)
    
    trialTs = spike_ts(Lia);
    if length(trialTs)==0; continue; end; 
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
    wave = internalPL2Waves(neurons(i_neuron).pl2File,neurons(i_neuron).ch_name , neurons(i_neuron).ch_idx);
    % find waveforms index per trials, sperate by pre-/post- stimulus
    wf1 = mean(wave.Waves(preWfIdx,:),1);
    wf2 = mean(wave.Waves(postWfIdx,:),1);
    cr = xcorr(wf1,wf2,0,'coeff');
    latency = find(pwn<0.01,1);
    if isempty(latency)
       latency = 99; % need a placeholder value to advance through 
    end
    neurons(i_neuron).latency = latency;
    neurons(i_neuron).cr = cr;    
        
    if latency <= 5 
        neurons(i_neuron).tagged = 1; %if less than 10ms, it's tagged
    else
        neurons(i_neuron).tagged = 0;
    end

    if ~is_gen_fig
        continue
    end

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
    title(neurons(i_neuron).name);
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
    
    [~,name,~] = fileparts(neurons(i_neuron).pl2File);
    figStimFileName = [name neurons(i_neuron).name];
    set(gcf, 'Renderer', 'painters'),
    saveas(gcf,figStimFileName,'svg');
    %pause
    close(gcf)
    

end

end


if strcmp(data_type,'pharm')
    
    save('dataSt_pharm.mat', 'dataSt');
    save('neurons_pharm.mat', 'neurons');

elseif strcmp(data_type,'tagging')

    save('dataSt_tagging.mat', 'dataSt');
    save('neurons_tagging.mat', 'neurons');
end


end