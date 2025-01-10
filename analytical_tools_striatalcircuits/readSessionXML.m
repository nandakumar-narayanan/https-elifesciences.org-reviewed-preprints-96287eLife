function info = readSessionXML(dataDir)

if isfolder(dataDir)
    xmlfile = fullfile(dataDir, 'settings.xml');
elseif isfile(dataDir)
    xmlfile = dataDir;
    [dataDir] = fileparts(dataDir);
end
    
DOMnode = xmlread(xmlfile);
xRoot = DOMnode.getDocumentElement;

% get Creation date
RecordDateTime = char(xRoot.getElementsByTagName('INFO').item(0).getElementsByTagName('DATE').item(0).getTextContent);
GUIVersion = char(xRoot.getElementsByTagName('INFO').item(0).getElementsByTagName('VERSION').item(0).getTextContent);

xSignalChain = xRoot.getElementsByTagName('SIGNALCHAIN').item(0);
% get all processors
processors = struct([]);
xProcessors = xSignalChain.getElementsByTagName('PROCESSOR');
for j = 0:xProcessors.getLength-1
    processors(j+1,1).nodeId = str2double(xProcessors.item(j).getAttribute('NodeId'));
    processors(j+1,1).name = char(xProcessors.item(j).getAttribute('name'));
end

% get source Continous recording save settings
xRecControlInfo = xRoot.getElementsByTagName('CONTROLPANEL').item(0).getAttributes;
for j = 0:xRecControlInfo.getLength-1
    subFieldName = char(xRecControlInfo.item(j).getNodeName);
    RecordParams.(subFieldName) =  char(xRecControlInfo.item(j).getValue);
end

% Sources/Rhythm FPGA or Sources/File Reader, processor item 0 should be sources
xProcessor = xSignalChain.getElementsByTagName('PROCESSOR').item(0); 
ContinousChannels = struct([]);
ContinousRecordParams =  struct([]);
if xProcessor.getElementsByTagName('CHANNEL_INFO').getLength > 0 % Sources/Rhythm FPGA
    % get source Continous channel Info
    xCh = xProcessor.getElementsByTagName('CHANNEL_INFO').item(0).getElementsByTagName('CHANNEL');
    for j = 0:xCh.getLength-1
        ContinousChannels(j+1,1).name = char(xCh.item(j).getAttribute('name'));
        ContinousChannels(j+1,1).number = str2double(xCh.item(j).getAttribute('number'));
        ContinousChannels(j+1,1).gain = str2double(xCh.item(j).getAttribute('gain'));
        ContinousChannels(j+1,1).isLFP = strncmp(ContinousChannels(j+1,1).name,'CH',2);
    end
    % get source Continous recording parameters
    xEditorInfo = xProcessor.getElementsByTagName('EDITOR').item(0).getAttributes;
    for j = 0:xEditorInfo.getLength-1
        subFieldName = char(xEditorInfo.item(j).getNodeName);
        value = str2double(char(xEditorInfo.item(j).getValue));
        if ~isnan(value)
            ContinousRecordParams(1).(subFieldName) =  value;
        else
            ContinousRecordParams(1).(subFieldName) =  char(xEditorInfo.item(j).getValue);
        end
    end
end

% get Events channel Info
xCh = xProcessor.getElementsByTagName('EVENTCHANNEL');
EventChannels = struct;
for j = 0:xCh.getLength-1
    EventChannels(j+1,1).name = char(xCh.item(j).getAttribute('name'));
    EventChannels(j+1,1).number = str2double(xCh.item(j).getAttribute('number'));
end
EventChannelNum = xCh.getLength;

% get Spike recording parameters
electrodes = struct([]);
if any(strcmp({processors.name}, 'Filters/Spike Sorter')) || any(strcmp({processors.name}, 'Filters/Spike Detector'))
    if any(strcmp({processors.name}, 'Filters/Spike Sorter'))
        spkFilterIdx = find(strcmp({processors.name}, 'Filters/Spike Sorter'),1) -1;
    elseif any(strcmp({processors.name}, 'Filters/Spike Detector'))
        spkFilterIdx = find(strcmp({processors.name}, 'Filters/Spike Detector'),1) -1;
    end
    xElectrodes = xProcessors.item(spkFilterIdx).getElementsByTagName('ELECTRODE');
    for j = 0:xElectrodes.getLength-1
        electrodes(j+1,1).name = char(xElectrodes.item(j).getAttribute('name'));
        electrodes(j+1,1).numChannels = str2double(xElectrodes.item(j).getAttribute('numChannels'));
        electrodes(j+1,1).prePeakSamples = str2double(xElectrodes.item(j).getAttribute('prePeakSamples'));
        electrodes(j+1,1).postPeakSamples = str2double(xElectrodes.item(j).getAttribute('postPeakSamples'));
        for k = 0:xElectrodes.item(j).getElementsByTagName('SUBCHANNEL').getLength-1
            electrodes(j+1,1).number(k+1) = str2double(xElectrodes.item(j).getElementsByTagName('SUBCHANNEL').item(k).getAttribute('ch'));
            electrodes(j+1,1).thresh(k+1) = str2double(xElectrodes.item(j).getElementsByTagName('SUBCHANNEL').item(k).getAttribute('thresh'));
            electrodes(j+1,1).isActive(k+1) = str2double(xElectrodes.item(j).getElementsByTagName('SUBCHANNEL').item(k).getAttribute('isActive'));
        end
    end
end

sampleRateIdx = [ 1 1.25 1.5 2 2.5 3 3.33 4 5 6.25 8 10 12.5 15 20 25 30 ];
info.dataPath = dataDir; %
info.RecordDateTime = RecordDateTime;
info.GUIVersion = GUIVersion;
info.Processors = processors;
info.RecordParams = RecordParams;

% temporary fix, eg E:\PlexonData\OpenEphys\TEST_2017-03-01_11-31-48_16
% OpenEphys 0.4.1 bug: sometime cnt xml have wrong info.
if length(ContinousChannels) == 19 % 16 WB + 3 AUX
    if all(strncmp('CH',{ContinousChannels.name},2)) && sum(strncmp('AUX',{ContinousChannels.name},3)) == 0 % no AUX channels
        for auxId = 17:19
            ContinousChannels(auxId).name = sprintf('AUX%d',auxId-16);
            ContinousChannels(auxId).isLFP = false;
        end
    end
end

% 5/16/2018 fix for OpenEphys > 0.4.3 bug - XML do not have proper EVENT channels
% Force EVENT channel info, since they are always there.
if length(EventChannels) == 1
    for j = 0:15
        EventChannels(j+1,1).name = num2str(j);
        EventChannels(j+1,1).number = j;
    end
    EventChannelNum = 16;
end

info.ContinousChannels = ContinousChannels;
info.SpikeChannels = electrodes; 
info.EventChannels = EventChannels;

info.ContinousRecordParams = ContinousRecordParams;
if ~isempty(ContinousChannels)
    info.ContinousChannelNum = length(ContinousChannels);
    info.LFPChannelNum = sum(cell2mat({info.ContinousChannels.isLFP}));
    info.Gain = info.ContinousChannels(1).gain;
    info.SampleRate = sampleRateIdx(info.ContinousRecordParams.SampleRate)*1000;
end
if ~isempty(electrodes) 
    info.SpikeRecordParam.NumSpike = electrodes(1).prePeakSamples + electrodes(1).postPeakSamples;
    info.SpikeRecordParam.NumPreSpike = electrodes(1).prePeakSamples;
end
info.SpikeChannelNum = length(electrodes);
info.EventChannelNum = EventChannelNum;

info.RecordEngine = info.RecordParams.recordEngine;
info.RecordDT = datetime(info.RecordDateTime,'InputFormat','d MMM yyyy H:mm:ss');
info.ADCRefVolt = 1.225;
info.ADCGain = 192;
info.CoeffToConvertToUnits = info.ADCRefVolt/info.ADCGain/double(intmax('int16'));
