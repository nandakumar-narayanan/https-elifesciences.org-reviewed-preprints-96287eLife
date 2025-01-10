function events = get_mpc_bin_event_pl2_v2(pl2File)

% Binary event
event_record_type = {...
'1 0 0 0 0', 'SESSION START'; ...
'0 1 0 0 0', 'Back Nosepoke Light ON'; ...
'0 0 1 0 0', 'Back Nosepoke Light OFF'; ...
'0 0 0 1 0', 'ITI START'; ...
'0 0 0 0 1', 'ITI END'; ...
'1 1 0 0 0', 'CUES ON'; ...
'1 0 1 0 0', 'CUES OFF'; ...
'1 0 0 1 0', 'REWARD DISPENSE'; ...
'1 0 0 0 1', 'TRIAL END'; ...
'1 1 1 0 0', 'LEFT RESPONSE'; ...
'1 1 0 1 0', 'LEFT RELEASE'; ...
'1 1 0 0 1', 'RIGHT RESPONSE'; ...
'1 1 1 1 0', 'RIGHT RELEASE'; ...
'1 1 1 1 1', 'BACK RESPONSE'; ...
'0 1 1 0 0', 'BACK RELEASE'; ...
'0 1 0 1 0', 'REWARD RESPONSE/CHECKING FEEDER'; ...
'0 1 0 0 1', 'REWARD RELEASE/REMOVE HEAD FROM FEEDER'...; ...
};

% Z^zRPBreakBACK; Z^zPlexRESPBACK; -> Z^zPlexCUESON;, Z^zPlexNPLBACKOFF;

% convert binary event to decimal
event_type_dec = cellfun(@(x) bin2dec(fliplr(x)), event_record_type(:,1));
event_record_type(1:17,3) = num2cell(event_type_dec);

% read pl2 event data
% pl2File = datapath;
pl2Info = PL2GetFileIndex(pl2File);
% get events
eventChannels = pl2Info.EventChannels;
%numEvnetChannels = numel(eventChannels);

ttlEventIdx = cellfun(@(x) strncmp(x.Name,'EVT',3), eventChannels) & cellfun(@(x) x.NumEvents > 0, eventChannels);
eventChannel = cellfun(@(x) x.Channel, eventChannels);
ttlChannel = eventChannel(ttlEventIdx);
eventChIdx = find(ttlEventIdx);
events = struct;
if ~length(ttlChannel) == 0
    for i=1:length(ttlChannel)
        evt = PL2EventTs(pl2File, eventChIdx(i));
        events.(sprintf('evt%d',ttlChannel(i))) = evt.Ts;
    end
else
    events = struct();
end

eventData = struct;
evt_data_idx = 1;
for i=1:length(ttlChannel)
    evt = PL2EventTs(pl2File, eventChIdx(i));
    evt_num = length(evt.Ts);
    eventData.channel(evt_data_idx:evt_data_idx+evt_num-1) = eventChIdx(i);
    eventData.timestamps(evt_data_idx:evt_data_idx+evt_num-1) = evt.Ts;
    if eventChIdx(i) > 8
        eventData.eventId(evt_data_idx:evt_data_idx+evt_num-1) = 1;
    else
        eventData.eventId(evt_data_idx:evt_data_idx+evt_num-1) = 0;
    end
    evt_data_idx  = evt_data_idx + evt_num;
end
[~,sort_idx] = sort(eventData.timestamps);
eventData.timestamps = eventData.timestamps(sort_idx)';
eventData.channel = eventData.channel(sort_idx)';
eventData.eventId = eventData.eventId(sort_idx)';
eventData.channel(eventData.channel>8) = eventData.channel(eventData.channel>8)-8;
eventData.channel = eventData.channel-1;
sample_rate = 30000;
eventData.timestamps = round(eventData.timestamps.*sample_rate);
ch_uni = unique(eventData.channel);
% event channel 4 6 7 is not used for binary event % event channel 4 used for external TTL(Laser/LED)
ch_uni = double(ch_uni(ch_uni<6 & ch_uni~= 4)); 

% get event channel rising/falling 
all_ch_evt = zeros(sum(eventData.channel<6 & eventData.channel~= 4 & eventData.eventId == 1),4);
ts_idx = 1;
for i_ch = 1:length(ch_uni)
    rising_ts = eventData.timestamps(eventData.eventId == 1 & eventData.channel == ch_uni(i_ch));
    falling_ts = eventData.timestamps(eventData.eventId == 0 & eventData.channel == ch_uni(i_ch));
    all_evt = [repmat(ch_uni(i_ch), numel(falling_ts), 1) double(falling_ts) double(rising_ts) double(rising_ts) - double(falling_ts)];
    all_ch_evt(ts_idx:ts_idx+size(all_evt,1)-1,:) = all_evt;
    ts_idx = ts_idx+size(all_evt,1);
end
[~,skey]= sort(all_ch_evt(:,2));
all_ch_evt = all_ch_evt(skey,:); % sort by rising tick

all_ch_evt(all_ch_evt(:,1) == 5,1) = 4; % replace channel 5 to 4
all_ch_evt(:,1) = all_ch_evt(:,1)+1; % channel range from 0-4 to 1-5

% find splitted time-stamp by 1 tick and replace them with correct tick
uni_ts = unique(all_ch_evt(:,2:3));
split_ts_idx = find(diff(uni_ts) == 1);
for i_ts = 1:length(split_ts_idx)
    split_idx = all_ch_evt(:,2) == uni_ts(split_ts_idx(i_ts)+1);
    all_ch_evt(split_idx,2) = uni_ts(split_ts_idx(i_ts));
    split_idx = all_ch_evt(:,3) == uni_ts(split_ts_idx(i_ts)+1);
    all_ch_evt(split_idx,3) = uni_ts(split_ts_idx(i_ts));    
end

% track event channel status
eventType = '00000';
sample_rate = 30000;
for i_evt = 1:length(uni_ts)
    ts_on_idx = all_ch_evt(uni_ts(i_evt) == all_ch_evt(:,2),1);
    ts_off_idx = all_ch_evt(uni_ts(i_evt) == all_ch_evt(:,3),1);
    eventType(ts_on_idx) = '1';
    eventType(ts_off_idx) = '0';
    event_tt(i_evt,1) = double(uni_ts(i_evt))/sample_rate;
    event_tt(i_evt,2) = bin2dec(fliplr(eventType));
    event_tt(i_evt,3) = double(uni_ts(i_evt));
end

events = struct;
uni_type = unique(event_tt(:,2));
uni_type = uni_type(uni_type>0);
for i_type = 1:length(uni_type)
    if any(event_type_dec == uni_type(i_type))
        events.(sprintf('evt%d',uni_type(i_type))).ts = event_tt(uni_type(i_type) == event_tt(:,2),1);
        events.(sprintf('evt%d',uni_type(i_type))).type = event_record_type{event_type_dec == uni_type(i_type),2};
    end
end

led_ttl = eventData.timestamps(eventData.eventId == 1 & eventData.channel == 4);
if ~isempty(led_ttl)
    events.evt32.ts = double(led_ttl)./sample_rate;
    events.evt32.type = 'LED TTL';
end

% verbose output
eventnames = fieldnames(events);
for i_type = 1:length(eventnames)
    fprintf('%s: %s - %d events\n', eventnames{i_type}, events.(eventnames{i_type}).type,length(events.(eventnames{i_type}).ts))
end


