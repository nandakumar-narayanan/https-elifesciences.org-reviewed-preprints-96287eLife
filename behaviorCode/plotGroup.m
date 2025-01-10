function GroupAnSt = plotGroup(TrialAnSt, group)

if isempty(group)
    return
end

GroupAnSt = struct;

groupFields = fieldnames(group);
numOfGroup = size(groupFields,1);
clear groupDataTableIdx
for groupIdx = 1:numOfGroup
    dayNanimalGroup = group.(groupFields{groupIdx});
    dayAnimal_cnt = 1;
    for i_dayAnimal = 1:size(dayNanimalGroup,2)
        dayNanimal = dayNanimalGroup{i_dayAnimal};
        date = dayNanimal{1};
        animal = dayNanimal{2};
        if isempty(date)
            for dayIdx = 1:size(TrialAnSt,2)
                if ~isempty(TrialAnSt(dayIdx).(animal))
                    groupDataTableIdx{dayAnimal_cnt,groupIdx}={animal,dayIdx};
                    dayAnimal_cnt = dayAnimal_cnt+1;
                end
            end
            
        else
            dateNum = date2TableNum(animal, date, TrialAnSt);
            groupDataTableIdx{dayAnimal_cnt,groupIdx}={animal,dateNum};
            dayAnimal_cnt = dayAnimal_cnt+1;
        end
       
    end
end


% Group pooled data analysis
% normalized mean, smoothed mean, pooled raw timing
centers = 0:0.1:20;
for groupIdx = 1:numOfGroup
    % within group, find unique animal , find all date of unique animal
    cell_data_idx = cellfun(@(x) ~isempty(x), groupDataTableIdx(:,groupIdx));
    ids = cellfun(@(x) x{1},groupDataTableIdx(cell_data_idx,groupIdx),'UniformOutput',false);
    uni_ids = unique(ids);
%     cdf_switch_depart = zeros(numel(centers)-1,size(uni_ids,1));
%     pdf_switch_depart = zeros(numel(centers)-1,size(uni_ids,1));
    cdf_switch = [];
    pdf_switch = [];
    for uniNameIdx = 1:size(uni_ids,1)
        groupMember = find(strcmp(uni_ids{uniNameIdx},ids));
        totSwitchTrDe = []; totSwitchTrAr = [];
        pooled_short_rsp = []; pooled_long_rsp = [];
        totTraversalTimeSwitch = []; totRspDuration = []; totReward = [];
        % nose poke duration
        % traversal time between the short and long response port        
        for groupMemberIdx = 1:length(groupMember)
            i_day = groupDataTableIdx{groupMember(groupMemberIdx),groupIdx}{2};
            animal = groupDataTableIdx{groupMember(groupMemberIdx),groupIdx}{1};
            tmp_struct = TrialAnSt(i_day).(animal);
            if numel([tmp_struct.SwitchDepart]) < 5 % this restriction was only applied to opto data?
                continue
            end
            if ~isfield(tmp_struct,'trialType') && isfield(tmp_struct,'opto')
                [tmp_struct.trialType] = deal(tmp_struct.opto);
                TrialAnSt(i_day).(animal) = tmp_struct;
            elseif ~isfield(tmp_struct,'trialType')
                [tmp_struct.trialType] = deal(0);
                TrialAnSt(i_day).(animal) = tmp_struct;
            end
            valid_left_idx = cellfun(@(x) numel(x), {tmp_struct.leftRelTimeTrial}) == cellfun(@(x) numel(x), {tmp_struct.leftRspTimeTrial});
            valid_right_idx = cellfun(@(x) numel(x), {tmp_struct.rightRelTimeTrial}) == cellfun(@(x) numel(x), {tmp_struct.rightRspTimeTrial});            
            uni_type = unique([tmp_struct.trialType]);
            SwitchTrDe = []; SwitchTrAr = []; TraversalTimeSwitch = []; RspDuration = []; Reward = [];
            for i_type = 1:length(uni_type)
                trial_type_idx = uni_type(i_type) == [tmp_struct.trialType];
                switch_idx = cellfun(@(x) ~isempty(x), {tmp_struct.SwitchDepart});
                % Pulls out the long trials that should contain a switch
                sd = [tmp_struct(trial_type_idx).SwitchDepart];
                sa = [tmp_struct(trial_type_idx).SwitchArrival];
                SwitchTrDe(1:numel(sd),i_type) = sd;
                SwitchTrAr(1:numel(sd),i_type) = sa;
                TraversalTimeSwitch(1:numel(sa),i_type) = sa-sd;

                sa = [tmp_struct(trial_type_idx & valid_left_idx).leftRelTimeTrial] - [tmp_struct(trial_type_idx & valid_left_idx).leftRspTimeTrial];
                sb = [tmp_struct(trial_type_idx & valid_right_idx).rightRelTimeTrial] - [tmp_struct(trial_type_idx & valid_right_idx).rightRspTimeTrial];
                sc = [sa sb];
                RspDuration(1:numel(sc),i_type) = sc;
                sa = [tmp_struct(trial_type_idx & switch_idx).reward];
                Reward(1:numel(sa),i_type) = sa;
            end
            totSwitchTrDe = [totSwitchTrDe; SwitchTrDe]; totSwitchTrAr = [totSwitchTrAr; SwitchTrAr];
            totTraversalTimeSwitch  = [totTraversalTimeSwitch; TraversalTimeSwitch]; 
            totRspDuration = [totRspDuration; RspDuration]; totReward = [totReward; Reward];
            pooled_short_rsp = [pooled_short_rsp, [tmp_struct.ShortRsp]];
            pooled_long_rsp = [pooled_long_rsp, [tmp_struct.LongRsp]];
        end
        totSwitchTrDe = totSwitchTrDe(:,~all(totSwitchTrDe == 0));
        totSwitchTrAr = totSwitchTrAr(:,~all(totSwitchTrAr == 0));

        stat = struct;
        for i_var = 1:2
            num_type = size(totSwitchTrDe,2);
            for i_type = 1:num_type
                if i_var == 1
                    st = totSwitchTrDe;
                elseif i_var ==2
                    st = totSwitchTrAr;
                end
                if isempty(nonzeros(st(:,i_type)))
                    continue
                end
                % Compiles all departure and arrival across animals and days
                [N1,~] = histcounts(nonzeros(st(:,i_type)),centers,'Normalization','cdf');
                [f1,~] = ksdensity(nonzeros(st(:,i_type)),centers(2:end),'Bandwidth',0.7);
                cdf_switch(:,i_type,i_var,uniNameIdx) = N1;
                pdf_switch(:,i_type,i_var,uniNameIdx) = f1;
                % Calculates and builds summary stats table for each animal across days
                [stat(i_type).avg, STD, stat(i_type).N, SEM, stat(i_type).med, Q1, Q2, IntQRange, Limen, WR] = getStat(st(:,i_type));

            end
        end         
        totSwitchTrDe(totSwitchTrDe == 0) = NaN;
        totTraversalTimeSwitch(totTraversalTimeSwitch == 0) = NaN;
        totRspDuration(totRspDuration==0) = NaN;
        totReward(totReward==0) = NaN;
        GroupAnSt.(groupFields{groupIdx}).ids{uniNameIdx} = uni_ids{uniNameIdx};
        GroupAnSt.(groupFields{groupIdx}).totReward{uniNameIdx} = totReward;
        GroupAnSt.(groupFields{groupIdx}).totRspDuration{uniNameIdx} = totRspDuration;
        GroupAnSt.(groupFields{groupIdx}).totTraversalTime{uniNameIdx} = totTraversalTimeSwitch;
        GroupAnSt.(groupFields{groupIdx}).totSwitchTrDe{uniNameIdx} = totSwitchTrDe;
        GroupAnSt.(groupFields{groupIdx}).totSwitchTrAr{uniNameIdx} = totSwitchTrAr;
        GroupAnSt.(groupFields{groupIdx}).pooled_short_rsp{uniNameIdx} = pooled_short_rsp;
        GroupAnSt.(groupFields{groupIdx}).pooled_long_rsp{uniNameIdx} = pooled_long_rsp;          
    end
    GroupAnSt.(groupFields{groupIdx}).cdf_switch = cdf_switch;
    GroupAnSt.(groupFields{groupIdx}).pdf_switch = pdf_switch;

end

return

cmap = lines(numOfGroup);
trial_time = centers(2:end);
figure('Units', 'Normalized', 'OuterPosition', [0.01 0.02 0.98, 0.96])
subplot(1,2,1)
for groupIdx = 1:numOfGroup
    data = GroupAnSt.(groupFields{groupIdx}).cdf_switch_depart;
    cdf50 = [];
    for i_animal = 1:size(data,2)
        idx1 = find(data(:,i_animal)<0.5,1,'last');
        idx2 = find(data(:,i_animal)>=0.5,1,'first');
        cdf50(i_animal) = interp1(data([idx1 idx2],i_animal),trial_time([idx1 idx2]),0.5);
    end
    plotband(trial_time,mean(data,2),std(data,[],2)/sqrt(size(data,2)),cmap(groupIdx,:))
    text(0.25, 0.95-0.05*groupIdx, sprintf('%s: %0.2f +- %0.2f',groupFields{groupIdx},mean(cdf50),std(cdf50)/sqrt(size(data,2))), 'Color', 'k', 'FontSize', 12);
    data_mean = mean(data,2);
    idx1 = find(data_mean<0.5,1,'last');
    idx2 = find(data_mean>=0.5,1,'first');
    cdf50_mean = interp1(data_mean([idx1 idx2]),trial_time([idx1 idx2]),0.5);
    plot([0 cdf50_mean cdf50_mean], [0.5 0.5 0],'Color',cmap(groupIdx,:))
    GroupAnSt.(groupFields{groupIdx}).cdf50 = cdf50;
end
hold off
grid on
xlabel('Trial Duration'); ylabel({'Cumulative Probability of a Switch'});
set(gca, 'Color', 'none', 'box','off'); 
ylim([0 1])


subplot(1,2,2)
hold on
for groupIdx = 1:numOfGroup
    data = GroupAnSt.(groupFields{groupIdx}).pdf_switch_depart;
    plotband(centers(2:end),mean(data,2),std(data,[],2)/sqrt(size(data,2)),cmap(groupIdx,:))
end
hold off
xlabel('Trial Duration'); ylabel('Probability Distrubtion Function');
set(gca, 'Color', 'none', 'box','off'); 


for groupIdx = 1:numOfGroup
    group_f = groupFields{groupIdx};
    data_cell = GroupAnSt.(group_f).totSwitchTrDe;

    subplot_idx = 1;
    numplot_y = ceil(sqrt(length(data_cell))); %  positions and sizes for subplot
    numplot_x = ceil(length(data_cell)/numplot_y);
    figure('Units', 'Normalized', 'OuterPosition', [0.01 0.02 0.98, 0.96])
    stat = [];
    for i_animal = 1:length(data_cell)
        data = data_cell{i_animal}';
        [h,p,stats] = chi2gof(data);
        pd = fitdist(data,'Normal');
        fwhm = pd.sigma * 2*sqrt(2*log(2));
        x_pdf = [0:0.1:14];
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
        title(group_f)
        xlabel('Trial Duration'); ylabel({'Switch Departure PDF'});
        set(gca, 'Color', 'none', 'box','off');     
        stat(1:6,i_animal) = [mean(data), median(data), p, fwhm, pd.sigma, iqr(pd)];
        subplot_idx = subplot_idx+1;
    end
    GroupAnSt.(groupFields{groupIdx}).stats = stat;
end

end

function dateNum = date2TableNum(animal, date, anSt)
if isfield(anSt,animal)
    for dayIdx = 1:size(anSt,2)
        if ~isempty(anSt(dayIdx).(animal))
            if strcmp(date,anSt(dayIdx).(animal)(1).mpc.StartDate)
                dateNum = dayIdx;
            end
        end
    end
else
    return
end
end


% builds summary stats table
function [tAvg, tSTD, tN, tSEM, tMedian, tQ1, tQ2, tIntQRange, tLimen, tWR] = getStat(x)
    x = nonzeros(x);
    tAvg = mean(x); 
    tSTD = std(x); 
    tN = numel(x); 
    tSEM = tSTD/(sqrt(tN));
    tMedian = median(x); 
    tQ1 = quantile(x, 0.25); 
    tQ2 = quantile(x, 0.75); 
    tIntQRange = tQ2 - tQ1; tLimen = tIntQRange/2; tWR = tLimen/tMedian;
    % bstats = [tAvg, tSTD, tN, tSEM, tMedian, tQ1, tQ2, tIntQRange, tLimen, tWR];
end
