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
    for i_dayAnimal = 1:size(dayNanimalGroup,2)
        dayNanimal = dayNanimalGroup{i_dayAnimal};
        date = dayNanimal{1};
        animal = dayNanimal{2};
        dateNum = date2TableNum(animal, date, TrialAnSt);
        groupDataTableIdx{i_dayAnimal,groupIdx}={animal,dateNum};
    end
end

% Group pooled data analysis
% normalized mean, smoothed mean, pooled raw timing
centers = 0:0.1:20;
edges_ansr = linspace(0,20,21);
for groupIdx = 1:numOfGroup
    % within group, find unique animal , find all date of unique animal
    cell_data_idx = cellfun(@(x) ~isempty(x), groupDataTableIdx(:,groupIdx));
    ids = cellfun(@(x) x{1},groupDataTableIdx(cell_data_idx,groupIdx),'UniformOutput',false);
    uni_ids = unique(ids);
    cdf_switch_depart = zeros(numel(centers)-1,size(uni_ids,1));
    pdf_switch_depart = zeros(numel(centers)-1,size(uni_ids,1));
    pdf_shortRsp = zeros(numel(centers)-1,size(uni_ids,1));
    pdf_longRsp = zeros(numel(centers)-1,size(uni_ids,1));
    ansr_switch_depart = zeros(numel(edges_ansr),size(uni_ids,1)); %%%%% - Matt would recommend taking this out, I would use a PDF of swtich times instead
    
    %LOOPS TO SET UP DATA FOR EACH ANIMAL + DAY
    for uniNameIdx = 1:size(uni_ids,1)
        groupMember = find(strcmp(uni_ids{uniNameIdx},ids));
        totSwitchTrDe = []; totSwitchTrAr = [];
        pooled_short_rsp = []; pooled_long_rsp = []; pooled_long_rew=[]; %unsure if we want to include reward here
        pooled_rt = []; pooled_short_rel = []; pooled_long_rel= [];
        short_rerponse_durations = []; long_response_durations = [];
        pooled_tt = [];
        
        for groupMemberIdx = 1:length(groupMember)
            i_day = groupDataTableIdx{groupMember(groupMemberIdx),groupIdx}{2};
            animal = groupDataTableIdx{groupMember(groupMemberIdx),groupIdx}{1};
            tmp_struct = TrialAnSt(i_day).(animal);
            % Pulls out the long trials that should contain a switch
            SwitchTrDe = [tmp_struct.SwitchDepart]; SwitchTrAr = [tmp_struct.SwitchArrival];
            totSwitchTrDe = [totSwitchTrDe, SwitchTrDe]; totSwitchTrAr = [totSwitchTrAr, SwitchTrAr];
            pooled_short_rsp = [pooled_short_rsp, [tmp_struct.ShortRsp]];
            pooled_short_rel = [pooled_short_rel, [tmp_struct.ShortRel]];
            pooled_long_rsp = [pooled_long_rsp, [tmp_struct.LongRsp]];
            pooled_long_rel = [pooled_long_rel, [tmp_struct.LongRel]];
            % this is the trial on which medPC was terminated - MSN is
            % written to terminate on beam break.
            %should be either one or the other for every mouse unless
            %manually terminated at 90 minutes.
                if length(pooled_short_rsp) > length(pooled_short_rel)
                    pooled_short_rsp = pooled_short_rsp(1:end-1);
                elseif length(pooled_long_rsp) > length(pooled_long_rel)
                    pooled_long_rsp = pooled_long_rsp(1:end-1);
                end
            short_response_durations = [pooled_short_rel] - [pooled_short_rsp];
            long_response_durations = [pooled_long_rel] - [pooled_long_rsp];
            pooled_tt = totSwitchTrAr - totSwitchTrDe;
            pooled_long_rew = [pooled_long_rew, [tmp_struct.reward]];
            pooled_rt = [pooled_rt, [tmp_struct.RT]];
        end
        [N,edges] = histcounts(totSwitchTrDe,centers,'Normalization','cdf');
        [f,xi] = ksdensity(totSwitchTrDe,centers(2:end),'Bandwidth',0.6);
        cdf_switch_depart(:,uniNameIdx) = N;
        pdf_switch_depart(:,uniNameIdx) = f;
        
        [fsh,xi] = ksdensity(pooled_short_rsp,centers(2:end),'Bandwidth',0.6);
        pdf_shortRsp(:,uniNameIdx) = fsh;
        
        [flo,xi] = ksdensity(pooled_long_rsp,centers(2:end),'Bandwidth',0.6);
        pdf_longRsp(:,uniNameIdx) = flo;
        
        
        %% If you want to keep this, I would recommend changing the smooth - I know that I was oversmoothing data early on
        ansr_hist = hist(totSwitchTrDe, edges_ansr);
        ansr_norm_hist = ((ansr_hist - (min(ansr_hist)))) / ((max(ansr_hist)) - (min(ansr_hist)));
        ansr_norm_hist = smooth(ansr_norm_hist, 1.25);
        ansr_switch_depart(:, uniNameIdx) = ansr_norm_hist;
        %%
        
        GroupAnSt.(groupFields{groupIdx}).totSwitchTrDe{uniNameIdx} = totSwitchTrDe;
        GroupAnSt.(groupFields{groupIdx}).median_switch_dep{uniNameIdx} = median(totSwitchTrDe);
        GroupAnSt.(groupFields{groupIdx}).mean_switch_dep{uniNameIdx} = mean(totSwitchTrDe);
        GroupAnSt.(groupFields{groupIdx}).totSwitchTrAr{uniNameIdx} = totSwitchTrAr;
        GroupAnSt.(groupFields{groupIdx}).median_switch_arr{uniNameIdx} = median(totSwitchTrAr);
        GroupAnSt.(groupFields{groupIdx}).mean_switch_arr{uniNameIdx} = mean(totSwitchTrAr);
        GroupAnSt.(groupFields{groupIdx}).pooled_short_rsp{uniNameIdx} = pooled_short_rsp;
        GroupAnSt.(groupFields{groupIdx}).pooled_long_rsp{uniNameIdx} = pooled_long_rsp;
        GroupAnSt.(groupFields{groupIdx}).pooled_long_rew{uniNameIdx} = pooled_long_rew;
        GroupAnSt.(groupFields{groupIdx}).median_short_rsp{uniNameIdx} = median(pooled_short_rsp);
        GroupAnSt.(groupFields{groupIdx}).median_long_rsp{uniNameIdx} = median(pooled_long_rsp);
        GroupAnSt.(groupFields{groupIdx}).pooled_rt{uniNameIdx} = pooled_rt;
        GroupAnSt.(groupFields{groupIdx}).median_rt{uniNameIdx} = median(pooled_rt);
        GroupAnSt.(groupFields{groupIdx}).pooled_tt{uniNameIdx}= pooled_tt;
        GroupAnSt.(groupFields{groupIdx}).median_tt{uniNameIdx} = median(pooled_tt);
        GroupAnSt.(groupFields{groupIdx}).srd{uniNameIdx} = short_response_durations;
        GroupAnSt.(groupFields{groupIdx}).median_srd{uniNameIdx} = median(short_response_durations);
        GroupAnSt.(groupFields{groupIdx}).lrd{uniNameIdx} = long_response_durations;
        GroupAnSt.(groupFields{groupIdx}).median_lrd{uniNameIdx} = median(long_response_durations);
    end
    GroupAnSt.(groupFields{groupIdx}).cdf_switch_depart = cdf_switch_depart;
    GroupAnSt.(groupFields{groupIdx}).pdf_switch_depart = pdf_switch_depart;
    GroupAnSt.(groupFields{groupIdx}).pdf_shortRsp = pdf_shortRsp;
    GroupAnSt.(groupFields{groupIdx}).pdf_longRsp = pdf_longRsp;
    GroupAnSt.(groupFields{groupIdx}).ansr_switch_depart = ansr_switch_depart;
    
end

%COLORSCHEME FOR FIGURES
% black (w/in subjects control); blue (D1-cre experimental); red (D2-cre
% experimental); grey (control fluorophore)

if ismember("Sulpiride", groupFields)
    cmap = [0 0 0; .9 .1 .1; .1 .1 .75; .5 .5 .5];
elseif ismember("SCH23390", groupFields)
    cmap = [0 0 0; .1 .1 .75; .5 .5 .5];
end

colormap = cmap;

trial_time = centers(2:end);
figure();

for groupIdx = 1:numOfGroup
    data = GroupAnSt.(groupFields{groupIdx}).cdf_switch_depart;
    cdf50 = [];
    for i_animal = 1:size(data,2)
        idx1 = find(data(:,i_animal)<0.5,1,'last');
        idx2 = find(data(:,i_animal)>=0.5,1,'first');
        cdf50(i_animal) = interp1(data([idx1 idx2],i_animal),trial_time([idx1 idx2]),0.5);
    end
    plotband(trial_time,mean(data,2),std(data,[],2)/sqrt(size(data,2)),cmap(groupIdx,:))
    text(0.25, 1.025-0.075*groupIdx, sprintf('%s: %0.2f +- %0.2f',groupFields{groupIdx},mean(cdf50),std(cdf50)/sqrt(size(data,2))), 'Color', cmap(groupIdx,:), 'FontSize', 12);
    data_mean = mean(data,2);  % mean on CDF and statistical mean are different.
    idx1 = find(data_mean<0.5,1,'last');
    idx2 = find(data_mean>=0.5,1,'first');
    cdf50_mean = interp1(data_mean([idx1 idx2]),trial_time([idx1 idx2]),0.5);
    plot([0 cdf50_mean cdf50_mean], [0.5 0.5 0],'Color',cmap(groupIdx,:))
    GroupAnSt.(groupFields{groupIdx}).cdf50 = cdf50;
end
hold off
grid on
xlabel('Trial Duration'); ylabel({'Switch Liklihood'});
text(1, 1.05-0.1*(groupIdx+1), sprintf('%s = %d/group', 'n', length(data(1,:))), 'Color', 'k', 'FontSize', 16);
set(gca, 'Color', 'none', 'box','off');
ylim([0 1])


%% Same comment as above, I would recommend taking this out.
figure()

hold on
for groupIdx = 1:numOfGroup
    data = GroupAnSt.(groupFields{groupIdx}).ansr_switch_depart;
    plotband(edges_ansr,smooth(mean(data,2),3),std(data,[],2)/sqrt(size(data,2)),cmap(groupIdx,:))
end
hold off
xlabel('Trial Duration (s)'); ylabel('Normalized Time to Switch'); xline(6, 'LineWidth', 2);
yticks(0:0.25:1); ylim([0 0.75]);
set(gca, 'Color', 'none', 'box','off');
%%

for groupIdx = 1:numOfGroup
    group_f = groupFields{groupIdx};
    data_cell = GroupAnSt.(group_f).totSwitchTrDe;
    numplot_y = ceil(sqrt(length(data_cell))); %  positions and sizes for subplot
    numplot_x = ceil(length(data_cell)/numplot_y);
    stat = [];
    sw_stats = {};
    for i_animal = 1:length(data_cell)
        data = data_cell{i_animal}';
        %suppress warnings
        warning('off')
        [h,p,stats] = chi2gof(data);
        warning('on')
        pd = fitdist(data,'Normal');
        fwhm = pd.sigma * 2*sqrt(2*log(2));
        x_pdf = [0:0.1:24];
        y = pdf(pd,x_pdf);
        cv = ((std(data))/(mean(data)))*100;
        Q1 = quantile(data, [0.25]); Q2 = quantile(data, [0.75]); IQR = Q2 - Q1;
        stat(1:8,i_animal) = [mean(data), median(data), p, fwhm, pd.sigma, iqr(pd), cv, IQR];
        sw_stats{i_animal} = swft(data); %SHAPIRO WILK LOOKING FOR NORMALITY
    end
    GroupAnSt.(groupFields{groupIdx}).stats = stat;
    GroupAnSt.(groupFields{groupIdx}).sw_stats = sw_stats;
    
end

comparison_name = fieldnames(GroupAnSt);
comparison_name = strcat(comparison_name{1}, " v. ", comparison_name{2});

%paired t-test for mean switch times
[h_switch,p_switch,ci_switch,stats_switch] = ttest([GroupAnSt.(groupFields{1}).mean_switch_dep{:}], [GroupAnSt.(groupFields{2}).mean_switch_dep{:}]);%

fprintf('\n %s - Mean Switch Latency (seconds): SALINE %0.2f(%.2f-%.2f) vs DRUG: %0.2f(%.2f-%.2f), t(%0.2f)=%2.1f, p = %2.5f \n',...
    comparison_name,...
    mean([GroupAnSt.(groupFields{1}).mean_switch_dep{:}]), quantile([GroupAnSt.(groupFields{1}).mean_switch_dep{:}],0.25),quantile([GroupAnSt.(groupFields{1}).mean_switch_dep{:}], .75),...
    mean([GroupAnSt.(groupFields{2}).totSwitchTrDe{:}]), quantile([GroupAnSt.(groupFields{2}).totSwitchTrDe{:}],0.25),quantile([GroupAnSt.(groupFields{2}).mean_switch_dep{:}], .75),...
    stats_switch.tstat, stats_switch.df ,p_switch);

fprintf('\n Effect size - Cohend = %0.3f',...
    abs(cohend([GroupAnSt.(groupFields{2}).mean_switch_dep{:}], [GroupAnSt.(groupFields{1}).mean_switch_dep{:}])))

% PLOTTING TO SHOW NORMAL DISTRIBUTION OVERALL

[sw_dep_histw_saline, vinterval,wt_ex_sw_dep_saline] = histwcc([GroupAnSt.(groupFields{1}).totSwitchTrDe], 0:0.40:24);
[sw_dep_histw_drug, vinterval,wt_ex_sw_dep_drug] = histwcc([GroupAnSt.(groupFields{2}).totSwitchTrDe], 0:0.40:24);

[pdf_sw_dep_saline,~] = ksdensity(cell2mat([GroupAnSt.(groupFields{1}).totSwitchTrDe]),centers(2:end),'Bandwidth',0.6,'Weights',wt_ex_sw_dep_saline);
[pdf_sw_dep_drug,~] = ksdensity(cell2mat([GroupAnSt.(groupFields{2}).totSwitchTrDe]),centers(2:end),'Bandwidth',0.6,'Weights',wt_ex_sw_dep_drug);

figure();
hold on;
histogram('BinEdges',vinterval,'BinCounts',sw_dep_histw_drug,'LineStyle','none','FaceColor',colormap(2,:))
plot(centers(2:end),pdf_sw_dep_drug,'Color',colormap(1,:),'LineWidth',1);
histogram('BinEdges',vinterval,'BinCounts',sw_dep_histw_saline,'LineStyle','none','FaceColor',[colormap(1,:)+.1])
plot(centers(2:end),pdf_sw_dep_saline,'Color',[colormap(1,:)+.1],'LineWidth',1);
hold off
xlabel('Response Time'); ylabel('Liklihood of Response');
xlim([0 20]); xticks([0 6 12 18]); ylim([0 0.25]);
set(gca, 'Color', 'none', 'box','off');
title('Pooled Switch Histograms')
text(0.5, 0.20, strcat('mean switch time DRUG - ' , num2str(mean([GroupAnSt.(groupFields{2}).mean_switch_dep{:}])), 'secs'), 'FontSize', 12);
text(0.5, 0.17, strcat('mean switch time SALINE - ' , num2str(mean([GroupAnSt.(groupFields{1}).mean_switch_dep{:}])), 'secs'), 'FontSize', 12);
hold off;

% % %
% % % % IF DATA IS NOT NORMALLY DISTRIBUTED, CAN TURN TO WILCOXON RANLSUM
% % [p_switch,h_switch,stats_switch] = ranksum([GroupAnSt.(groupFields{1}).median_switch_dep{:}], [GroupAnSt.(groupFields{2}).median_switch_dep{:}]);
% % [p_cv,h_cv,stats_cv] = ranksum(GroupAnSt.(groupFields{1}).stats(7,:), GroupAnSt.(groupFields{2}).stats(7,:));
% % fprintf('%s - Median Switch (seconds) SALINE: %0.3f(%.3f-%.3f) vs DRUG: %0.3f(%.3f-%.3f), zval=%2.3f, p = %2.5f \n',...
% %     comparison_name,...
% %     mean([GroupAnSt.(groupFields{1}).totSwitchTrDe{:}]), quantile([GroupAnSt.(groupFields{1}).totSwitchTrDe{:}],0.25),quantile([GroupAnSt.(groupFields{1}).totSwitchTrDe{:}], .75),...
% %     mean([GroupAnSt.(groupFields{2}).totSwitchTrDe{:}]), quantile([GroupAnSt.(groupFields{2}).totSwitchTrDe{:}],0.25),quantile([GroupAnSt.(groupFields{2}).totSwitchTrDe{:}], .75),...
% %     stats_switch.zval, p_switch);
% % %


x = [1:2];
figure();
subplot(1,3,2);
hold on;
for i_An = 1:length(GroupAnSt.(groupFields{1}).stats(1,:))
    plot(x, [GroupAnSt.(groupFields{1}).stats(1,i_An) GroupAnSt.(groupFields{2}).stats(1,i_An)],...
        '-o', 'Color', colormap(1,:), 'MarkerEdgeColor', colormap(1,:), 'MarkerFaceColor', colormap(1,:))
end
errorbar(x(1)-0.15, mean(cellfun(@mean, GroupAnSt.Saline.totSwitchTrDe)), std(cellfun(@mean, GroupAnSt.Saline.totSwitchTrDe)),...
    '-s', 'Color', colormap(1,:), 'MarkerFaceColor', colormap(1,:), 'MarkerSize', 6, 'LineWidth', 3, 'CapSize', 10)
errorbar(x(2)+0.15, mean(cellfun(@mean,GroupAnSt.(groupFields{2}).totSwitchTrDe)), std(cellfun(@mean,GroupAnSt.(groupFields{2}).totSwitchTrDe)),...
    '-s', 'Color', colormap(2,:), 'MarkerFaceColor', colormap(1,:), 'MarkerSize', 6, 'LineWidth', 3, 'CapSize', 10)
set(gcf, 'Position',  [100, 50, 1200, 500]); xlim([.75 2.25]);
ylim([6 12]); yticks([6:1:12]); ylabel('Average Switch Time (s)');
set(gca,'XTick',[1 2]); xticklabels(["SALINE"; "DRUG"]);
hold off;


% AUSTIN ADDED PAIRED LINES PLOT FOR CV 3-6-2022  -- REMOVE FOR PRINT

subplot(1,3,3);
hold on;
for i_An = 1:length(GroupAnSt.(groupFields{1}).stats(1,:))
    plot(x, [GroupAnSt.(groupFields{1}).stats(7,i_An) GroupAnSt.(groupFields{2}).stats(7,i_An)],...
        '-o', 'Color', colormap(1,:), 'MarkerEdgeColor', colormap(1,:), 'MarkerFaceColor', colormap(1,:))
end
errorbar(x(1)-0.15, mean(GroupAnSt.(groupFields{1}).stats(7,:)), std(GroupAnSt.(groupFields{1}).stats(7,:)),...
    '-s', 'Color', colormap(1,:), 'MarkerFaceColor', colormap(1,:), 'MarkerSize', 6, 'LineWidth', 3, 'CapSize', 10)
errorbar(x(2)+0.15, mean(GroupAnSt.(groupFields{2}).stats(7,:)), std(GroupAnSt.(groupFields{1}).stats(7,:)),...
    '-s', 'Color', colormap(2,:), 'MarkerFaceColor', colormap(1,:), 'MarkerSize', 6, 'LineWidth', 3, 'CapSize', 10)
set(gcf, 'Position',  [100, 50, 1200, 500]); xlim([.75 2.25]);
ylim([10 50]); yticks([10:5:50]); ylabel('Coefficient of Variation'); ;
set(gca,'XTick',[1 2]); xticklabels(["SALINE"; "DRUG"]);
hold off;

%TABLES FOR MIXED LOGISTIC MODEL BUILT HERE

swTABLE_SALINE = [];
swTABLE_DRG = [];
counter = 1;
T = [];


for i_An = 1:length(GroupAnSt.(groupFields{1}).totSwitchTrDe)
    
    SW_TIME = GroupAnSt.(groupFields{1}).totSwitchTrDe{i_An}';
    DRUG = fieldnames(GroupAnSt); DRUG = repmat(string(DRUG{1}), length(SW_TIME), 1);
    ANIMAL_ID = repmat(string(uni_ids{i_An}), length(SW_TIME), 1);
    
    swTABLE_SALINE = table(SW_TIME, DRUG, ANIMAL_ID);
    
    SW_TIME = GroupAnSt.(groupFields{2}).totSwitchTrDe{i_An}';
    DRUG = fieldnames(GroupAnSt); DRUG = repmat(string(DRUG{2}), length(SW_TIME), 1);
    ANIMAL_ID = repmat(string(uni_ids{i_An}), length(SW_TIME), 1);
    
    swTABLE_DRUG = table(SW_TIME, DRUG, ANIMAL_ID);
    
    T = [T; swTABLE_SALINE; swTABLE_DRUG];
    
    
end

%PRINTING OUT FULL MODEL STATISTICS FOR NOW
behavior_model = fitglme(T, 'SW_TIME~DRUG+(1|ANIMAL_ID)')
% creates a mixed effects linear regression
anova(behavior_model)

%PAIRED WILCOXON TESTS FOR MOVEMENT RELATED METRICS
%CAN CREATE FIGURES FOR INDIVIDUAL METRICS IF DESIRED.
%UNLIKE SWITCH, THESE ARE USUALLY GAMMA DISTRIBUTED FOR EVERY ANIMAL

%GAMMA DIST

%%
% Matt: changed ranksum to signrank for 
% Matt: added in 'method', 'approximate' to each signrank test - function did not work when I took two mice out of the SCH group
%%


[p_srsp,h_srsp,stats_srsp] = signrank([GroupAnSt.(groupFields{1}).median_short_rsp{:}], [GroupAnSt.(groupFields{2}).median_short_rsp{:}], 'method', 'approximate');
cles_srsp = CLES([GroupAnSt.(groupFields{1}).median_short_rsp{:}; GroupAnSt.(groupFields{2}).median_short_rsp{:}]');

fprintf('\n %s - Median Short Response Time (seconds) SALINE: %0.3f(%.3f-%.3f) vs DRUG: %0.3f(%.3f-%.3f), zval=%2.3f, p = %2.5f, CLES = %2.3f \n',...
    comparison_name,...
    mean([GroupAnSt.(groupFields{1}).median_short_rsp{:}]), quantile([GroupAnSt.(groupFields{1}).median_short_rsp{:}],0.25),quantile([GroupAnSt.(groupFields{1}).median_short_rsp{:}], .75),...
    mean([GroupAnSt.(groupFields{2}).median_short_rsp{:}]), quantile([GroupAnSt.(groupFields{2}).median_short_rsp{:}],0.25),quantile([GroupAnSt.(groupFields{2}).median_short_rsp{:}], .75),...
    stats_srsp.zval, p_srsp, cles_srsp);

%GAMMA DIST

[p_lrsp,h_lrsp,stats_lrsp] = signrank([GroupAnSt.(groupFields{1}).median_long_rsp{:}], [GroupAnSt.(groupFields{2}).median_long_rsp{:}], 'method', 'approximate');
cles_lrsp = CLES([GroupAnSt.(groupFields{1}).median_long_rsp{:}; GroupAnSt.(groupFields{2}).median_long_rsp{:}]');


fprintf('\n %s - Median Long Response Time (seconds) SALINE: %0.3f(%.3f-%.3f) vs DRUG: %0.3f(%.3f-%.3f), zval=%2.3f, p = %2.5f, CLES = %2.3f \n',...
    comparison_name,...
    mean([GroupAnSt.(groupFields{1}).median_long_rsp{:}]), quantile([GroupAnSt.(groupFields{1}).median_long_rsp{:}],0.25),quantile([GroupAnSt.(groupFields{1}).median_short_rsp{:}], .75),...
    mean([GroupAnSt.(groupFields{2}).median_long_rsp{:}]), quantile([GroupAnSt.(groupFields{2}).median_long_rsp{:}],0.25),quantile([GroupAnSt.(groupFields{2}).median_short_rsp{:}], .75),...
    stats_lrsp.zval, p_lrsp, cles_lrsp);

[p_rt,h_rt,stats_rt] = signrank([GroupAnSt.(groupFields{1}).median_rt{:}], [GroupAnSt.(groupFields{2}).median_rt{:}], 'method', 'approximate');
cles_rt = CLES([GroupAnSt.(groupFields{1}).median_rt{:}; GroupAnSt.(groupFields{2}).median_rt{:}]');

%TT + RT 

fprintf('\n %s - Median RT (seconds) SALINE: %0.3f(%.3f-%.3f) vs DRUG: %0.3f(%.3f-%.3f), zval=%2.3f, p = %2.5f, CLES = %2.3f \n',...
    comparison_name,...
    mean([GroupAnSt.(groupFields{1}).median_rt{:}]), quantile([GroupAnSt.(groupFields{1}).median_rt{:}],0.25),quantile([GroupAnSt.(groupFields{1}).median_rt{:}], .75),...
    mean([GroupAnSt.(groupFields{2}).median_rt{:}]), quantile([GroupAnSt.(groupFields{2}).median_rt{:}],0.25),quantile([GroupAnSt.(groupFields{2}).median_rt{:}], .75),...
    stats_rt.zval, p_rt, cles_rt);

[p_tt,h_tt,stats_tt] = signrank([GroupAnSt.(groupFields{1}).median_tt{:}], [GroupAnSt.(groupFields{2}).median_tt{:}], 'method', 'approximate');
cles_tt = CLES([GroupAnSt.(groupFields{1}).median_tt{:}; GroupAnSt.(groupFields{2}).median_tt{:}]');

fprintf('\n %s - Median TT (seconds) SALINE: %0.3f(%.3f-%.3f) vs DRUG: %0.3f(%.3f-%.3f), zval=%2.3f, p = %2.5f, CLES = %2.3f \n',...
    comparison_name,...
    mean([GroupAnSt.(groupFields{1}).median_tt{:}]), quantile([GroupAnSt.(groupFields{1}).median_tt{:}],0.25),quantile([GroupAnSt.(groupFields{1}).median_tt{:}], .75),...
    mean([GroupAnSt.(groupFields{2}).median_tt{:}]), quantile([GroupAnSt.(groupFields{2}).median_tt{:}],0.25),quantile([GroupAnSt.(groupFields{2}).median_tt{:}], .75),...
    stats_tt.zval, p_tt, cles_tt);

[p_srd,h_srd,stats_srd] = signrank([GroupAnSt.(groupFields{1}).median_srd{:}], [GroupAnSt.(groupFields{2}).median_srd{:}], 'method', 'approximate');
cles_srd = CLES([GroupAnSt.(groupFields{1}).median_srd{:}; GroupAnSt.(groupFields{2}).median_srd{:}]');

% GAMMA DIST

fprintf('\n %s - Median Short Response Duration (seconds) SALINE: %0.3f(%.3f-%.3f) vs DRUG: %0.3f(%.3f-%.3f), zval=%2.3f, p = %2.5f, CLES=%2.3f \n',...
    comparison_name,...
    mean([GroupAnSt.(groupFields{1}).median_srd{:}]), quantile([GroupAnSt.(groupFields{1}).median_srd{:}],0.25),quantile([GroupAnSt.(groupFields{1}).median_srd{:}], .75),...
    mean([GroupAnSt.(groupFields{2}).median_srd{:}]), quantile([GroupAnSt.(groupFields{2}).median_srd{:}],0.25),quantile([GroupAnSt.(groupFields{2}).median_srd{:}], .75),...
    stats_srd.zval, p_srd, cles_srd);

%GAMMA DIST

[p_lrd,hlrd,stats_lrd] = signrank([GroupAnSt.(groupFields{1}).median_lrd{:}], [GroupAnSt.(groupFields{2}).median_lrd{:}], 'method', 'approximate');
cles_lrd = CLES([GroupAnSt.(groupFields{1}).median_lrd{:}; GroupAnSt.(groupFields{2}).median_lrd{:}]');

fprintf('\n %s - Median Long Response Duration (seconds) SALINE: %0.3f(%.3f-%.3f) vs DRUG: %0.3f(%.3f-%.3f), zval=%2.3f, p = %2.5f, CLES = %2.3f \n',...
    comparison_name,...
    mean([GroupAnSt.(groupFields{1}).median_lrd{:}]), quantile([GroupAnSt.(groupFields{1}).median_lrd{:}],0.25),quantile([GroupAnSt.(groupFields{1}).median_lrd{:}], .75),...
    mean([GroupAnSt.(groupFields{2}).median_lrd{:}]), quantile([GroupAnSt.(groupFields{2}).median_lrd{:}],0.25),quantile([GroupAnSt.(groupFields{2}).median_lrd{:}], .75),...
    stats_lrd.zval, p_lrd, cles_lrd);

%these were the only metrics consistently different, so here are response histograms

[short_rsps_histw_saline, vinterval,wt_ex_short_saline] = histwcc([GroupAnSt.(groupFields{1}).pooled_short_rsp], 0:0.40:24);
[short_rsps_histw_drug, vinterval,wt_ex_short_drug] = histwcc([GroupAnSt.(groupFields{2}).pooled_short_rsp], 0:0.40:24);

[long_rsps_histw_saline, vinterval,wt_ex_long_saline] = histwcc([GroupAnSt.(groupFields{1}).pooled_long_rsp], 0:0.40:24);
[long_rsps_histw_drug, vinterval,wt_ex_long_drug] = histwcc([GroupAnSt.(groupFields{2}).pooled_long_rsp], 0:0.40:24);

[pdf_shortrsps_saline,~] = ksdensity(cell2mat([GroupAnSt.(groupFields{1}).pooled_short_rsp]),centers(2:end),'Bandwidth',0.6,'Weights',wt_ex_short_saline);
[pdf_longrsps_saline,~] = ksdensity(cell2mat([GroupAnSt.(groupFields{1}).pooled_long_rsp]),centers(2:end),'Bandwidth',0.6,'Weights',wt_ex_long_saline);

[pdf_shortrsps_drug,~] = ksdensity(cell2mat([GroupAnSt.(groupFields{2}).pooled_short_rsp]),centers(2:end),'Bandwidth',0.6,'Weights',wt_ex_short_drug);
[pdf_longrsps_drug,~] = ksdensity(cell2mat([GroupAnSt.(groupFields{2}).pooled_long_rsp]),centers(2:end),'Bandwidth',0.6,'Weights',wt_ex_long_drug);

figure();
hold on;
histogram('BinEdges',vinterval,'BinCounts',short_rsps_histw_saline,'LineStyle','none','FaceColor',colormap(1,:))
plot(centers(2:end),pdf_shortrsps_saline,'Color',colormap(1,:),'LineWidth',1);
histogram('BinEdges',vinterval,'BinCounts',long_rsps_histw_saline,'LineStyle','none','FaceColor',[colormap(1,:)+.15])
plot(centers(2:end),pdf_longrsps_saline,'Color',[colormap(1,:)+.15],'LineWidth',1);
hold off
xlabel('Response Time'); ylabel('Liklihood of Response');
xlim([0 20]); xticks([0 6 12 18]); ylim([0 0.25]);
set(gca, 'Color', 'none', 'box','off');
title('SALINE Pooled Rsp Histograms')
text(0.5, 0.20, strcat('mdn short rsp ' , num2str(mean([GroupAnSt.(groupFields{1}).median_short_rsp{:}])), 'secs'), 'FontSize', 12);
text(0.5, 0.17, strcat('mdn long rsp ' , num2str(mean([GroupAnSt.(groupFields{1}).median_long_rsp{:}])), 'secs'), 'FontSize', 12);
hold off;

figure();
hold on;
histogram('BinEdges',vinterval,'BinCounts',short_rsps_histw_drug,'LineStyle','none','FaceColor',colormap(2,:))
plot(centers(2:end),pdf_shortrsps_drug,'Color',colormap(2,:),'LineWidth',1);
histogram('BinEdges',vinterval,'BinCounts',long_rsps_histw_drug,'LineStyle','none','FaceColor',[colormap(2,:)+.1])
plot(centers(2:end),pdf_longrsps_drug,'Color',[colormap(2,:)+.1],'LineWidth',1);
hold off
xlabel('Response Time'); ylabel('Liklihood of Response');
xlim([0 20]); xticks([0 6 12 18]); ylim([0 0.25]);
set(gca, 'Color', 'none', 'box','off');
title('DRUG Pooled Rsp Histograms')
text(0.5, 0.20, strcat('mdn short rsp ' , num2str(mean([GroupAnSt.(groupFields{2}).median_short_rsp{:}])), 'secs'), 'FontSize', 12);
text(0.5, 0.17, strcat('mdn long rsp ' , num2str(mean([GroupAnSt.(groupFields{2}).median_long_rsp{:}])), 'secs'), 'FontSize', 12);

hold off;

end

function dateNum = date2TableNum(animal, date,anSt)
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