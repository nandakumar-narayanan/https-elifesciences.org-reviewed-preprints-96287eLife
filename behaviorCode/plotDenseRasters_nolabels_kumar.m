function SessionAnSt = plotSession(TrialAnSt, savefolder_name, ind, no_plot)
if nargin < 4
    no_plot = false;
end

animalNames = fieldnames(TrialAnSt);
if nargout > 1
    SessionAnSt = struct;
    for i_animal = 1:length(animalNames)
        animal = animalNames{i_animal};
        trialWithDataIdx = find(~cellfun('isempty', {TrialAnSt.(animal)}));
        numOfDay = numel(trialWithDataIdx);
        for i_day = 1:numOfDay
            tmp_struct = TrialAnSt(trialWithDataIdx(i_day)).(animal);
            SwitchTrDe = [tmp_struct.SwitchDepart]; SwitchTrAr = [tmp_struct.SwitchArrival];
            ShortRsp = [tmp_struct.ShortRsp]; LongRsp = [tmp_struct.LongRsp];
            % SessionAnSt.(animal).Stats = tSummaryStatsTable;
            ssAnSt.SwitchTrDe = SwitchTrDe;
            ssAnSt.SwitchTrAr = SwitchTrAr;
            ssAnSt.ShortRsp = ShortRsp;
            ssAnSt.LongRsp = LongRsp; 
            SessionAnSt(trialWithDataIdx(i_day)).(animal) = ssAnSt;   
        end
    end
end
% 
if no_plot
    return
end
%% Plots each session for each animal (across days)

for i_animal = 1:length(animalNames)
    animal = animalNames{i_animal};
%     sessionplot(TrialAnSt, animal); 
%     f1_FileName = sprintf('CDF_%s.png', animal); f1_fullFileName = fullfile(savefolder_name, f1_FileName); saveas(gcf, f1_fullFileName);
%     close(gcf)    
%     figure('Units', 'Normalized', 'OuterPosition', [0 0.04 1, 0.96])

    rasterplot(TrialAnSt, animal, ind);
%     f1_FileName = sprintf('raster_%s.png', animal); f1_fullFileName = fullfile(savefolder_name, f1_FileName); saveas(gcf, f1_fullFileName);
%     close(gcf);
end

end

function sessionplot(TrialAnSt, animal)
centers = 0:0.1:24;
cmap = lines(2);

trialWithDataIdx = find(~cellfun('isempty', {TrialAnSt.(animal)}));
numOfDay = numel(trialWithDataIdx);
figure('Units', 'Normalized', 'OuterPosition', [0 0.04 1, 0.96])
subplot_idx = 1;
for i_day = 1:numOfDay
    
    tmp_struct = TrialAnSt(trialWithDataIdx(i_day)).(animal);
    numplot_y = ceil(sqrt(numOfDay)); %  positions and sizes for subplot
    numplot_x = ceil(numOfDay/numplot_y)*2;

    % Pulls out the long trials that should contain a switch
    SwitchTrDe = [tmp_struct.SwitchDepart]; SwitchTrAr = [tmp_struct.SwitchArrival];
    ShortRsp = [tmp_struct.ShortRsp]; LongRsp = [tmp_struct.LongRsp];

    % Calculates and builds summary stats table
    [Avg, STD, N, SEM, Median, Q1, Q2, IntQRange, Limen, WR] = getStat(SwitchTrDe);
    % subplots the daily CDFs
    if ~isempty(SwitchTrDe) && ~isempty(SwitchTrAr)
        subplot(numplot_y, numplot_x, subplot_idx);
        hold on; ecdf(SwitchTrDe); ecdf(SwitchTrAr); 
        yticks(0:0.25:1); ylabel({'Cumulative Probability of a Switch'});
        xlim([0 20]); xticks([0 6 12 18]); xlabel({'Trial Duration (s)'});
        txt1 = ['N=' num2str(N)]; txt2 = ['Median=' num2str(Median)]; txt3 = ['Avg=' num2str(Avg) ' +/- ' num2str(SEM)]; txt4 = ['Int. Quart. Range=' num2str(IntQRange)];
        text(1, 0.9, txt1); text(1, 0.8, txt2); text(1, 0.7, txt3); text(1, 0.6, txt4);
        hold off;
        grid on;
        set(gca, 'Color', 'none', 'box','off'); 
       % title(sprintf('%s %s',animal, tmp_struct(1).mpc.StartDate));

        [pdf_ShortRsp,~] = ksdensity(ShortRsp,centers(2:end),'Bandwidth',0.6);
        [pdf_LongRsp,~] = ksdensity(LongRsp,centers(2:end),'Bandwidth',0.6);
        subplot(numplot_y, numplot_x, subplot_idx+1);
        hold on
        histogram(ShortRsp,0:0.5:24,'Normalization','probability','LineStyle','none','FaceColor',cmap(1,:))
        plot(centers(2:end),pdf_ShortRsp,'Color',cmap(1,:),'LineWidth',2)
        histogram(LongRsp,0:0.5:24,'Normalization','probability','LineStyle','none','FaceColor',cmap(2,:))
        plot(centers(2:end),pdf_LongRsp,'Color',cmap(2,:),'LineWidth',2)
        hold off
        xlabel('Trial Duration(s)'); ylabel('PDF of Short and Long Rsp');
        xlim([0 24]); xticks([0 6 12 18]);
        set(gca, 'Color', 'none', 'box','off'); 
        title(sprintf('%s %s',animal, tmp_struct(1).mpc.StartDate));

        subplot_idx = subplot_idx+2;
    end
end


end

function rasterplot(TrialAnSt, animal, ind)
    trialWithDataIdx = find(~cellfun('isempty', {TrialAnSt.(animal)}));
    numOfDay = numel(trialWithDataIdx);
%     hold on; 
%     subplot(6, 4, ind)
    subplot(6, 2, ind);
    subplot_idx = [1 2];
    numplot_y = ceil(sqrt(numOfDay));
    numplot_x = ceil(numOfDay/numplot_y);
    for i_day = 1:numOfDay
        
        trialData = TrialAnSt(trialWithDataIdx(i_day)).(animal);
        
        
        fields = {'ShortRsp' 'LongRsp' 'SwitchDepart' 'SwitchArrival'};
%         fields = {'leftRspTimeTrial' 'rightRspTimeTrial' 'SwitchDepart' 'SwitchArrival'};
  %      colormap = [1 .5 0.5; .5 .5 1; 1 0 0; 0 0 1];
        
         colormap = [0.7 0.7 0.7; 0 0 0; .2 .8 .2; 0 0.5 1];
         %    {'ShortRsp'}    {'LongRsp'}    {'SwitchDepart'}    {'SwitchArrival'}
        lwidth = [2 2 3 3];
        for i_trial_type = 1:1
            if i_trial_type == 1
                trialTypeIdx = [trialData.programmedDuration] == 18000;
%                 subplot(6, 4, ind)
%             elseif i_trial_type == 2
%                 trialTypeIdx = [trialData.programmedDuration] == 6000;
%                 subplot(6, 4, ind + 1);
            end
            
            hold on
            for i_rsp_type = 1:4
                [rx, ry] = get_field_raster(trialData, fields{i_rsp_type}, trialTypeIdx);
                pH = plot(rx, ry,'Color',colormap(i_rsp_type,:),'LineWidth',lwidth(i_rsp_type)); 
               
            end
%             plot([18 18],[0 ceil(max(max(ry)))],'LineWidth',1,'Color',[0.8 0.8 0.8]);
            plot([6 6],[0 ceil(max(max(ry)))],'LineWidth',1,'Color',[0.8 0.8 0.8]);
%             hold off
            
            set(gca, 'Color', 'none', 'box','off'); 
%             xlabel('Trial Time (seconds)'); ylabel('Trial Number');
            if i_trial_type == 1
                xlim([0 18]);
                axis xy; axis off; 
           %      title(sprintf('%s -%s Long Trial',animal,trialData(1).mpc.StartDate));
            elseif i_trial_type == 2
                xlim([0 10]);
                 axis xy; axis off;  
%                 title('Short Trial');
            end            
        end
        subplot_idx = subplot_idx+3;
        
    end
%     hold on

end

% builds summary stats table
function [tAvg, tSTD, tN, tSEM, tMedian, tQ1, tQ2, tIntQRange, tLimen, tWR] = getStat(x)
    tAvg = mean(x); tSTD = std(x); tN = size(x); tN(1) = []; tSEM = tSTD/(sqrt(tN));
    tMedian = quantile(x, 0.5); tQ1 = quantile(x, 0.25); tQ2 = quantile(x, 0.75); 
    tIntQRange = tQ2 - tQ1; tLimen = tIntQRange/2; tWR = tLimen/tMedian;
    % bstats = [tAvg, tSTD, tN, tSEM, tMedian, tQ1, tQ2, tIntQRange, tLimen, tWR];
end


function [rx, ry] = get_field_raster(trialData, fieldname, trialTypeIdx)
    RPByTrialNum = takefield(trialData, fieldname);
    RPByTrialNum(RPByTrialNum==0) = NaN;
    RPByTrialNum = RPByTrialNum(:,trialTypeIdx);
    [rx, ry] = rasterize(RPByTrialNum,1); 
end

function [x1, y1] = rasterize(data,trialStartNum,lineLength)

if nargin < 3
    lineLength = 1;
end

trialgap = 0.1;
yyStart = trialStartNum-1-trialgap;
numData = size(data,1);
x1 = NaN(3*size(data,1),size(data,2));
y1 = x1;
x1(1:3:3*numData,:) = data;
x1(2:3:3*numData,:) = data;
y1(1:3:3*numData,:) = repmat((1:size(data,2))+ trialStartNum + yyStart,size(data,1),1);
y1(2:3:3*numData,:) = y1(1:3:end,:)+lineLength-trialgap;

end

function stmat = takefield(st,fieldname)

numelst = cellfun(@numel,{st.(fieldname)});
stmat = zeros(max(numelst),length(st));
for i = 1:length(st)
    stmat(1:numelst(i),i) = st(i).(fieldname);

end

end

