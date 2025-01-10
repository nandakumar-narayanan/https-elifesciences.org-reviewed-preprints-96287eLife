function stats = EnsembleAnalysisSwitch(neurons, kernel_bw, bin_size)

for i = 1:length(neurons); num_reward_trial(i) = length(neurons(i).events.switchTrialInit); end;
 [N,edges] = histcounts(num_reward_trial,1:max(num_reward_trial),'Normalization','cdf');
 max_trial_num = edges(find(N>0.1,1)); % >90% neuron number

 max_trial_num = 20; % 14 gets every trial; but 20 is a more reasonable # for classification.
interval=18; 
 pre_post = 6;
llookback = pre_post;
llookforward = interval + pre_post;  
intervals = -llookback:bin_size:llookforward;
time = intervals(1:end-1);

% fitcnb treat NaN as missing data, so we can put neurons with smaller trial num
data = NaN(numel(time),max_trial_num,length(neurons));



for i_neuron = 1:length(neurons)
   
    trialStart = neurons(i_neuron).events.switchTrialInit;   
    a = peSpike(neurons(i_neuron).spikeTS, trialStart, intervals);
    a = a(~all(isnan(a),2),:); % remove trial with no spike
    a = a(1:min([max_trial_num size(a,1)]),:); % <= max trial  
    b = zeros(numel(time),size(a,1));
    for i_trial = 1:size(a,1)
        b(:,i_trial) = ksdensity(a(i_trial,:),time,'Bandwidth',kernel_bw,'Support',[intervals(1)-0.1 intervals(end)+0.1],'BoundaryCorrection','reflection');
    end
    data(:,1:size(b,2),i_neuron) = b;                                                                                                     
end
data = zscore(data);
fprintf('time %d x trial %d x neuron %d, blank: %d\n',size(data), length(find(isnan(data)))) % time x trial x neuron

numtrial=size(data,2); 
bayesdata = reshape(data,[],size(data,3));
t_time = repmat(time,1,numtrial);
t_trial_num = reshape(repmat(1:numtrial,numel(time),1),1,[]);
shuff_time = cell2mat(cellfun(@(x) x(randperm(numel(x))),repmat({time},1,numtrial),'UniformOutput',false));
goodtimes = time>=0 & time<=interval;  

earlytimes = time>=0 & time<=6;  % Reviewer wanted 0 - 6 seconds
midtimes = time>=6 & time<=12;  % Reviewer wanted 6 - 18 seconds
latetimes = time>=12 & time<=18;  % Reviewer wanted 6 - 18 seconds


confusion_mat = zeros(numel(time),numel(time));
confusion_mat_shuff = zeros(numel(time),numel(time));
error = zeros(sum(goodtimes),1);
shufferror = zeros(sum(goodtimes),1);
t = tic;
% while we can do cross validation with 'Leaveout' in the fitcnb, parfor with manual cross validation is much faster.
parfor testtrial = 1:numtrial
    test_trial_idx = t_trial_num == testtrial;
    test_data= bayesdata(test_trial_idx,:);
    train_pred = bayesdata(~test_trial_idx,:);
    train_time = t_time(~test_trial_idx);
    
    NBModel = fitcnb(train_pred,train_time); %,'Distribution','kernel','KSWidth',3 make the model including trainingdata and class labels (i.e., what time in corresponding row was)
    ShuffModel = fitcnb(train_pred,shuff_time(~test_trial_idx));
    
    predictLabels1 = predict(NBModel,test_data);
    predictLabels2 = predict(ShuffModel,test_data);

    rsquared = corrcoef(predictLabels1(goodtimes),time(goodtimes));
    r2(testtrial) = rsquared(1,2)^2;  
    shuffsquared = corrcoef(predictLabels2(goodtimes),time(goodtimes)); 
    shuffr2(testtrial) = shuffsquared(1,2)^2;  
    
    rsquared_0_6 = corrcoef(predictLabels1(earlytimes),time(earlytimes));
    r2_0_6(testtrial) = rsquared_0_6(1,2)^2;  
    
    rsquared_6_12 = corrcoef(predictLabels1(midtimes),time(midtimes));
    r2_6_12(testtrial) = rsquared_6_12(1,2)^2;  
    
    rsquared_12_18 = corrcoef(predictLabels1(latetimes),time(latetimes));
    r2_12_18(testtrial) = rsquared_12_18(1,2)^2;  
        
    confusion_mat = confusion_mat + confusionmat(time',predictLabels1);
    confusion_mat_shuff = confusion_mat_shuff + confusionmat(time',predictLabels2);
    % compute error
    error = error + abs(predictLabels1(goodtimes)-time(goodtimes)');
    shufferror = shufferror + abs(predictLabels2(goodtimes)-time(goodtimes)');
end

confusion_mat = confusion_mat/numtrial;
confusion_mat_shuff = confusion_mat_shuff/numtrial;
error = error/numtrial;
shufferror = shufferror/numtrial;

smoothf = 1.6; % in sec
for count=1:numel(time)
    y_s = gausssmooth(confusion_mat(count,:),smoothf/bin_size);
    y(count,:) = y_s/max(y_s);
    y_s = gausssmooth(confusion_mat_shuff(count,:),smoothf/bin_size);
    y2(count,:) = y_s/max(y_s);
end

trial_num_neuron = squeeze(sum(~isnan(data(1,:,:))));
missing_trial_ratio = sum(max(trial_num_neuron) - trial_num_neuron)/(max(trial_num_neuron)*size(data,3));

stats.r2 = r2;   
stats.r2_0_6 = r2_0_6;  
stats.r2_6_12 = r2_6_12;  
stats.r2_12_18 = r2_12_18;  
stats.shuffr2 = shuffr2;
stats.y = y;
stats.y2 = y2;
stats.confusion = confusion_mat; 
stats.confusionShuff = confusion_mat_shuff;
stats.data = data; 
stats.times = time; 
stats.goodtimes = goodtimes;

stats.trialerror = error; 
stats.shufftrialerror = shufferror; 
stats.missing_trial_ratio = missing_trial_ratio;

