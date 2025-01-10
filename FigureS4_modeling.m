clear
load D1_D2_MSNs_x_SpkLogistic_convSpks_1000ms.mat
%%
mouse = 4; %iER4
trial = 3; % 3, 5, 6, 7,8, 10,11,14, 17, 18, 20

mouseName = D2_LogisticConv{1,mouse+1}; 
numNeurons = length(D2_LogisticConv{2,mouse+1});

xt = D2_LogisticConv{9,mouse+1};

trialInterval =  D2_LogisticConv{10,mouse+1}; % [0 18] in seconds
time = trialInterval(1):0.001:trialInterval(2);

switchTime = D2_LogisticConv{5,mouse+1};
predictedTime = D2_LogisticConv{11,mouse+1};

figure(42)
subplot(221);
plot(time,xt(trial,:),'r','Linewidth',2)
hold on
plot(time, 0.5*ones(size(time)),':k','Linewidth',1)

text(12,0.1,['PredictedTime: ' num2str(predictedTime(trial),'%1.2f') ' s'])
text(12,0.3,['SwitchTime: ' num2str(switchTime(trial),'%1.2f') ' s'])

xlabel('Time (s)')
ylabel('Ensemble activity x(t)')
title(['D2-mouse ' mouseName ...
        ' (N=' num2str(numNeurons) ' neurons) - trial ' num2str(trial) ...
        ' -- Network Activity']); 

%%
subplot(223);
mouse = 1; % dEL7
trial = 14; % 4, 8, 9,11, 14, 15, 17, 18, 20, 21, 22

mouseName = D1_LogisticConv{1,mouse+1}; 
numNeurons = length(D1_LogisticConv{2,mouse+1});

xt = D1_LogisticConv{9,mouse+1};

trialInterval =  D1_LogisticConv{10,mouse+1}; % [0 18] in seconds
time = trialInterval(1):0.001:trialInterval(2);

switchTime = D1_LogisticConv{5,mouse+1};
predictedTime = D1_LogisticConv{11,mouse+1};


plot(time,xt(trial,:),'b','Linewidth',2)
hold on
plot(time, 0.5*ones(size(time)),':k','Linewidth',1)

text(12,0.6,['PredictedTime: ' num2str(predictedTime(trial),'%1.2f') ' s'])
text(12,0.7,['SwitchTime: ' num2str(switchTime(trial),'%1.2f') ' s'])

xlabel('Time (s)')
ylabel('Ensemble activity x(t)')
title(['D1-mouse ' mouseName ...
        ' (N=' num2str(numNeurons) ' neurons) - trial ' num2str(trial) ...
        ' -- Network Activity']); 


load asims_D1_D2_MSNs_Accuracy_convSpks_1000ms.mat

d1accVals = round(cell2mat(D1_LogisticConv_SimsResults(2,2:6))); 
              % [100 41 46 0 100]; 
d2accVals = round(cell2mat(D2_LogisticConv_SimsResults(2,2:5))); 
              % [4 12 92 95]

d1NumberTrials = [15 6 6 1 13];
d2NumberTrials = [3 4 12 13];

%%
d1Poisson = [cell2mat(D1_LogisticConv_SimsResults(3,2))
    cell2mat(D1_LogisticConv_SimsResults(3,3))
    cell2mat(D1_LogisticConv_SimsResults(3,4))
    cell2mat(D1_LogisticConv_SimsResults(3,5))
    cell2mat(D1_LogisticConv_SimsResults(3,6))];

d2Poisson = [cell2mat(D2_LogisticConv_SimsResults(3,2))
    cell2mat(D2_LogisticConv_SimsResults(3,3))
    cell2mat(D2_LogisticConv_SimsResults(3,4))
    cell2mat(D2_LogisticConv_SimsResults(3,5))];

%%

subplot(2,2,[2 4])
scatter(0.15+ones(1,100)*d1NumberTrials(1),d1Poisson(1,:),[],'k','o','MarkerFaceColor','k')
hold on
scatter(-0.15+ones(1,100)*d1NumberTrials(2),d1Poisson(2,:),[],'k','o','MarkerFaceColor','k')
scatter(0.15+ones(1,100)*d1NumberTrials(3),d1Poisson(3,:),[],'k','o','MarkerFaceColor','k')
scatter(0.15+ones(1,100)*d1NumberTrials(4),d1Poisson(4,:),[],'k','o','MarkerFaceColor','k')
scatter(0.15+ones(1,100)*d1NumberTrials(5),d1Poisson(5,:),[],'k','o','MarkerFaceColor','k')

scatter(-0.15+ones(1,100)*d2NumberTrials(1),d2Poisson(1,:),[],[0.8 0.8 0.8],'o','MarkerFaceColor',[0.8 0.8 0.8])
hold on
scatter(-0.15+ones(1,100)*d2NumberTrials(2),d2Poisson(2,:),[],[0.8 0.8 0.8],'o','MarkerFaceColor',[0.8 0.8 0.8])
scatter(-0.15+ones(1,100)*d2NumberTrials(3),d2Poisson(3,:),[],[0.8 0.8 0.8],'o','MarkerFaceColor',[0.8 0.8 0.8])
scatter(-0.15+ones(1,100)*d2NumberTrials(4),d2Poisson(4,:),[],[0.8 0.8 0.8],'o','MarkerFaceColor',[0.8 0.8 0.8])

%%
D1ACC = [1 1 1 1 1; d1NumberTrials; d1accVals]; 
D2ACC = [2 2 2 2;d2NumberTrials; d2accVals];
scatter(D1ACC(2,:), D1ACC(3,:),100, 'bo', ...
    'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.5)
hold on
scatter(D2ACC(2,:), D2ACC(3,:),100, 'ro', ...
    'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.5)
xlabel('Ensemble Size (Neurons)'); 
ylabel('Logistic Regression Accuracy'); 

