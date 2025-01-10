function [r2, accuracy,predY, C, r2shuffled, accuracyShuffled] =  classifyTrials(data, class)

layer_size = [35 20 10]; %[10]  % set up layers
dataShuffled = data(:, randperm(size(data,2))); % make shuffled data

clear predY  % clear the prediction variable
fprintf('  Classifying:')
N = size(data,2); alltrials = 1:N;
% Loop through each data point
for i_trial = 1:N
    fprintf('.')
    % Take ith point as test set
    testX = data(:,i_trial);
    testXshuffled = dataShuffled(:,i_trial); 
    
    % Take remaining points as training set
   trainX = data(:,find(1:N~=i_trial));
   trainXshuffed = dataShuffled(:,find(1:N~=i_trial));

   try;
    % Train a model on training set
    nbModel = fitcnet(trainX',class(find(1:N~=i_trial)),"LayerSizes",layer_size);
    %nbModel = fitcnb(trainX',class(find(1:N~=i_trial)));  tried NBM, not
    %as good

    
    % Make prediction on test point        
    predY(i_trial) = predict(nbModel,testX');

   catch
    predY(i_trial) = -1;  %sometimes the model doesn't converge
   end
    
   try
    nbModelShuffled = fitcnet(trainXshuffed',class(find(1:N~=i_trial)),"LayerSizes",layer_size);
    predShuffled(i_trial) = predict(nbModelShuffled,testXshuffled');
   catch
    predShuffled(i_trial) = -1;
   end       
    
end


C = confusionmat(predY, class);
r = corrcoef(predY, class); r2= r(1,2)^2; 
% Calculate average 
% Get total number of samples
total = sum(sum(C));

% Get number of correct predictions on diagonal
correct = sum(diag(C));

% Calculate accuracy
accuracy = correct/total;

Cshuffled = confusionmat(predShuffled, class); correctShuffled = sum(diag(Cshuffled)); accuracyShuffled = correctShuffled/total;
r = corrcoef(predShuffled, class); r2shuffled= r(1,2)^2; 
fprintf('%0.2f, %0.2f', accuracy, accuracyShuffled) % report results



