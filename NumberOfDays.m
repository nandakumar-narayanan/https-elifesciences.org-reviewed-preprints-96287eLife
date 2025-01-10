groupNames = fieldnames(group);
animalID=[]; 
for i = 1:length(groupNames)
    counter = 1;
    for k = 1:length(group.(groupNames{i}))
        animalID{counter} = cell2mat(group.(groupNames{i}){1,k}(2)); counter=counter+1;        
    end
end

 animals = unique(animalID);
 for i = 1:length(animals)
        fileNum(i) = sum(cell2mat([strfind(animalID, animals{i})])); 
 end

 fprintf('# Files %0.2f (%0.2f-%0.2f)', median(fileNum), quantile(fileNum,0.25),quantile(fileNum,0.75))

 
