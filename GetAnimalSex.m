function Sex = GetAnimalSex(ID)

Females = {'PB1', 'PB2','PB5','PB6','IOL2','SL3', 'SL9','SR1','SR2','CL3', 'CR5','CR6', 'EP8', 'ER5', 'icOL1', 'icOL6', 'iEL3','iER1', 'dEL7', 'dEL8', 'iEL3', 'iER1', 'iER4', 'dcOR1', 'dcOR5', 'dcOR6'}; 
Males = {'PB3','PB7','PB8','PB11','PB12','IOL1','IOR3','SL10','SL6','SL8','CL6','CR2','CR3','icOR2','icOR3','icOR5','dcOL2','dcOL3','ER3','ER4','dER1','dER6','iEL5','iER2','iER4','dER1','dER2','dER9','iEL5'};

Sex=0; 
for i = 1:length(Females)
    if contains(ID, Females{i})
        Sex=1;
    end
end





