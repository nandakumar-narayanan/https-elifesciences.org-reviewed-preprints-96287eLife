

T = readtable('OptoBehavioralTableV3.csv');
D2creMice = (strcmp(T.GENOTYPE,'OptoD2Halo'));
T = T(D2creMice,:);

LaserStatus = 'ctrl';  
[dataLaserOff,m_LaserOff,md_LaserOff,cv_LaserOff] = AnimalData(T,LaserStatus);


groups = round(rand(1,length(dataLaserOff))*9)+1; % Sets up 10 groups for 10-fold CV
for i_CV = 1:10; 
    pdL1_Off = fitdist(dataLaserOff(groups~=i_CV),'Gamma');  % Finds groups that are not equal to each fold
    [dataCDF] = cdfcalc(dataLaserOff(groups==i_CV)); % And now finds groups that are equal to each fold
    xx1 = linspace(0,1,size(dataCDF,1));
    modelCDF = cdf(pdL1_Off,xx1);
    dataVsModel = fitlm(dataCDF, modelCDF);

    R(i_CV) = dataVsModel.Rsquared.Adjusted;
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------------------------------------------------------------%%
function [x,mx,mdx,cvx] = AnimalData(T,LaserStatus)
% 
x = nan(size(T,1),1);

%%%%% animal = unique(T.ANIMALID(strcmp(T.LASER,LaserStatus))) ;
selection = (strcmp(T.LASER,LaserStatus)) ;

if ~isempty(selection)
   x(1:sum(selection)) = T(selection,1).RSPTIME; 
end
x = x(~isnan(x));

mx=mean(x);
mdx=median(x);
cvx = std(x)./mean(x); 
%
end
%-------------------------------------------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%