
function EffectSize = cles(data)
%data is a 2d array with paired median values in columns 
mu_1 = mean(data(:,1)); mu_2 = mean(data(:,2));
var_1 = var(data(:,1)); var_2 = var(data(:,2));

EffectSize = normcdf((mu_1-mu_2)/sqrt(var_1+var_2));

end 