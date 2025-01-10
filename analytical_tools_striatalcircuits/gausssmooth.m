function [data] = gausssmooth(dataIn, FWHM); 

if nargin ==1; FWHM = 15; end; 

sig = FWHM/sqrt(8*log(2));

% % % Slow, but it's done right with a true gaussian kernel             
y = dataIn; sy = zeros(size(y));  x = 1:length(y); 
for xi = 1:length(y);
  kerny_i =  exp(-(x-xi).^2/(2*sig^2));
  kerny_i = kerny_i / sum(kerny_i);
  sy(xi) = sum(y.*kerny_i);
end

data= sy; 
