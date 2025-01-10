function [histw, vinterval,wt_ex] = histwcc(vvc, vinterval)

ww = cellfun(@(x) size(x,2),vvc);
wt = 1./ww .* min(ww);
wt_ex = [];
for i_ani = 1:length(vvc)
    wt_ex  = [wt_ex; repmat(wt(i_ani),ww(i_ani),1)];
end
[histw, vinterval] = histwc(cell2mat(vvc), wt_ex', vinterval);
histw = histw ./ sum(histw);
histw = histw(1:end-1);
wt_ex = wt_ex';

end

function [histw, vinterval] = histwc(vv, ww, vinterval)
%   minV  = min(vv);
%   maxV  = max(vv);
%   delta = (maxV-minV)/nbins;
%   vinterval = linspace(minV, maxV, nbins)-delta/2.0;
  nbins = length(vinterval);
  histw = zeros(nbins, 1);
  for i=1:length(vv)
    ind = find(vinterval < vv(i), 1, 'last' );
    if ~isempty(ind)
      histw(ind) = histw(ind) + ww(i);
    end
  end
end