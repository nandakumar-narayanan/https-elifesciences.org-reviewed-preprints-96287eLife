function [dataCondi, ci] = gksmooth(periEvtSPK, interval, FWHM)

% Estimate neuronal spike firing rate with gaussian kernel smoothing
% data = gksmooth(periEvtSPK, interval, FWHM)

% Inputs
%   - periEvtSPK: raw spike timestamps in seconds required 
%                 cell type for multiple condition - output from peSpike.
%   - interval: time vector of peri-event time, required (eg: interval = -1:0.1:+1)
%   - FWHM: optional, bandwidth of the kernel-smoothing window, default: 2*interval bin size (eg: FWHM = 0.1)

% Outputs
%   - data: estimated neuronal spike firing rate in Hz

% default bandwidth for kernel
if nargin < 3
    binsize = interval(2)-interval(1);
    FWHM = binsize*2;
end

% single condition
if ~iscell(periEvtSPK)
    periEvtSPK = {periEvtSPK};
end

% kernel density estimation
numCond = length(periEvtSPK);
dataCondi = cell(numCond,1);
for condIdx = 1:numCond
    if ~isempty(periEvtSPK{condIdx})
        periEvtSPKCond = periEvtSPK{condIdx};
        cellTime = periEvtSPKCond(~isnan(periEvtSPKCond));
        nelements = histc(cellTime, interval);
        f = ksdensity(cellTime, interval, 'Bandwidth', FWHM);
        data = f * sum(nelements) / size(periEvtSPKCond, 1);
        dataCondi{condIdx} = data;
    else
        dataCondi{condIdx} = [];
    end
end
if numCond == 1
    dataCondi = dataCondi{1};
end

% confidence intervals calculation by boot straping 
% number of bootstrap samples: neural data fit to be rms noise 0.2 (100 trials: 1500, 1000 trials 180)
if nargout > 1
    ci = cell(numCond,1);
    for condIdx = 1:numCond
        if ~isempty(periEvtSPK{condIdx})
            periEvtSPKCond = periEvtSPK{condIdx};
            numTrial = size(periEvtSPKCond,1);
            nbs = round(4.009e+04 .* numTrial.^ -0.7264);
            periEvtSPKCond = periEvtSPKCond(sum(~isnan(periEvtSPKCond),2)>0,:);
            if size(periEvtSPKCond,1) > 5
                nelements = histc(periEvtSPKCond(~isnan(periEvtSPKCond)), interval);
                defun = @(x) ksdensity(x(~isnan(x)),interval,'Bandwidth', FWHM); 
                opt = statset('UseParallel',false);
                B = bootci(nbs,{defun,periEvtSPKCond},'alpha',0.05,'type','bca','Options',opt);
                ci{condIdx} = B .* sum(nelements) ./ numTrial;
            else
                ci{condIdx} = [];
            end
        else
            ci{condIdx} = [];
        end
    end
    if numCond == 1
        ci = ci{1};
    end
end


end

