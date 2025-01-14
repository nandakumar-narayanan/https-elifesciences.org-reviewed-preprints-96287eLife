function [p, I] = saltP(spt_baseline,spt_test,dt,wn)
%SALT   Stimulus-associated spike latency test.
%   [P I] = SALT(SPT_BASELINE,SPT_TEST,DT,WN) calculates a modified version
%   of Jensen-Shannon divergence (see Endres and Schindelin, 2003) for
%   spike latency histograms.
%
%   Input arguments:
%       SPT_BASELINE - Discretized spike raster for stimulus-free baseline
%           period. N x M binary matrix with N rows for trials and M 
%           columns for spikes. Spike times have to be converted to a
%           binary matrix with a temporal resolution provided in DT. The
%           baseline segment has to excede the window size (WN) multiple
%           times, as the length of the baseline segment divided by the
%           window size determines the sample size of the null
%           distribution (see below).
%       SPT_TEST - Discretized spike raster for test period, i.e. after
%           stimulus. N x M binary matrix with N rows for trials and M 
%           columns for spikes. Spike times have to be converted to a
%           binary matrix with a temporal resolution provided in DT. The
%           test segment has to excede the window size (WN). Spikes out of
%           the window are disregarded.
%       DT - Time resolution of the discretized spike rasters in seconds.
%       WN - Window size for baseline and test windows in seconds
%           (optional; default, 0.001 s).
%
%   Output arguments:
%       P - Resulting P value for the Stimulus-Associated spike Latency
%           Test.
%       I - Test statistic, difference between within baseline and 
%           test-to-baseline information distance values. 
%
%   Briefly, the baseline binned spike raster (SPT_BASELINE) is cut to 
%   non-overlapping epochs (window size determined by WN) and spike latency
%   histograms for first spikes are computed within each epoch. A similar
%   histogram is constructed for the test epoch (SPT_TEST). Pairwise
%   information distance measures are calculated for the baseline
%   histograms to form a null-hypothesis distribution of distances. The
%   distances of the test histogram and all baseline histograms are
%   calculated and the median of these values is tested against the
%   null-hypothesis distribution, resulting in a p value (P).
%
%   Reference:
%   Endres DM, Schindelin JE (2003) A new metric for probability
%   distributions. IEEE Transactions on Information Theory 49:1858-1860.
%
%   See also JSDIV.
%
%   Please cite:
%
%   Kvitsiani D*, Ranade S*, Hangya B, Taniguchi H, Huang JZ, Kepecs A (2013)
%   Distinct behavioural and network correlates of two interneuron types in
%   prefrontal cortex. Nature 498:363?6.
%
%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   16-Dec-2012

% Input argument check
% error(nargchk(3,4,nargin))
% if nargin < 4
%     wn = 0.01;   % window size in s
% end
wn = wn * 1000;   % convert window size to ms
dt = dt * 1000;   % convert time resolution of bin raster to ms

% Latency histogram - baseline
[~, st] = size(spt_baseline);  % number of trials and number of baseline (pre-stim) data points
nmbn = round(wn/dt);   % number of bins for latency histograms
edges = 0:nmbn+1;  % bin boundaries
nm = floor(st/nmbn);   % size of the null dist. sample
nhlsi = zeros(nmbn+1,nm);    % preallocate latency histograms
for t = 1:nm   % loop through baseline windows
    [h,cmin2] = max(spt_baseline(:,(t-1)*nmbn+1:t*nmbn),[],2); % column wise find('first')
    nhlsi(:,t) = histcounts(cmin2.*h,edges,'Normalization','probability');
end

kn = nm + 1;   % number of all windows (nm baseline win. + 1 test win.)
% ISI histogram - test
[h,cmin2] = max(spt_test(:,1:nmbn),[],2);
nhlsi(:,kn) = histcounts(cmin2.*h,edges,'Normalization','probability');

% JS-divergence
jsd = nan(kn,kn);
for k1 = 1:kn
    D1 = nhlsi(:,k1);  % 1st latency histogram 
    for k2 = k1+1:kn
        D2 = nhlsi(:,k2);   % 2nd latency histogram
        jsd(k1,k2) = sqrt(JSdiv(D1,D2)*2);  % pairwise modified JS-divergence (real metric!)
    end
end

% Calculate p-value and information difference
[p, I] = makep(jsd,kn);

% -------------------------------------------------------------------------
function [p_value, Idiff] = makep(kld,kn)
% Calculates p value from distance matrix.

pnhk = kld(1:kn-1,1:kn-1);
nullhypkld = pnhk(~isnan(pnhk));   % nullhypothesis
testkld = median(kld(1:kn-1,kn));  % value to test
sno = length(nullhypkld(:));   % sample size for nullhyp. distribution
p_value = length(find(nullhypkld>=testkld)) / sno;
Idiff = testkld - median(nullhypkld);   % information difference between baseline and test latencies

% -------------------------------------------------------------------------
function D = JSdiv(P,Q)
%JSDIV   Jensen-Shannon divergence.
%   D = JSDIV(P,Q) calculates the Jensen-Shannon divergence of the two 
%   input distributions.

% JS-divergence
M = (P + Q) / 2;
D1 = KLdist(P,M);
D2 = KLdist(Q,M);
D = (D1 + D2) / 2;

% -------------------------------------------------------------------------
function D = KLdist(P,Q)
%KLDIST   Kullbach-Leibler distance.
%   D = KLDIST(P,Q) calculates the Kullbach-Leibler distance (information
%   divergence) of the two input distributions.

% KL-distance
P2 = P(P.*Q>0);     % restrict to the common support
Q2 = Q(P.*Q>0);
P2 = P2 / sum(P2);  % renormalize
Q2 = Q2 / sum(Q2);

D = sum(P2.*log(P2./Q2));