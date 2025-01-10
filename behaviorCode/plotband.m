function plotband(x, means, variance, color)
% Usage: plotband(X,Y,E, [0 0 1]);

% for compatibility with Kumar's plot band, no X is defined
% plotband(Y,E, [0 0 1]);
if nargin == 3
    color = variance;
    variance = means;
    means = x;
    x = 1:length(means);
end

% for dimension consistency
x = x(:);
means = means(:);
if size(variance,2) == size(x,1)
    variance = variance';
end
px = [x; flipud(x)];

% single variance like std or sem to dual variance like confidence intervals
if size(variance,2) == 1 
    stdvar = variance;
    variance = zeros(size(variance,1),2);
    variance(:,1) = means+stdvar;
    variance(:,2) = means-stdvar;
end

patch(px, [variance(:,1); flipud(variance(:,2))], color, 'EdgeAlpha', 0.1, 'FaceAlpha', 0.3);
hold on
plot(x, means, 'LineWidth', 2, 'Color',color)

end
