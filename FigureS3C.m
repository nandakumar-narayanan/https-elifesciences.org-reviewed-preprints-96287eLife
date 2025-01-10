

t=readtable('PCAPharm.csv');

% Reviewers wanted PCA per animal
animals = unique(t.Animal);
clear data; 
figure(67); clf;  set(gcf, 'Color', 'White'); hold on; 
for i_animal = 1:length(animals)
subplot(1,11,i_animal); 
data{1} = t.PC1(t.Animal==i_animal&t.Type==1); 
data{2} = t.PC1(t.Animal==i_animal&t.Type==2); 
data{3} = t.PC1(t.Animal==i_animal&t.Type==3); 
jitterPlot(data, [0.1 0.1 0.1; 1 0 0; 0 0 1]);  
ylim([-4.5 6.5]);  % Make sure the range matches the min/max of PC1
box off; 
if i_animal==1; ylabel('PC1 Score'); else; set(gca, 'ytick', []);end;
set(gca, 'xtick', []); 
end



figure(68); 
subplot(121); hold on; 
scatter(t.PC1(t.Type==1),t.Slope(t.Type==1), 'k.', 'MarkerFaceColor', 'k'); lsline;
ylabel('GLM Slope'); xlabel('PC1');
xlabel('PC1 Score');
ylabel('Slope');

subplot(122); hold on; 
numClusters = 2;
points = [t.PC1(t.Type==1),t.Slope(t.Type==1)];
[clusterIdx, clusterCenters] = kmeans(points, numClusters);

pH = gscatter(points(:,1), points(:,2), clusterIdx, 'kk', '.*');

plot(clusterCenters(:,1), clusterCenters(:,2), 'k+', 'MarkerSize', 10, 'LineWidth', 2);
hold off;
title('K-means Clustering of (x, y) Points');
xlabel('PC1 Score');
ylabel('Slope');

t.Cluster = [clusterIdx; clusterIdx; clusterIdx];

figure(67); clf; set(gcf, 'color', 'White')
for i_animal = 1:length(animals)
subplot(1,11,i_animal); clear data; 
data{1} = t.PC1(t.Animal==i_animal&t.Type==1&t.Cluster==1); 
% data{2} = t.PC1(t.Animal==i_animal&t.Type==2&t.Cluster==1); 
% data{3} = t.PC1(t.Animal==i_animal&t.Type==3&t.Cluster==1); 
jitterPlot(data, [0.1 0.1 0.1; 1 0 0; 0 0 1]);  
ylim([-4.5 6.5]);  % Make sure the range matches the min/max of PC1
box off; 
if i_animal==1; ylabel('PC1 Score'); else; set(gca, 'ytick', []);end;
set(gca, 'xtick', 1, 'xticklabel', num2str(i_animal)); 
end

figure(66); clf; set(gcf, 'color', 'White')
for i_animal = 1:length(animals)
subplot(1,11,i_animal); clear data; 
data{1} = t.PC1(t.Animal==i_animal&t.Type==1&t.Cluster==2); 
% data{2} = t.PC1(t.Animal==i_animal&t.Type==2&t.Cluster==2); 
% data{3} = t.PC1(t.Animal==i_animal&t.Type==3&t.Cluster==2); 
jitterPlot(data, [0.1 0.1 0.1; 1 0 0; 0 0 1]);  
ylim([-4.5 6.5]);  % Make sure the range matches the min/max of PC1
box off; 
if i_animal==1; ylabel('PC1 Score'); else; set(gca, 'ytick', []);end;
set(gca, 'xtick', 1, 'xticklabel', num2str(i_animal)); 
end
