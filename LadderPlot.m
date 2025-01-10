function [ps] = LadderPlot(valuesOff, valuesOn, groupIdx, yaxislabel)

colormap = [0 0 0; .9 .1 .1; .1 .1 .75; .5 .5 .5; .4 .4 .4];
%Red is D2, Blue is D1 always, Grey is control fluorophore.
x = [1:2];
group_names = {'G1', 'G2'};
hold on; 


for i_An = 1:length(valuesOn);
hold on; 
plot(x, [median(valuesOff{i_An}) median(valuesOn{i_An})],'Color', colormap(1,:))
h = plot(1,median(valuesOff{i_An}),'wo'); 
set(h, 'MarkerFaceColor', colormap(1,:)); set(h, 'MarkerSize', 10)
h = plot(2,median(valuesOn{i_An}),'wo'); 
set(h, 'MarkerFaceColor', colormap(groupIdx+1,:)); set(h, 'MarkerSize', 10) ; 
end  


MqOff = median(cellfun(@median,valuesOff)); 
YqOff = quantile(cellfun(@median,valuesOff),[0.25 0.75]); 
MqOn = median(cellfun(@median,valuesOn)); 
YqOn = quantile(cellfun(@median,valuesOn),[0.25 0.75]); 


line([0.7 1.3], [MqOff MqOff], 'LineWidth', 5, 'Color', [0 0 0]);
line([1.7 2.3], [MqOn MqOn], 'LineWidth', 5, 'Color', [0 0 0]);

xlim([0.5 2.5]); set(gca, 'xtick', [1:2], 'xticklabel', group_names);
ylabel(yaxislabel); 

ps = signrank(cellfun(@median,valuesOn), cellfun(@median,valuesOff));
cd = cohend(cellfun(@median,valuesOn), cellfun(@median,valuesOff)); 

fprintf('\n >>> %s: Median (IQR):\n', yaxislabel); 
if MqOn > 1; 
fprintf('OFF: %.1f (%.1f-%.1f); \n ON: %.1f (%.1f-%.1f),\n p = %0.3f, Cohen d=%1.1f', MqOff, YqOff(1), YqOff(2),MqOn, YqOn(1), YqOn(2), ps, cd)  
else
fprintf('OFF: %0.2f (%0.2f-%0.2f); ON: %0.2f (%0.2f-%0.2f),\n p = %0.5f, Cohen d=%1.1f', MqOff, YqOff(1), YqOff(2),MqOn, YqOn(1), YqOn(2), ps, cd)  
end
fprintf('\n\n'); 
csvwrite(['./stats/' yaxislabel 'D' num2str(groupIdx) 'Off.csv'], cellfun(@median,valuesOff))
csvwrite(['./stats/' yaxislabel 'D' num2str(groupIdx) 'On.csv'], cellfun(@median,valuesOn))


end





