%% sub functions
function jitterPlot(data, color)

if nargin==1; colorTable=[0.5 0.5 0.5]; 
else colorTable = color; end; 

colorTable = [colorTable;colorTable;     colorTable;     colorTable;     ];     
colorTable = [colorTable;colorTable;     colorTable;     colorTable;     ];     
colorTable = [colorTable;colorTable;     colorTable;     colorTable;     ];     

hold on
for count = 1:length(data)
    
        s = length(data{count}); % doing this in case you have a situation where a subject had no data.     
        se=nanstd(data{count})./sqrt(s);
        for i_line = 1:length(data{count}); 
            
                    pH = plot([count-0.15 + 0.3*i_line/length(data{count})],data{count}(i_line),'ko');   
                    set(pH, 'MarkerFaceColor', colorTable(count,:)); 
        end;

        means = nanmean(data{count});
        medians = nanmedian(data{count});

        pH = line([count-0.3 count+0.3], [medians medians]); set(pH, 'Color', [0 0 0], 'LineWidth', 4);


        
end




