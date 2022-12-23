function sems = REL_helper_BarPlots(dataType,datasets,subPlot,noGroup,Lines, XLabels, YLabel,yLims,figTitle,printMode,printPath)

f = figure;
f.Position = [0 0 800 600];
hold on

for p = 1:numel(datasets)
    positions{p} = p;
end

colorZ = [0.4 0.4 0.4]; %[127,201,127]./255;

if noGroup
    
    bar(cell2mat(positions), mean(cell2mat(datasets),'omitnan'), ...
        'facecolor', [1 1 1], 'edgecolor', colorZ, ...
        'linewidth', 2);
    
    
    
    if Lines
        
        
        for j = 1:numel(datasets)
            
            sem= std(datasets{j}, 'omitnan')/sqrt(size(datasets{j},1), 'omitnan');
            e(j) = line([positions{j} positions{j}], ...
                [nanmean(datasets{j})-sem nanmean(datasets{j})+sem], ...
                'linewidth', 2.5, 'color', colorZ);
            sems(j) = sem;
            sem = [];
        end
        
        for p = 1:size(datasets,2)-1
            
            for j = 1:length(datasets{1})
                line([positions{p}+0.05 positions{p+1}-0.05], [datasets{p}(j) datasets{p+1}(j)], ...
                    'color', [colorZ 0.3], 'linewidth', 2.5);
            end
        end
        
    else
        
        tmp = cell2mat(datasets);
        
        
        for i = 1:size(tmp,2)
            
            scatter(ones(size(tmp,1),1).*(positions{i}+(rand(size(tmp,1),1)-0.5)/3),tmp(:,i), 60 ,[colorZ*1.5],'filled', ...
                'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
            
        end
        
        
        for j = 1:numel(datasets)
            
            sem= std(datasets{j})/sqrt(size(datasets{j},1));
            e(j) = line([positions{j} positions{j}], ...
                [nanmean(datasets{j})-sem nanmean(datasets{j})+sem], ...
                'linewidth', 2.5, 'color', colorZ);
            
            sems(j) = sem;
            sem = [];
        end
        
        
    end
    
else
    
    datasets_tmp = mean(cell2mat(datasets),'omitnan');
    
    datasets2= [datasets_tmp(1) datasets_tmp(2);datasets_tmp(1) datasets_tmp(3); datasets_tmp(1) datasets_tmp(4); ...
        datasets_tmp(1) datasets_tmp(5); datasets_tmp(1) datasets_tmp(6); datasets_tmp(1) datasets_tmp(7)];
    
    bar_handle = bar(1:size(datasets2,1), datasets2, ...
        'facecolor', [1 1 1], 'edgecolor', [0 0 0], ...
        'linewidth', 2);
    
    set(bar_handle(1),'FaceColor',[1,1,1],'EdgeColor', [0.5 0.5 0.5]);
    set(bar_handle(2),'FaceColor',[1,1,1],'EdgeColor', [0.3 0.3 0.3]);
    
    datasets_raw = {datasets{1},datasets{2}, ...
        datasets{1},datasets{3}, ...
        datasets{1}, datasets{4},...
        datasets{1}, datasets{5}, ...
        datasets{1}, datasets{6}...,
        datasets{1} datasets{7}};
    
    
    positions = {1,2,3,4,5,6};
    
    counter   = 0;
    
    
    
    for j = 1:numel(datasets_raw)
        
        sem= std(datasets_raw{j}, 'omitnan')/sqrt(size(datasets_raw{j},1));
        
        if ismember(j,[1 2])
            posit = 1;
        elseif ismember(j,[3 4])
            posit = 2;
        elseif ismember(j,[5 6])
            posit = 3;
        elseif ismember(j,[7 8])
            posit = 4;
        elseif ismember(j,[9 10])
            posit = 5;
        elseif ismember(j,[11 12])
            posit = 6;
            
        end
        
        if mod(j,2) ==0
            pp = [posit+0.15 posit+0.15];
        else
            pp = [posit-0.15 posit-0.15];
            
        end
        
        e(j) = line(pp, ...
            [nanmean(datasets_raw{j})-sem nanmean(datasets_raw{j})+sem], ...
            'linewidth', 2,'color', colorZ);
    end
    
    
    for p = [1:2:12]
        counter   = counter+1;
        for j = 1:length(datasets_raw{1})
            line([positions{counter}-0.09 positions{counter}+0.09], [datasets_raw{p}(j) datasets_raw{p+1}(j)], ...
                'color', [0.5 0.5 0.5 0.5], 'linewidth', 2);
        end
    end
    
end


ylim(yLims)
xticks(1:10);
xticklabels(XLabels)
set(gca,'TickLabelInterpreter','none')


ax = gca;
ax.FontSize = 14;
ylabel(YLabel, 'Interpreter', 'none')

title([figTitle  ' ,N = ' num2str(size(datasets{1},1)) ])
xticklabels(XLabels)
set(gca,'TickLabelInterpreter','none')


ax = gca;
ax.FontSize = 14;
ylabel(YLabel, 'Interpreter', 'none')

title([figTitle  ' ,N = ' num2str(size(datasets{1},1)) ])



if printMode
    
    print([ printPath filesep ...
        'Figure_' dataType '_' num2str(subPlot) figTitle 'N' num2str(size(datasets{1},1))  '_Grup_' num2str(noGroup) '.tif'],'-dtiff')
    
    print([ printPath filesep ...
        'Figure_' dataType '_' num2str(subPlot) figTitle 'N' num2str(size(datasets{1},1))  '_Grup_' num2str(noGroup) '.svg'],'-dsvg')
    
    
end


end