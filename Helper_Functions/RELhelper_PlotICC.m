function RELhelper_PlotICC(datasets1, datasets2, xLabels, colors, varargin)

switch numel(varargin)
    
    case 0
        
        titleTxt = [];
        subjects = [];
        subtitleTxt = [];
        
        
    case 1
        
        titleTxt = varargin{1};
        subjects = [];
        subtitleTxt = [];
        
    case 2
        
        titleTxt = varargin{1};
        subjects =  varargin{2};
        subtitleTxt = [];
        
    case 3
        
        titleTxt = varargin{1};
        subjects =  varargin{2};
        subtitleTxt = varargin{3};
        
        
end


fiGure = figure;
fiGure.Position = [0 0 800 800];
hold on

ylabel('ICC value');
yLimVal = [0 1];
ICCLims = [0.4 0.6 0.75 1 ];
xlim([0 numel(xLabels)+1]);
xticks([1:numel(xLabels)]);
ylim(yLimVal)
yticks(ICCLims)
xticklabels(xLabels)
set(gca,'TickLabelInterpreter','none')
ax = gca;
ax.FontSize = 20;
%title([titleTxt ' N = ' num2str(numel(subjects)) ])
%subtitle(subtitleTxt)


area([0:numel(xLabels)+1],ICCLims(1)*ones(1,numel([0:numel(xLabels)+1])), 'EdgeColor', [0.5 0.5 0.5], 'FaceColor' , [0.5 0.5 0.5], 'FaceAlpha', 0.05)
area([0:numel(xLabels)+1],ICCLims(2)*ones(1,numel([0:numel(xLabels)+1])), 'EdgeColor', [0.5 0.5 0.5], 'FaceColor' , [0.5 0.5 0.5], 'FaceAlpha', 0.05)
area([0:numel(xLabels)+1],ICCLims(3)*ones(1,numel([0:numel(xLabels)+1])), 'EdgeColor', [0.5 0.5 0.5], 'FaceColor' , [0.5 0.5 0.5], 'FaceAlpha', 0.05)
area([0:numel(xLabels)+1],ICCLims(4)*ones(1,numel([0:numel(xLabels)+1])), 'EdgeColor', [0.5 0.5 0.5], 'FaceColor' , [0.5 0.5 0.5], 'FaceAlpha', 0.05)


for i = 1:size(datasets2,1)
    
    plot([i i],datasets2(i,:)','color',colors,'LineWidth',3);
    hold on
    
end

for i=1:numel(datasets1)
    
    plot(i, datasets1(i), 'ko', 'LineWidth', 3, 'color',colors, 'MarkerSize', 20, 'MarkerFaceColor', colors)
    
end


end