%3. Results

% In this code, results and figures associated with following section are created
% Section 3.3. Within-segment functional connectivity

% Please see:
% manuscript:  "" 
% and for the following dataset: https://openneuro.org/datasets/ds004386 


% Merve Kaptan, mkaptan@cbs.mpg.de
% 22.12.2022

%------------------------------------------
% Calculate the results & make figures

clc
clear all
close all
%% setup paths and directories

dataDir = '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives'
cd(dataDir)
subjects = dir('sub-ZS*');
% exclude problematic subjects 9, 18, 30!!! Because of bad ECG
subjects([9,18,30]) = [];

saveFolderName = '_Denoised_Max_Segmental'

subjects = {subjects.name};
saveDir = ['/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/ICCandCon' saveFolderName '/']
printPath = ['/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/Figures' saveFolderName '/'];

codeDir = '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/Code'
addpath(genpath(codeDir))

load('/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/subjectsSegmentalLevels.mat')
Levels = SubjectsSegLevels;

if ~exist(printPath)
    
    mkdir(printPath)
    
end

if ~exist(saveDir)
    
    mkdir(saveDir)
    
end

order = {'slice-wise [mean] + avg & simple corr'}

saveICC = 1;
saveCorrs = 1;
CalculateICC =1;
printMode = 1;
plotType = 'avgOverSes';
groupLevelICC = 1;

fileNames  = {'denoised_max'}

signalFolder = 'signalFiles';

analysisSpace = 'ns' %native space
sessions = {'manual_denoised', 'auto_denoised'}

if isequal(analysisSpace, 'ts')
    aspace = 'templspace'  %'templspace'  %or subjspace
elseif isequal(analysisSpace, 'ns')
    aspace = 'nativespace';
end



colorDD = [102,194,165]./255;
colorVV = [252,141,98]./255 ;
colorwDV =  [141,160,203]./255;
colorbDV = [231,138,195]./255;

legendNames = {'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'T1'};


o = 1;
f = 1;

%%
load([saveDir filesep order{o} '_' aspace '_' fileNames{f} 'ICC.mat']);

datasets1 = ICCvalue(1,:);
datasets1 = fliplr(datasets1);

datasets2 = squeeze(ICCvalue_CI(1,:,:));
datasets2 = flipud(datasets2);

colors = colorDD;

RELhelper_PlotICC(datasets1, datasets2, legendNames, colors, 'DD',subjects, order{1})


if printMode
    
    print([ printPath filesep ...
        'Figure_' sessions{:}  '_ICCSegment_DD_' num2str(numel(subjects)) '_' fileNames{f} '.tif'],'-dtiff')
    
    
end

datasets1 = ICCvalue(2,:);
datasets1 = fliplr(datasets1);

datasets2 = squeeze(ICCvalue_CI(2,:,:));
datasets2 = flipud(datasets2);

colors = colorVV;

RELhelper_PlotICC(datasets1, datasets2, legendNames, colors, 'VV',subjects, order{1})


if printMode
    
    print([ printPath filesep ...
        'Figure_' sessions{:}  '_ICCSegment_VV_' num2str(numel(subjects)) '_' fileNames{f} '.tif'],'-dtiff')
    
end

%%

datasets1 = ICCvalue(3,:);
datasets1 = fliplr(datasets1);

datasets2 = squeeze(ICCvalue_CI(3,:,:));
datasets2 = flipud(datasets2);

colors =  [141,160,203]./255;

RELhelper_PlotICC(datasets1, datasets2, legendNames, colors, 'wDV',subjects, order{1})


if printMode
    
    print([ printPath filesep ...
        'Figure_' sessions{:}  '_ICCSegment_wDV_' num2str(numel(subjects)) '_' fileNames{f} '.tif'],'-dtiff')
    
    
end

%%
datasets1 = ICCvalue(4,:);
datasets1 = fliplr(datasets1);

datasets2 = squeeze(ICCvalue_CI(4,:,:));
datasets2 = flipud(datasets2);

colors =  [231,138,195]./255;

RELhelper_PlotICC(datasets1, datasets2, legendNames, colors, 'bDV',subjects, order{1})


if printMode
    
    print([ printPath filesep ...
        'Figure_' sessions{:}  '_ICCSegment_bDV_' num2str(numel(subjects)) '_' fileNames{f} '.tif'],'-dtiff')
    
    
end



%%
load([saveDir filesep order{o} '_' aspace  '_' fileNames{f} 'corrs.mat']);

DD = (DD_session1 + DD_session2)./2;

DD = fliplr(DD);

fiGure = figure;
fiGure.Position = [0 0 800 800];
hold on
box off
YLims = [-.1 .6]
ylabel('r');
XLabels= legendNames;
xlim([0 numel(XLabels)+1]);
xticks([1:numel(XLabels)]);
ylim(YLims)
%yticks(-0.5:0.2:1)
ax = gca;
ax.FontSize = 15;
subtitle(order{o})
yline(0, 'LineWidth', 2)

h1 = boxplot(DD, ...
    'Colors', colorDD ,...
    'Symbol','.r', ...
    'Widths', 0.3);
set(h1,'linew',3)

lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');


MeanDD  = mean(DD, 'omitnan')
for ml = 1:numel(MeanDD)
    
    line([ml-0.15 ml+0.15], [MeanDD(ml) MeanDD(ml)], ...
        'Color', colorDD, 'LineWidth',3, 'LineStyle','-')
end


xticklabels(legendNames)


positions = [1:numel(legendNames)]


for l = 1:numel(legendNames)
    
    scatter(ones(1,size(DD,1)).*(positions(l)+(rand(1,size(DD,1))-0.5)/10),DD(:,l),[] ,colorDD ,'filled', ...
        'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.4)
    
end

ylim(YLims)
grid on
box off

if printMode
    
    print([ printPath filesep ...
        'Figure_' sessions{:}  '_ConnSegment_DD_' num2str(numel(subjects)) '_' fileNames{f} '.tif'],'-dtiff')
    
    
end


%%
VV = (VV_session1 + VV_session2)./2;

VV = fliplr(VV);


fiGure = figure;
fiGure.Position = [0 0 800 800];
hold on
box off
YLims = [-.1 .6]
ylabel('r');
XLabels= legendNames;
xlim([0 numel(XLabels)+1]);
xticks([1:numel(XLabels)]);
ylim(YLims)
%yticks(-0.5:0.2:1)
ax = gca;
ax.FontSize = 15;
subtitle(order{o})
yline(0, 'LineWidth', 2)

h1 = boxplot(VV, ...
    'Colors', colorVV ,...
    'Symbol','.r', ...
    'Widths', 0.3);
set(h1,'linew',3)

lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');



MeanVV  = mean(VV, 'omitnan');
for ml = 1:numel(MeanVV)
    
    line([ml-0.15 ml+0.15], [MeanVV(ml) MeanVV(ml)], ...
        'Color', colorVV, 'LineWidth',3, 'LineStyle','-')
end


xticklabels(legendNames)


positions = [1:numel(legendNames)]


for l = 1:numel(legendNames)
    
    scatter(ones(1,size(VV,1)).*(positions(l)+(rand(1,size(VV,1))-0.5)/10),VV(:,l),[] ,colorVV ,'filled', ...
        'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.4)
    
end
ylim(YLims)
grid on
box off
if printMode
    
    print([ printPath filesep ...
        'Figure_' sessions{:}  '_ConnSegment_VV_' num2str(numel(subjects)) '_' fileNames{f} '.tif'],'-dtiff')
    
    
end


%%

load([saveDir filesep order{o} '_' aspace  '_' fileNames{f} 'corrs.mat']);

wDV = (within_session1 + within_session2)./2;

wDV = fliplr(wDV);


fiGure = figure;
fiGure.Position = [0 0 800 800];
hold on
box off
YLims = [-.1 .6]
ylabel('r');
XLabels= legendNames;
xlim([0 numel(XLabels)+1]);
xticks([1:numel(XLabels)]);
ylim(YLims)
%yticks(-0.5:0.2:1)
ax = gca;
ax.FontSize = 15;
subtitle(order{o})
yline(0, 'LineWidth', 2)

h1 = boxplot(wDV, ...
    'Colors', colorwDV ,...
    'Symbol','.r', ...
    'Widths', 0.3);
set(h1,'linew',3)

lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');

MeanwDV  = mean(wDV, 'omitnan');
for ml = 1:numel(MeanwDV)
    
    line([ml-0.15 ml+0.15], [MeanwDV(ml) MeanwDV(ml)], ...
        'Color', colorwDV, 'LineWidth',3, 'LineStyle','-')
end

xticklabels(legendNames)


positions = [1:numel(legendNames)]


for l = 1:numel(legendNames)
    
    scatter(ones(1,size(wDV,1)).*(positions(l)+(rand(1,size(wDV,1))-0.5)/10),wDV(:,l),[] ,colorwDV ,'filled', ...
        'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.4)
    
end

grid on
box off
if printMode
    
    print([ printPath filesep ...
        'Figure_' sessions{:}  '_ConnSegment_wDV_' num2str(numel(subjects)) '_' fileNames{f} '.tif'],'-dtiff')
    
    
end




%%

bDV = (between_session1 + between_session2)./2;

bDV = fliplr(bDV);


fiGure = figure;
fiGure.Position = [0 0 800 800];
hold on
box off
YLims = [-.1 .6]
ylabel('r');
XLabels= legendNames;
xlim([0 numel(XLabels)+1]);
xticks([1:numel(XLabels)]);
ylim(YLims)
%yticks(-0.5:0.2:1)
ax = gca;
ax.FontSize = 15;
subtitle(order{o})
yline(0, 'LineWidth', 2)

h1 = boxplot(bDV, ...
    'Colors', colorbDV ,...
    'Symbol','.r', ...
    'Widths', 0.3);
set(h1,'linew',3)

lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');

MeanbDV  = mean(bDV, 'omitnan');
for ml = 1:numel(MeanbDV)
    
    line([ml-0.15 ml+0.15], [MeanbDV(ml) MeanbDV(ml)], ...
        'Color', colorbDV, 'LineWidth',3, 'LineStyle','-')
end

xticklabels(legendNames)


positions = [1:numel(legendNames)]


for l = 1:numel(legendNames)
    
    scatter(ones(1,size(bDV,1)).*(positions(l)+(rand(1,size(bDV,1))-0.5)/10),bDV(:,l),[] ,colorbDV ,'filled', ...
        'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.4)
    
end

grid on
box off


print([ printPath filesep ...
    'Figure_' sessions{:}  '_ConnSegment_bDV_' num2str(numel(subjects)) '_' fileNames{f} '.tif'],'-dtiff')


[row,col] =find(isnan(DD));

DD(row,:) = [];


csvwrite('tmp.csv', DD);
palm -i tmp.csv -twotail -quiet;
t_values = load('palm_dat_tstat_c1.csv');
p_values_c = load('palm_dat_tstat_fwep_c1.csv');

save([saveDir filesep 'corrPALM_DD_' fileNames{f} '.mat'],'t_values', 'p_values_c')
unix('rm -f tmp.csv');
unix('rm -f palm_*.csv');

[row,col] =find(isnan(VV));

VV(row,:) = [];


csvwrite('tmp.csv', VV);
palm -i tmp.csv -twotail -quiet;
t_values = load('palm_dat_tstat_c1.csv');
p_values_c = load('palm_dat_tstat_fwep_c1.csv');

save([saveDir filesep 'corrPALM_VV_' fileNames{f} '.mat'],'t_values', 'p_values_c')
unix('rm -f tmp.csv');
unix('rm -f palm_*.csv');

[row,col] =find(isnan(wDV));

wDV(row,:) = [];

csvwrite('tmp.csv', wDV);
palm -i tmp.csv -twotail -quiet;
t_values = load('palm_dat_tstat_c1.csv');
p_values_c = load('palm_dat_tstat_fwep_c1.csv');

save([saveDir filesep 'corrPALM_wDV_' fileNames{f} '.mat'],'t_values', 'p_values_c')
unix('rm -f tmp.csv');
unix('rm -f palm_*.csv');

[row,col] =find(isnan(bDV));

bDV(row,:) = [];


csvwrite('tmp.csv', bDV);
palm -i tmp.csv -twotail -quiet;
t_values = load('palm_dat_tstat_c1.csv');
p_values_c = load('palm_dat_tstat_fwep_c1.csv');

save([saveDir filesep 'corrPALM_bDV_' fileNames{f} '.mat'],'t_values', 'p_values_c')
unix('rm -f tmp.csv');
unix('rm -f palm_*.csv');
