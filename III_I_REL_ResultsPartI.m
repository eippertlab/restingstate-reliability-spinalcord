%3. Results

% In this code, results and figures associated with following section are created
% Section 3.1.Replication and extension of previous resting-state functional connectivity results

% Please see:
% manuscript:  ""
% and for the following dataset: https://openneuro.org/datasets/ds004386


% Merve Kaptan, mkaptan@cbs.mpg.de
% 22.12.2022

%------------------------------------------
% Calculate the results & make figures

clc; clear all; close all
%% Setting up directories

saveFolderName = '_Denoised_NoThermal';  %extension of directory
saveDir = ['/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/ICCandCon' saveFolderName '/']; %directory to save
printPath = '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/FiguresManuscript/Figure1/'
codeDir = '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/Code'
fileNames  = {'max'}; % filename
connections = {'DD', 'VV', 'wDV', 'bDV'}; % which connections to run the test for
order = {'slice-wise [mean] + avg & simple corr'}; % type of FC calculation

addpath(genpath(codeDir))

cd(dataDir)
subjects = dir('sub-ZS*');

% exclude problematic subjects 9, 18, 30!!! Because of bad ECG
subjects([9,18,30]) = [];
subjects = {subjects.name};

analysisSpace = 'ns'; %native space

if isequal(analysisSpace, 'ts')
    aspace = 'templspace'
elseif isequal(analysisSpace, 'ns')
    aspace = 'nativespace';
end


if ~exist(printPath)
    
    mkdir(printPath)
    
end

if ~exist(saveDir)
    
    mkdir(saveDir)
    
end

printMode = 1;
o = 1;


colorDD = [102,194,165]./255;
colorVV = [252,141,98]./255 ;
colorwDV = [141,160,203]./255;
colorbDV = [231,138,195]./255;
%% Permuation test for functional connectivity strength
% FC values were calculated as part of II_IV_I_REL_Stats_FCandREL.m

for f = 1:numel(fileNames)
    
    load([saveDir filesep order{o} '_' aspace  '_' fileNames{f} 'corrs.mat']);
    
    DD = (DD_session1+DD_session2)./2;
    VV = (VV_session1+VV_session2)./2;
    wDV = (within_session1+within_session2)./2;
    bDV = (between_session1+between_session2)./2;
    
    corrDD(:,f) = DD;
    corrVV(:,f) = VV;
    corrwDV(:,f) = wDV;
    corrbDV(:,f) = bDV;
    
    clear DD_session1 DD_session2 VV_session1 VV_session2 within_session1 within_session2 ...
        between_session1 between_session2
    
end

for  c = 1:numel(connections)
    
    csvwrite('tmp.csv', eval(['corr' connections{c} ]));
    palm -i tmp.csv -twotail -quiet;
    t_values = load('palm_dat_tstat_c1.csv');
    p_values_c = load('palm_dat_tstat_fwep_c1.csv');
    
    save([saveDir filesep 'corrPALM_' connections{c}  '.mat'],'t_values', 'p_values_c')
    unix('rm -f tmp.csv');
    unix('rm -f palm_*.csv');
end


fiGure = figure;
fiGure.Position = [0 0 800 800];
hold on
box off

ylabel('r');
XLabels= fileNames;
xlim([0 numel(fileNames)+1]);
xticks([1:numel(fileNames)]);
ylim(YLims)
%yticks(-0.5:0.2:1)
ax = gca;
ax.FontSize = 15;
subtitle(order{o})
yline(0, 'LineWidth', 2)

positions = [1,2,3,4]


h1 = boxplot(DD, ...
    'Positions', ones(1,numel(DD)).*positions(1), ...
    'Colors', colorDD ,...
    'Symbol','.r', ...
    'Widths', 0.3);
set(h1,'linew',3)
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');
line([positions(1)-0.15 positions(1)+0.15], [mean(DD) mean(DD)],'Color', colorDD, 'LineWidth',3, 'LineStyle','-')



h2 = boxplot(VV, ...
    'Positions', ones(1,numel(VV)).*positions(2), ...
    'Colors', colorVV ,...
    'Symbol','.r', ...
    'Widths', 0.3);
set(h2,'linew',3)
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');
line([positions(2)-0.15 positions(2)+0.15], [mean(VV) mean(VV)],'Color', colorVV, 'LineWidth',3, 'LineStyle','-')


h3 = boxplot(wDV, ...
    'Positions', ones(1,numel(wDV)).*positions(3), ...
    'Colors', colorwDV ,...
    'Symbol','.r', ...
    'Widths', 0.3);
set(h3,'linew',3)
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');
line([positions(3)-0.15 positions(3)+0.15], [mean(wDV) mean(wDV)],'Color', colorwDV, 'LineWidth',3, 'LineStyle','-')

h4 = boxplot(bDV, ...
    'Positions', ones(1,numel(bDV)).*positions(4), ...
    'Colors', colorbDV ,...
    'Symbol','.r', ...
    'Widths', 0.3);
set(h4,'linew',3)
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');
line([positions(4)-0.15 positions(4)+0.15], [mean(bDV) mean(bDV)],'Color', colorbDV, 'LineWidth',3, 'LineStyle','-')


scatter(ones(size(DD)).*(positions(1)+(rand(size(DD))-0.5)/10),DD,[] ,colorDD ,'filled', 'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.4)
scatter(ones(size(VV)).*(positions(2)+(rand(size(VV))-0.5)/10),VV,[] ,colorVV ,'filled','MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.4)
scatter(ones(size(wDV)).*(positions(3)+(rand(size(wDV))-0.5)/10),wDV,[] ,colorwDV ,'filled','MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.4)
scatter(ones(size(bDV)).*(positions(4)+(rand(size(bDV))-0.5)/10),bDV,[] ,colorbDV ,'filled', 'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.4)


xticks(positions)
xticklabels({'Dorsal', 'Ventral', 'Within', 'Between'})
box off
ax = gca;
ax.FontSize = 20;
grid on

if printMode
    
    print([ printPath filesep ...
        'Figure_' sessions{:}  '_Connectivity' num2str(numel(subjects)) '_' fileNames{1} '.tif'],'-dtiff')
    
    
end
%%  Reliability and CIs
% ICCs results are calculated and saved as a part of II_IV_I_REL_Stats_FCandREL.m

load([saveDir filesep order{o} '_' aspace  '_' sessions{:}  '_' fileNames{f} 'ICC.mat']);
RELhelper_PlotICC(ICCvalue, ICCvalue_CI, {'Dorsal', 'Ventral', 'Within', 'Between'}, [0.4 0.4 0.4], 'Max Reliability', subjects)

if printMode
    
    print([ printPath filesep ...
        'Figure_' sessions{:}  '_Connectivity_ICC' num2str(numel(subjects)) '_' fileNames{1} '.tif'],'-dtiff')
    
end
%% tSNR
for sub = 1:numel(subjects)
    
    cd(fullfile(dataDir,subjects{sub},'func', 'signalFiles_tsnr'))
    
    for f = 1:numel(fileNames)
        
        for ses = 1:numel(sessions)
            
            [tmpLD,tmpRD, tmpLV, tmpRV,~] = REL_helper_loadRois(' ',sessions{ses}, fileNames{f}, 'tsnr', 'notWhole');
            TSNR_ld(sub,f,ses) = tmpLD;
            TSNR_rd(sub,f,ses) = tmpRD;
            TSNR_lv(sub,f,ses) = tmpLV;
            TSNR_rv(sub,f,ses) = tmpRV;
            clear tmpLD tmpRD tmpLV tmpRV
            
        end
        
    end
end

TSNR_ld = mean(squeeze(TSNR_ld),2);
TSNR_lv = mean(squeeze(TSNR_lv),2);
TSNR_rd = mean(squeeze(TSNR_rd),2);
TSNR_rv = mean(squeeze(TSNR_rv),2);

datasets{1} = TSNR_ld;
datasets{2} = TSNR_rd;
datasets{3} = TSNR_lv;
datasets{4} = TSNR_rv;

dataType = 'tsnr';
XLabels = fileNames;
subPlot = 1;
YLabel = 'Gray matter tSNR';
yLims = [0 25];
figTitle = 'mean tSNR_rois_figure1'
noGroup = 1;
Lines = 0;

meanDtsnr = mean((TSNR_ld + TSNR_rd)./2)
meanVtsnr = mean((TSNR_lv + TSNR_rv)./2)

REL_helper_BarPlots(dataType,datasets,subPlot, noGroup,Lines ,XLabels, YLabel,yLims,figTitle,printMode,printPath)