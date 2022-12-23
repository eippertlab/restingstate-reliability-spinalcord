%3. Results

% In this code, results and figures associated with the sub-sections of following section are created
% Section 3.2.	 Impact of noise sources on resting-state functional connectivity and its reliability

% Please see:
% manuscript:  "" 
% and for the following dataset: https://openneuro.org/datasets/ds004386 

% Merve Kaptan, mkaptan@cbs.mpg.de
% 22.12.2022
%------------------------------------------
% Calculate the results & make figures
clc; clear all; close all
%% set up the directories
%PALM toolbox is needed to run non-parametric permutation tests
saveFolderName = '_Denoised_NoThermal' %extension of directory 
saveDir = ['/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/ICCandCon' saveFolderName '/']; %directory to save
savePath = '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/Physio/';
printPath = '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/FiguresManuscript/Figure2/'
printMode = 1;
codeDir = '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/Code'
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
%%  3.2.1. Physiological noise and amplitude of functional connectivity

% Calculate stats
o = 1
YLims = [-0.1 0.2]
fileNames   = {'min','csf', 'moco', 'breath', 'cardiac', 'pnm','max'}; % filenames


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

for c = 1:numel(connections)
    
    csvwrite('tmp.csv', eval(['corr' connections{c} ]));
    palm -i tmp.csv -twotail -quiet;
    t_values = load('palm_dat_tstat_c1.csv');
    p_values_c = load('palm_dat_tstat_fwep_c1.csv');
    
    save([saveDir filesep 'corrPALM_' connections{c}  '.mat'],'t_values', 'p_values_c')
    unix('rm -f tmp.csv');
    unix('rm -f palm_*.csv');
    
    
    tmp = eval(['corr' connections{c} '- ' 'corr' connections{c} '(:,1)' ]);
    tmp = tmp(:,2:end);
    
    csvwrite('tmp.csv', tmp)
    
    palm -i tmp.csv -twotail -quiet;
    t_values = load('palm_dat_tstat_c1.csv');
    p_values_c = load('palm_dat_tstat_fwep_c1.csv');
    
    save([saveDir filesep 'corrPALM_' connections{c} '_againstBaseline.mat'],'t_values', 'p_values_c')
    unix('rm -f tmp.csv');
    unix('rm -f palm_*.csv');
    
end

% Figure for FC
REL_helperPlotConnectivity(subjects, saveDir,YLims ,fileNames, analysisSpace, ...
    order,sessions,0,0,dataDir,'avgOverSes', printPath, printMode)



% Reliability
connNames = {'DD', 'VV', 'wDV', 'bDV'};

for c = 1:numel(connNames)
    
   [icc, icc_ci] =  REL_helperPlotFCReliability(subjects, saveDir, fileNames, analysisSpace, sessions, ...
        order,connNames{c}, printPath, printMode)
end

%% 3.2.2.  Physiological noise and reliability of functional connectivity
printPath = '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/FiguresManuscript/Figure3/'

load([savePath filesep 'physioRel_N' num2str(numel(subjects)) '.mat'])
xLabels = {'Refrms', 'Dvars', 'Avg HP', 'Std HP', 'Avg BP', 'Std BP'}

% Scatter plots

colors = [.4 .4 .4]
figure('Units','normalized','Position',[0 0 1 1])
hold on
subplot(2,4,2); hold on;title('dVars')
axis square

scatter(dVars(:,1), dVars(:,2),40, colors,'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
xlim([250 580])
ylim([250 580])
x = linspace(250,550);
y = x;
plot(x,y);
xticks([250:100:580])
yticks([250:100:580])

subplot(2,4,1); hold on;title('refrms')
axis square

scatter(refrms(:,1), refrms(:,2),40, colors,'filled', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
xlim([0.2 0.32])
ylim([0.2 0.32])
x = linspace(0.2,0.32);
y = x;
plot(x,y);
xticks([0.2:0.05:0.32])
yticks([0.2:0.05:0.32])

subplot(2,4,3); hold on;title('IBImean')
axis square

scatter(IBImean(:,1), IBImean(:,2),40, colors,'filled', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
xlim([650 1200])
ylim([650 1200])
x = linspace(650,1200);
y = x;
plot(x,y);
xticks([650:150:1250])
yticks([650:150:1250])

subplot(2,4,4); hold on;title('SDNN')
axis square

scatter(SDNN(:,1), SDNN(:,2),40, colors,'filled', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
xlim([10 150])
ylim([10 150])
x = linspace(10,150);
y = x;
plot(x,y);
xticks([10:50:150 150])
yticks([10:50:150 150])

subplot(2,4,5); hold on;title('rmsRR')
axis square

scatter(rmsRR(:,1), rmsRR(:,2),40, colors,'filled', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
xlim([10 150])
ylim([10 150])
x = linspace(10,150);
y = x;
plot(x,y);
xticks([10:50:150 150])
yticks([10:50:150 150])


subplot(2,4,6); hold on;title('Avg BP')
axis square

scatter(AvgBP(:,1), AvgBP(:,2),40, colors,'filled', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
xlim([2 8])
ylim([2 8])
x = linspace(0,8);
y = x;
plot(x,y);
xticks([2:2:8])
yticks([2:2:8])


subplot(2,4,7); hold on;title('Std BP')
axis square

scatter(StdBP(:,1), StdBP(:,2),40, colors,'filled', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
xlim([0 2.5])
ylim([0 2.5])
xlim([0 2.5])
ylim([0 2.5])
x = linspace(0,2.5);
y = x;
plot(x,y);
xticks([0:.5:2.5])
yticks([0:.5:2.5])
    
print([printPath filesep ...
        'Figure3A_Physio_N_' num2str(numel(subjects)) '_Data.tif'],'-dtiff')
    
% ICC plots
RELhelper_PlotICC(datasets1, datasets2, xLabels, [0.4 0.4 0.4], 'Physio Reliability', subjects)

% TSNR 
for sub = 1:numel(subjects)
    
    cd(fullfile(dataDir,subjects{sub},'func', 'signalFiles_tsnr'))
        
    for f = 1:numel(fileNames)
        
        for ses = 1:numel(sessions)
            
           [~,~, ~,~ ,tmpTSNR] = REL_helper_loadRois(' ',sessions{ses}, fileNames{f}, 'tsnr', 'whole');
            
           
           TSNR(sub,f,ses) = tmpTSNR;
            clear tmpLD tmpRD tmpLV tmpRV
            
        end
        
    end
end

TSNR = mean(TSNR,3,'omitnan');

for d = 1:size(TSNR,2)
    
    datasets{d} = TSNR(:,d);

end

dataType = 'tsnr';
XLabels  = fileNames;
subPlot  = 1;
YLabel   = 'mean tSNR';
yLims    =  [10 22]  
figTitle = 'mean tSNR'
noGroup  = 1; 
Lines    = 0

REL_helper_BarPlots(dataType,datasets,subPlot, noGroup,Lines ,XLabels, YLabel,yLims,figTitle,printMode,printPath)

% Explained variance
for sub = 1:numel(subjects)
    
    cd(fullfile(dataDir,subjects{sub},'func', 'signalFiles_expvar'))
        
    for f = 1:numel(fileNames)
        
        for ses = 1:numel(sessions)
           

                [~,~, ~,~ ,tmpEXPVAR] = REL_helper_loadRois(' ',sessions{ses}, fileNames{f}, 'expvar', 'whole');
            
           
             EXPVAR(sub,f,ses) = tmpEXPVAR;
             clear tmpEXPVAR
        
            
        end
        
    end
end

EXPVAR = mean(EXPVAR,3,'omitnan');

for d = 1:size(EXPVAR,2)
    
    datasets{d} = EXPVAR(:,d);

end

dataType = 'expvar';
XLabels = fileNames;
subPlot = 1;
YLabel = 'variance (R^2)';
yLims = [0 0.3];
figTitle = 'Explained Variance-whole cord'
noGroup = 1
Lines = 0

REL_helper_BarPlots(dataType,datasets,subPlot, noGroup, Lines ,XLabels, YLabel,yLims,figTitle,printMode,printPath)

for f = 1:numel(fileNames)
   
    tmp = squeeze(EXPVAR(:,f,:));
    ICCs_expvar(f,1) = ICC(tmp,'2','1');    
    clear tmp
    
end

RELhelper_PlotICC(ICCs_expvar, ICCs_expvar_CI, fileNames, [0.4 0.4 0.4], 'ExpVar Reliability', subjects)
%% 3.2.3. Thermal noise
% Calculate stats
o = 1
YLims = [-0.1 0.6];
fileNames   = {'max', 'denoised_max', 'max_smooth', 'max_smooth4'};
saveFolderName = '_Denoised_ThermalvsSmooth'
saveDir = ['/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/ICCandCon' saveFolderName '/']; %directory to save
printPath = '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/FiguresManuscript/Figure4/'


% TSNR 
for sub = 1:numel(subjects)
    
    cd(fullfile(dataDir,subjects{sub},'func', 'signalFiles_tsnr'))
        
    for f = 1:numel(fileNames)
        
        for ses = 1:numel(sessions)
            
           [~,~, ~,~ ,tmpTSNR] = REL_helper_loadRois(' ',sessions{ses}, fileNames{f}, 'tsnr', 'whole');
            
           
           TSNR(sub,f,ses) = tmpTSNR;
            clear tmpLD tmpRD tmpLV tmpRV
            
        end
        
    end
end

TSNR = mean(TSNR,3,'omitnan');

for d = 1:size(TSNR,2)
    
    datasets{d} = TSNR(:,d);

end

dataType = 'tsnr';
XLabels  = fileNames;
subPlot  = 1;
YLabel   = 'mean tSNR';
yLims    =  [10 22]  
figTitle = 'mean tSNR'
noGroup  = 1; 
Lines    = 0

REL_helper_BarPlots(dataType,datasets,subPlot, noGroup,Lines ,XLabels, YLabel,yLims,figTitle,printMode,printPath)


% Calculate stats
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

for c = 1:numel(connections)
    
    csvwrite('tmp.csv', eval(['corr' connections{c} ]));
    palm -i tmp.csv -twotail -quiet;
    t_values = load('palm_dat_tstat_c1.csv');
    p_values_c = load('palm_dat_tstat_fwep_c1.csv');
    
    save([saveDir filesep 'corrPALM_' connections{c}  '.mat'],'t_values', 'p_values_c')
    unix('rm -f tmp.csv');
    unix('rm -f palm_*.csv');
    
    
    tmp = eval(['corr' connections{c} '- ' 'corr' connections{c} '(:,1)' ]);
    tmp = tmp(:,2:end);
    
    csvwrite('tmp.csv', tmp)
    
    palm -i tmp.csv -twotail -quiet;
    t_values = load('palm_dat_tstat_c1.csv');
    p_values_c = load('palm_dat_tstat_fwep_c1.csv');
    
    save([saveDir filesep 'corrPALM_' connections{c} '_againstBaseline.mat'],'t_values', 'p_values_c')
    unix('rm -f tmp.csv');
    unix('rm -f palm_*.csv');
    
end

% Figure for FC
REL_helperPlotConnectivity(subjects, saveDir,YLims ,fileNames, analysisSpace, ...
    order,sessions,0,0,dataDir,'avgOverSes', printPath, printMode)

% Reliability
connNames = {'DD', 'VV', 'wDV', 'bDV'};

for c = 1:numel(connNames)
    
   [icc, icc_ci] =  REL_helperPlotFCReliability(subjects, saveDir, fileNames, analysisSpace, sessions, ...
        order,connNames{c}, printPath, printMode)
end

% Calculate spatial smoothness using AFNI's 3dFWHMx function to estimate
% spatial smoothness

[eSmooth_autoMax, ~] = thermalDenoise_3dFWMHx(outDir,subjects,'func/auto_max.feat/stats/',{'res4d.nii.gz'},...
    0, ...
    0, 0, '', '', 0, '',0,0)

[eSmooth_autoThermal, ~] = thermalDenoise_3dFWMHx(outDir,subjects,'func/auto_denoised_max.feat/stats/',{'res4d.nii.gz'},...
    0, ...
    0, 0, '', '', 0,'',0,0)

[eSmooth_manualMax, ~] = thermalDenoise_3dFWMHx(outDir,subjects,'func/manual_max.feat/stats/',{'res4d.nii.gz'},...
    0, ...
    0, 0, '', '', 0,'',0,0)


[eSmooth_manualThermal, ~] = thermalDenoise_3dFWMHx(outDir,subjects,'func/manual_denoised_max.feat/stats/',{'res4d.nii.gz'},...
    0, ...
    0, 0, '', '', 0,'',0,0)


[eSmooth_AutoSmooth2, ~] = thermalDenoise_3dFWMHx(outDir,subjects,'func/auto_max_smooth.feat/stats/',{'res4d.nii.gz'},...
    0, ...
    0, 0, '', '', 0,'',0,0)


[eSmooth_AutoSmooth4, ~] = thermalDenoise_3dFWMHx(outDir,subjects,'func/auto_max_smooth4.feat/stats/',{'res4d.nii.gz'},...
    0, ...
    0, 0, '', '', 0,'',0,0)


[eSmooth_ManualSmooth2, ~] = thermalDenoise_3dFWMHx(outDir,subjects,'func/manual_max_smooth.feat/stats/',{'res4d.nii.gz'},...
    0, ...
    0, 0, '', '', 0,'',0,0)


[eSmooth_ManualSmooth4, ~] = thermalDenoise_3dFWMHx(outDir,subjects,'func/manual_max_smooth4.feat/stats/',{'res4d.nii.gz'},...
    0, ...
    0, 0, '', '', 0,'',0,0)


OnlyMax = (eSmooth_autoMax+ eSmooth_manualMax)./2;
Thermal = (eSmooth_autoThermal + eSmooth_manualThermal)./2;
Smooth2 = (eSmooth_AutoSmooth2 + eSmooth_ManualSmooth2)./2;
Smooth4 = (eSmooth_AutoSmooth4 + eSmooth_ManualSmooth4)./2;


data1 = [OnlyMax, Thermal, Smooth2, Smooth4];

figPosition  = [0 0 600 600];
figColor     = [1 1 1];
figXLabel    = [];
figYLabel    = 'Effective smoothness';
figXLim      = [0.5 3];
figYLim      = [1 6];

positions1 = [1 1.5 2 2.5];
colors1 = {[ 0.5 0.5 0.5 ]};

figure; hold on; set(gcf, 'color', figColor , 'position', figPosition); box off
for c = 1:4
    
    h=boxplot(data1(:,c), 'positions',positions1(c), 'Colors',colors1{1} ,'Widths', [0.1] ,'Symbol','.r');
    set(h,'LineWidth',2 )
    
    clear h g
    
end

plot([positions1(1)+0.08 positions1(2)-0.08], mean(data1(:,1:2)), '--','LineWidth',2, 'color',[0 0  0])
plot([positions1(2)+0.08 positions1(3)-0.08], mean(data1(:,2:3)), '--','LineWidth',2, 'color',[0 0  0])
plot([positions1(3)+0.08 positions1(4)-0.08], mean(data1(:,3:4)), '--','LineWidth',2, 'color',[0  0 0])

plotSpread([0.3 0.3 0.3], data1, 'spreadWidth', 0.2, 'binWidth', 0.2, 'xValues', positions1- 0.2);

xticks(positions1)
yticks(1:6)
set(gca,'xticklabel',{'Max', '+Thermal', '+Smooth2', '+Smooth4'}, 'xlim', figXLim, 'ylim', figYLim)
grid on
box off
xlabel(figXLabel)
ylabel(figYLabel)