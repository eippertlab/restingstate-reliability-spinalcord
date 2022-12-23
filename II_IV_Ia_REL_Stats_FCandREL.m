%2.4. Statistical analysis
%2.4.1. Functional connectivity and ICC calculation

% In this code, functional connectivity and its reliability were calculated for the following
% manuscript:  "" using the following dataset: https://openneuro.org/datasets/ds004386

% NOTE:
% run this script x4 times with different input data and save the results under different directories
% (see fileNames  saveFolderName  order  variables as defined below) 
%
% Merve Kaptan, mkaptan@cbs.mpg.de
% 22.12.2022

%%
clc
clear all
close all
%%
dataDir = '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives'
saveFolderName = '_Denoised_NoThermal'   % save results to this folder physiological 
                                          % noise correction
                                        
                    %'_Denoised_ThermalvsSmooth'    % save results to this folder for thermal noise
                    %'Denoised_Prewhiten'  % save results to this folder for comparison of prewhiten vs not 
                    %'Denoised_Partial'    % save results to this folder
                    %for comparison of partial vs pearson correlation

                                       
saveDir = ['/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/ICCandCon' saveFolderName '/'];
codeDir = '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/Code'
addpath(genpath(codeDir))



cd(dataDir)
subjects = dir('sub-ZS*');

subjects([9,18,30]) = []; %exclude problematic subjects 9, 18, 30!!! Because of bad ECG                                                  

subjects = {subjects.name};


if ~exist(saveDir)
    
    mkdir(saveDir)
    
end

save(fullfile(saveDir,filesep,'subjects.mat'))

%%
% type of connectivity calculation --> do slice-wise correlation
% calculations, take the mean. (simple corr --> Pearson correlation coefficient)
order = {'slice-wise [mean] + avg & simple corr'}
%{'slice-wise [mean] + avg & simple corr', 'slice-wise [mean] + avg & partial corr'} % use this for comparison of partial vs pearson correlation


saveICC = 1; % save ICC results or not
saveCorrs = 1; % save FC results or not
CalculateICC =1; % calculate ICC

fileNames = {'min','moco','csf', 'breath', 'cardiac', 'pnm','max'}; % physio noise correction to use for analysis
%{'max', 'denoised_max', 'max_smooth', 'max_smooth4'}; % use this for thermal noise section
%{'max', 'prewhiten_max'} %use this for supplementary results (prewhitening vs not)

if isequal(analysisSpace, 'ts')
    aspace = 'templspace'  %'templspace'  %or subjspace
elseif isequal(analysisSpace, 'ns')
    aspace = 'nativespace';
end

%%

for f = 1:numel(fileNames)
    
    for o = 1:numel(order)
        
        if isequal(order{o}, 'slice-wise [mean] & simple corr')
            
            [DD_session1,  DD_session2, VV_session1, VV_session2, ...
                LDV_session1,LDV_session2,RDV_session1,RDV_session2, ...
                within_session1, within_session2, between_session1, between_session2, ...
                ICCvalue, ICCvalue_CI ] = ...
                REL_helper_CalculateConnandRel('slicewise',dataDir,subjects,analysisSpace, sessions, 'pearson',...
                'mean',fileNames{f}, signalFolder, CalculateICC);
            
            if saveICC
                save([saveDir filesep order{o} '_' aspace '_' [sessions{:} ] '_' fileNames{f} 'ICC.mat'] ...
                    ,'ICCvalue', 'ICCvalue_CI')
            end
            
            if saveCorrs
                save([saveDir filesep order{o}  '_' aspace '_' fileNames{f} 'corrs.mat'] ...
                    ,'DD_session1',  'DD_session2', 'VV_session1', 'VV_session2', ...
                    'LDV_session1','LDV_session2','RDV_session1','RDV_session2', ...
                    'within_session1', 'within_session2', 'between_session1', 'between_session2')
            end
            
        elseif isequal(order{o}, 'slice-wise [mean] & partial corr')
            
            [DD_auto, DD_manual, VV_auto, VV_manual, ...
                LDV_auto,LDV_manual,RDV_auto,RDV_manual, ...
                within_auto, within_manual,...
                between_auto, between_manual, ICCvalue, ICCvalue_CI ] = ...
                REL_helper_CalculateConnandRel('slicewise',dataDir,subjects,analysisSpace, sessions, 'partial', ...
                'mean',fileNames{f}, signalFolder, CalculateICC);
            
            if saveICC
                save([saveDir filesep order{o} '_' aspace '_' [sessions{:} ] '_' fileNames{f} 'ICC.mat'] ...
                    ,'ICCvalue', 'ICCvalue_CI')
            end
            
            if saveCorrs
                save([saveDir filesep order{o}  '_' aspace '_' fileNames{f} 'corrs.mat'] ...
                    ,'DD_session1',  'DD_session2', 'VV_session1', 'VV_session2', ...
                    'LDV_session1','LDV_session2','RDV_session1','RDV_session2', ...
                    'within_session1', 'within_session2', 'between_session1', 'between_session2')
            end
            
        elseif isequal(order{o}, 'avg & simple corr')
            
            [DD_session1,  DD_session2, VV_session1, VV_session2, ...
                LDV_session1,LDV_session2,RDV_session1,RDV_session2, ...
                within_session1, within_session2, between_session1, between_session2, ...
                ICCvalue, ICCvalue_CI ] = ...
                REL_helper_CalculateConnandRel('avgTS',dataDir,subjects,analysisSpace, sessions, 'pearson', ...
                'mean',fileNames{f}, signalFolder, CalculateICC);
            
            if saveICC
                save([saveDir filesep order{o} '_' aspace '_' [sessions{:} ] '_' fileNames{f} 'ICC.mat'] ...
                    ,'ICCvalue', 'ICCvalue_CI')
            end
            
            if saveCorrs
                save([saveDir filesep order{o}  '_' aspace '_' fileNames{f} 'corrs.mat'] ...
                    ,'DD_session1',  'DD_session2', 'VV_session1', 'VV_session2', ...
                    'LDV_session1','LDV_session2','RDV_session1','RDV_session2', ...
                    'within_session1', 'within_session2', 'between_session1', 'between_session2')
            end
        elseif isequal(order{o}, 'avg & partial corr')
            
            [DD_session1,  DD_session2, VV_session1, VV_session2, ...
                LDV_session1,LDV_session2,RDV_session1,RDV_session2, ...
                within_session1, within_session2, between_session1, between_session2, ...
                ICCvalue, ICCvalue_CI ] = ...
                REL_helper_CalculateConnandRel('avgTS',dataDir,subjects,analysisSpace, sessions, 'partial',...
                'mean', fileNames{f}, signalFolder ,CalculateICC);
            
            if saveICC
                save([saveDir filesep order{o} '_' aspace '_' [sessions{:} ] '_' fileNames{f} 'ICC.mat'] ...
                    ,'ICCvalue', 'ICCvalue_CI')
            end
            
            if saveCorrs
                save([saveDir filesep order{o}  '_' aspace '_' fileNames{f} 'corrs.mat'] ...
                    ,'DD_session1',  'DD_session2', 'VV_session1', 'VV_session2', ...
                    'LDV_session1','LDV_session2','RDV_session1','RDV_session2', ...
                    'within_session1', 'within_session2', 'between_session1', 'between_session2')
            end
            
        elseif isequal(order{o}, 'slice-wise [95%] & partial corr')
            
            
            [DD_session1,  DD_session2, VV_session1, VV_session2, ...
                LDV_session1,LDV_session2,RDV_session1,RDV_session2, ...
                within_session1, within_session2, between_session1, between_session2, ...
                ICCvalue, ICCvalue_CI ] = ...
                REL_helper_CalculateConnandRel('slicewise',dataDir,subjects,analysisSpace, sessions, 'partial', ...
                'voxel', signalFolder ,CalculateICC);
            
            if saveICC
                save([saveDir filesep order{o} '_' aspace '_' [sessions{:} ] '_' fileNames{f} 'ICC.mat'] ...
                    ,'ICCvalue', 'ICCvalue_CI')
            end
            
            if saveCorrs
                save([saveDir filesep order{o}  '_' aspace '_' fileNames{f} 'corrs.mat'] ...
                    ,'DD_session1',  'DD_session2', 'VV_session1', 'VV_session2', ...
                    'LDV_session1','LDV_session2','RDV_session1','RDV_session2', ...
                    'within_session1', 'within_session2', 'between_session1', 'between_session2')
            end
        end
        
        
        
        if contains(order{o}, ' + avg')
            
            if isequal(order{o}, 'slice-wise [mean] + avg & simple corr')
                
                [DD_session1,  DD_session2, VV_session1, VV_session2, ...
                    LDV_session1,LDV_session2,RDV_session1,RDV_session2, ...
                    within_session1, within_session2, between_session1, between_session2, ...
                    ICCvalue, ICCvalue_CI ] = ...
                    REL_helper_CalculateConnandRel('slicewise', dataDir,subjects,analysisSpace, sessions, 'pearson', ...
                    'mean',fileNames{f},signalFolder, CalculateICC);
                
                
            elseif isequal(order{o}, 'slice-wise [mean] + avg & partial corr')
                
                [DD_session1,  DD_session2, VV_session1, VV_session2, ...
                    LDV_session1,LDV_session2,RDV_session1,RDV_session2, ...
                    within_session1, within_session2, between_session1, between_session2, ...
                    ~, ~] = ...
                    REL_helper_CalculateConnandRel('slicewise',dataDir,subjects,analysisSpace, sessions, 'partial', ...
                    'mean',fileNames{f},signalFolder, CalculateICC);
                
                
                
            elseif isequal(order{o}, 'slice-wise [95%] + avg & partial corr')
                
                [DD_auto, DD_manual, VV_auto, VV_manual, ...
                    LDV_auto,LDV_manual,RDV_auto,RDV_manual, ...
                    within_auto, within_manual,...
                    between_auto, between_manual, ICCvalue, ICCvalue_CI ] = ...
                    REL_helper_CalculateConnandRel('slicewise',dataDir,subjects,analysisSpace, sessions,...
                    'partial', 'voxel',fileNames{f}, signalFolder,CalculateICC);
                
                
            end
            
            
            
            
            DD_session1 = nanmean(DD_session1,2);
            DD_session2 = nanmean(DD_session2,2);
            VV_session1   = nanmean(VV_session1,2);
            VV_session2 = nanmean(VV_session2,2)
            LDV_session1  = nanmean(LDV_session1,2);
            LDV_session2 = nanmean(LDV_session2,2);
            RDV_session1   = nanmean(RDV_session1,2);
            RDV_session2   = nanmean(RDV_session2,2);
            within_session1   = nanmean(within_session1,2);
            within_session2 = nanmean(within_session2,2);
            between_session1   = nanmean(between_session1,2);
            between_session2 =  nanmean(between_session2,2);
            
            
            if CalculateICC
                X = [DD_session1, DD_session2];
                ICCvalue(1,1) =  ICC(X, '2', '1');
                ICCvalue_CI(1,:) = bootci(100000,@myICC,X);
                
                clear X
                X = [VV_session1, VV_session2];
                ICCvalue(2,1) =  ICC(X, '2', '1');
                ICCvalue_CI(2,:) = bootci(100000,@myICC,X);
                
                
                clear X
                X = [within_session1, within_session2];
                ICCvalue(3,1) =  ICC(X, '2', '1');
                ICCvalue_CI(3,:) = bootci(100000,@myICC,X);
                
                clear X
                X = [between_session1, between_session2];
                ICCvalue(4,1) =  ICC(X, '2', '1');
                ICCvalue_CI(4,:) = bootci(100000,@myICC,X);
            end
            
            if saveICC
                save([saveDir filesep order{o} '_' aspace '_' [sessions{:} ] '_' fileNames{f} 'ICC.mat'] ...
                    ,'ICCvalue', 'ICCvalue_CI')
            end
            
        end
        
        
        if saveCorrs
            save([saveDir filesep order{o}  '_' aspace '_' fileNames{f} 'corrs.mat'] ...
                ,'DD_session1',  'DD_session2', 'VV_session1', 'VV_session2', ...
                'LDV_session1','LDV_session2','RDV_session1','RDV_session2', ...
                'within_session1', 'within_session2', 'between_session1', 'between_session2')
        end
        
    end
    
    
    
end

