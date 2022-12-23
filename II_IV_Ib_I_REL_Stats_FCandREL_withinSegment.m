% 2.4. Statistical analysis
% 2.4.1.1. Within-segment functional connectivity
% In this code, functional connectivity (for each segment) and its reliability were calculated for the following
% manuscript:  "" using the following dataset: https://openneuro.org/datasets/ds004386

% Merve Kaptan, mkaptan@cbs.mpg.de
% 22.12.2022
%%
clc
clear all
close all
%%
dataDir = '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives'
cd(dataDir)
subjects = dir('sub-ZS*');

% exclude problematic subjects 9, 18, 30!!! Because of bad ECG
subjects([9,18,30]) = [];
subjects = {subjects.name};


saveFolderName = '_Denoised_Max_Segmental' 
saveDir = ['/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/ICCandCon' saveFolderName '/']

codeDir = '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/Code'
addpath(genpath(codeDir))

load('/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/subjectsSegmentalLevels.mat')
Levels = SubjectsSegLevels

if ~exist(saveDir)
    
    mkdir(saveDir)
    
end

save(fullfile(saveDir,filesep,'subjects.mat'))

order = {'slice-wise [mean] + avg & simple corr'} 

saveICC = 1;
saveCorrs = 1;
CalculateICC =1 ;

fileNames  = {'denoised_max'};

signalFolder = 'signalFiles';

save(fullfile(saveDir,filesep,'fileNames.mat'))
save(fullfile(saveDir,filesep,'signalFolder.mat'))

analysisSpace = 'ns' %native space
sessions = {'manual_denoised', 'auto_denoised'}

if isequal(analysisSpace, 'ts')
    aspace = 'templspace'  %'templspace'  %or subjspace
elseif isequal(analysisSpace, 'ns')
    aspace = 'nativespace';
end

%%
for f = 1:numel(fileNames)
    
    for o = 1
        
        if isequal(order{o}, 'slice-wise [mean] & simple corr')
            
            [DD_session1,  DD_session2, VV_session1, VV_session2, ...
                LDV_session1,LDV_session2,RDV_session1,RDV_session2, ...
                within_session1, within_session2, between_session1, between_session2, ...
                ICCvalue, ICCvalue_CI ] = ...
                REL_helper_CalculateConnandRel_Segmental('slicewise',dataDir,subjects,analysisSpace, sessions, 'pearson',...
                'mean',fileNames{f}, CalculateICC, Levels);
                        
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
                REL_helper_CalculateConnandRel_Segmental('slicewise',dataDir,subjects,analysisSpace, sessions, 'partial', ...
                'mean',fileNames{f}, CalculateICC,Levels);
            
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
                REL_helper_CalculateConnandRel_Segmental('avgTS',dataDir,subjects,analysisSpace, sessions, 'pearson', ...
                'mean',fileNames{f}, CalculateICC, Levels);
            
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
                REL_helper_CalculateConnandRel_Segmental('avgTS',dataDir,subjects,analysisSpace, sessions, 'partial',...
                'mean', fileNames{f}, CalculateICC, Levels);
            
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
                REL_helper_CalculateConnandRel_Segmental('slicewise',dataDir,subjects,analysisSpace, sessions, 'partial', ...
                'voxel', CalculateICC, Levels);
            
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
                    REL_helper_CalculateConnandRel_Segmental('slicewise', dataDir,subjects,analysisSpace, sessions, 'pearson', ...
                    'mean',fileNames{f},signalFolder, CalculateICC, Levels);
                
                
            elseif isequal(order{o}, 'slice-wise [mean] + avg & partial corr')
                
                [DD_session1,  DD_session2, VV_session1, VV_session2, ...
                    LDV_session1,LDV_session2,RDV_session1,RDV_session2, ...
                    within_session1, within_session2, between_session1, between_session2, ...
                    ~, ~] = ...
                    REL_helper_CalculateConnandRel_Segmental('slicewise',dataDir,subjects,analysisSpace, sessions, 'partial', ...
                    'mean',fileNames{f}, CalculateICC, Levels);
                
                
                
            elseif isequal(order{o}, 'slice-wise [95%] + avg & partial corr')
                
                [DD_auto, DD_manual, VV_auto, VV_manual, ...
                    LDV_auto,LDV_manual,RDV_auto,RDV_manual, ...
                    within_auto, within_manual,...
                    between_auto, between_manual, ICCvalue, ICCvalue_CI ] = ...
                    REL_helper_CalculateConnandRel_Segmental('slicewise',dataDir,subjects,analysisSpace, sessions,...
                    'partial', 'voxel',fileNames{f}, CalculateICC, Levels);
                
                
            end
            
            if saveCorrs
                save([saveDir filesep order{o}  '_' aspace '_' fileNames{f} 'corrs.mat'] ...
                    ,'DD_session1',  'DD_session2', 'VV_session1', 'VV_session2', ...
                    'LDV_session1','LDV_session2','RDV_session1','RDV_session2', ...
                    'within_session1', 'within_session2', 'between_session1', 'between_session2')
            
                save([saveDir filesep order{o}  '_' aspace '_' fileNames{f} 'ICC.mat'] ...
                    ,'ICCvalue',  'ICCvalue_CI')
            
            end
            
        end
        
        
        
    end
    
end