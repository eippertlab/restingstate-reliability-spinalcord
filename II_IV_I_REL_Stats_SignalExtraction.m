% 2.4. Statistical analysis
% No specific section for this the manuscript (link) --> but the time-courses from functional data 
% needs to be extracted before FC and ICC calculation

% Merve Kaptan, mkaptan@cbs.mpg.de
% 22.12.2022
%%
clc
clear all
close all
%% set directories, paths etc
outDir = '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/';     % derivatives (where to put output files)
codePath =  '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/Code/';
addpath(genpath(codePath))

cd(outDir)
subjects = dir('sub-ZS*');       % list of subjects
subjects([9,18,30]) = [];        % exclude subjects with ECG recording problems

sessions    = {'auto', 'manual'} ; % sessions

% all the files that were preprocessed
fileNames   = {'min','csf', 'moco', 'breath', 'cardiac', 'pnm','max', 'max_smooth','max_smooth4' , 'max_prewhiten', 'denoised_max'}

%% extract signal from preprocessed data for FC calculation
masks      = {'LD_thresholded_binarized', 'LV_thresholded_binarized', ...
    'RD_thresholded_binarized', 'RV_thresholded_binarized'}; % binary masks to use for signal extraction
txtOutDir   = 'signalFiles'; % output directory
niftiName   = 'res4d.nii.gz'; % input nifti
showAll     = 1; % extract the signal from individual voxels in the mask (1) or not
REL_helper_ExtractTimeSeries(outDir, subjects,niftiName,fileNames, masks, sessions, txtOutDir, showAll);
