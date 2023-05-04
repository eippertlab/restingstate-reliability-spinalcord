clc
clear all
close all
%%
inDir  = '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/rawdata/';                            % raw .nii data
outDir = '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/';                       % derivatives (where to put output files)
codePath =  '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/Code/Preprocessing';
addpath(genpath(codePath))
cd(inDir)
subjects = dir('sub-ZS*');       % list of subjects
subjects([9,18,30]) = [];        % exclude subjects with ECG recording problems
sessions = {'auto', 'manual'}
%%
for sub = 1:size(subjects,1)
    
    
    cd(fullfile(outDir,subjects(sub).name,'func'))
    cd(fullfile('auto_pam50_templates','template'))
    system(['fslmaths ' ...
        ' PAM50_wm ' ...
        ' -bin ' ...
        ' PAM50_wm_bin ' ])
    
    system(['fslmaths ' ...
        ' PAM50_gm ' ...
        ' -bin ' ...
        ' PAM50_gm_bin ' ])
    
    system(['fslmaths ' ...
        ' PAM50_wm_bin ' ...
        ' -sub ' ...
        ' PAM50_gm_bin ' ...
        ' PAM50_wm_bin_thresholded'])
    
    system(['fslmaths ' ...
        ' PAM50_wm_bin_thresholded ' ...
        ' -thr 0.1 ' ...
        ' PAM50_wm_bin_thresholded'])
    
    
    cd(fullfile(outDir,subjects(sub).name,'func'))
    cd(fullfile('manual_pam50_templates','template'))
    system(['fslmaths ' ...
        ' PAM50_wm ' ...
        ' -bin ' ...
        ' PAM50_wm_bin ' ])
    
    system(['fslmaths ' ...
        ' PAM50_gm ' ...
        ' -bin ' ...
        ' PAM50_gm_bin ' ])
    
    system(['fslmaths ' ...
        ' PAM50_wm_bin ' ...
        ' -sub ' ...
        ' PAM50_gm_bin ' ...
        ' PAM50_wm_bin_thresholded'])
    
    system(['fslmaths ' ...
        ' PAM50_wm_bin_thresholded ' ...
        ' -thr 0.1 ' ...
        ' PAM50_wm_bin_thresholded'])
    
    cd(fullfile(outDir,subjects(sub).name,'func'))
    system(['fslmeants ' ...
        ' -i auto_moco2 ' ...
        ' -m ' fullfile('auto_pam50_templates','template','PAM50_wm_bin_thresholded') ...
        ' -o  auto_moco2_WM.txt '])
    
    system(['fslmeants ' ...
        ' -i manual_moco2 ' ...
        ' -m ' fullfile('manual_pam50_templates','template','PAM50_wm_bin_thresholded') ...
        ' -o  manual_moco2_WM.txt '])
    
end
%%
for sub = 1:size(subjects,1)
    
    
    cd(fullfile(outDir,subjects(sub).name,'func'))
    system(['paste manual_moco_outlies.txt manual_moco2_WM.txt > manual_mocooutliers_WM.txt'])
    system(['paste auto_moco_outlies.txt auto_moco2_WM.txt > auto_mocooutliers_WM.txt'])
    
    
    
end
%%
clc
clear all
close all

cd /data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives2/ICCandCon_ReviewAnalysis
load('slice-wise [mean] + avg & simple corr_nativespace_wmcorrs.mat')


DD = mean((DD_session1+DD_session2)./2)
VV = mean((VV_session1+VV_session2)./2)
wDV = mean((within_session1+within_session2)./2)
bDV = mean((between_session1+between_session2)./2)

%%
clc
clear all
close all

cd /data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives2/ICCandCon_ReviewAnalysis
load('slice-wise [mean] + avg & simple corr_nativespace_max_wmcorrs.mat')


DD = mean((DD_session1+DD_session2)./2)
VV = mean((VV_session1+VV_session2)./2)
wDV = mean((within_session1+within_session2)./2)
bDV = mean((between_session1+between_session2)./2)

