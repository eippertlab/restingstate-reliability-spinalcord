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
    system(['sct_apply_transfo ' ...
        ' -i auto_min.feat' filesep 'stats' filesep 'tsnr.nii.gz ' ...
        ' -d ' fullfile(outDir,'PAM50', 'template', 'PAM50_t2s.nii.gz')...
        ' -w warps' filesep 'warp_auto_moco2_mean2PAM50_t2.nii.gz ' ...
        ' -o TSNR_auto_normalized.nii.gz' ])
    
    system(['sct_apply_transfo ' ...
        ' -i manual_min.feat' filesep 'stats' filesep 'tsnr.nii.gz ' ...
        ' -d ' fullfile(outDir,'PAM50', 'template', 'PAM50_t2s.nii.gz')...
        ' -w warps' filesep 'warp_manual_moco2_mean2PAM50_t2.nii.gz ' ...
        ' -o TSNR_manual_normalized.nii.gz' ])
    

end
%%
for sub = 1:size(subjects,1)
    
        
        filename_auto{sub} = [subjects(sub).name filesep 'func' filesep 'TSNR_auto_normalized.nii.gz ' ];
           
        filename_manual{sub} = [subjects(sub).name filesep 'func' filesep 'TSNR_manual_normalized.nii.gz ' ];
      
      
        
  
end
%%
cd(outDir)

system(['fslmerge -t MergedTSNR_auto.nii.gz '  cell2mat(join(filename_auto))]);
system(['fslmerge -t MergedTSNR_manual.nii.gz '  cell2mat(join(filename_manual))]);



system(['fslmaths MergedTSNR_auto -Tmean MergedTSNR_auto_mean '])
system(['fslmaths MergedTSNR_manual -Tmean MergedTSNR_manual_mean '])
%%
cd(outDir)
system('fslmaths MergedTSNR_auto_mean -mul PAM50/template/PAM50_gm_bin MergedTSNR_auto_mean_gm ')
system('fslmaths MergedTSNR_manual_mean -mul PAM50/template/PAM50_gm_bin MergedTSNR_manual_mean_gm ')








