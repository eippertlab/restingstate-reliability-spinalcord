% In this code, the preprocessing steps are listed for the following
% manuscript:  "" and for the following dataset:
% https://openneuro.org/datasets/ds004386
% In order to run this code the Spinal Cord Toolbox (version 4.2.2 or higher -
% note that for the manuscript version 4.2.2 was used)
% and FSL should be added to the bash profile of the user to call their
% functions from MATLAB.

% When manual processing was performed, the respective outputs can be found here: 
% https://openneuro.org/datasets/ds004386
% under the "derivatives" parent directory and subject specific
% subdirectories

% Merve Kaptan, mkaptan@cbs.mpg.de
% 22.12.2022

clc
clear all
close all
%%
inDir           = '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/rawdata/';                    % raw .nii data
outDir          = '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/';                % derivatives (where to put output files)
scttemplatepath = '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/template/SCT/';  % SCT template folder
funcsegQCdir    = '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/funcsegQC/';    % QC file for segmentation
templateDir     = '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/PAM50/spinal_levels/';
codePath        =  '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/Code/';
addpath(genpath(codePath))

if ~exist(funcsegQCdir)
    mkdir(funcsegQCdir)
end


if ~exist(outDir)
    mkdir(outDir)
end

cd(inDir)
subjects = dir('sub-ZS*');       % list of subjects
subjects([9,18,30]) = [];        % exclude subjects with ECG recording problems

sessions   = {'manual', 'auto'}   % 2 sessions
numSlices  = 24;                  % number of EPI slices

%set FSL environment
%Running FSL from matlab
setenv('FSLDIR', '/afs/cbs.mpg.de/software/fsl/5.0.11/ubuntu-xenial-amd64/');
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be
%
%% Step 1:
%  _   _   _   _   _   _   _   _   _   _  
% / \ / \ / \ / \ / \ / \ / \ / \ / \ / \
%( a | n | a | t | o | m | i | c | a | l |
% \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ 


% Segmentation of T2w image

for sub = 1:size(subjects,1)
    
    subid = subjects(sub).name;
    
    % create the output directory
    outdir = fullfile(outDir,subid,'anat');
    mkdir(outdir)
    
    % copy the T2-weighted anatomical image there and cd to that directory
    cd(fullfile(inDir,subid))
    copyfile([ 'anat' filesep '*_T2w.nii.gz*'], [outdir filesep]);
    cd(outdir)
    
    % make a tmp directory to keep the folder organization
    mkdir('tmp')
    copyfile('*_T2w*.nii.gz', ['tmp' filesep])
    cd tmp
    
    % first segment the T2 initially
    % then smooth this segmentation with 8mm kernel
    % and use the smoothed image for final segmentation
    system('sct_deepseg_sc -i    *_T2w*  -c t2 -brain 0 ' );
    system('sct_smooth_spinalcord -i *_T2w.nii.gz* -s  *seg.nii.gz* -smooth 0,0,8 ');
    system('sct_deepseg_sc -i *smooth* -c t2 -brain 0 ');
    
    % move the main segmentation to the parent folder
    movefile('*T2w_smooth_seg*', [outdir filesep] )
    
    cd(outdir)
    % remove unnecessary ones
    system('rm -rf tmp')
end

% Registration of T2w image to the template space
for sub = 1:size(subjects,1)
    
    subid = subjects(sub).name;
    
    % go to the output directory
    outdir = fullfile(outDir,subid,'anat');
    cd(outdir)
    
    % label vertebrae automatically
    system(['sct_label_vertebrae -i *T2w.nii.gz* -s ' outdir ...
        filesep '*T2w_smooth_seg.nii* -c t2 ']);
    % remove unnecessary labels
   system (['sct_label_utils -i *T2w_smooth_seg_labeled_discs.nii* '...
        '-keep 2,3,4,5,6,7,8,9 -o ' outdir  filesep 'disc_labels.nii.gz'])
    % note that the  labels were manually corrected if necessary (the name of
    % the file was not changed! The labels used for the normalization can be found 
    % under derivatives/subid/anat/ folder)
    
    % register T2w image to the template
     system(['sct_register_to_template -i *T2w.nii.gz* -s *_smooth_seg.nii* -ldisc disc_labels_edited.nii.gz ' ...
        ' -param step=1,type=seg,algo=slicereg,metric=MeanSquares,' ...
        'iter=10,smooth=2,gradStep=0.5,slicewise=0,smoothWarpXY=2,' ...
        'pca_eigenratio_th=1.6:step=2,type=seg,algo=bsplinesyn,metric=MeanSquares,' ...
        'iter=3,smooth=1,gradStep=0.5,slicewise=0,smoothWarpXY=2,pca_eigenratio_th=1.6 ' ...
        ' -c t2']);
    % rename the normalized T2w image
    system('mv anat2template.nii.gz T2w_normalized.nii.gz')
        
    % remove unnecessary files
    delete('*T2w.nii.gz*')
    delete('*.cache*')
    delete('straight_ref.nii.gz')
    
    % make a directory for warping fields for organizational purposes
    mkdir('warps')
    system(['mv *warp*.nii* warps' filesep])
    % rename the warp fields
    cd warps
    system('mv warp_template2anat.nii.gz warp_PAM502T2.nii.gz')
    system('mv warp_anat2template.nii.gz warp_T22PAM50.nii.gz')

end

%% STEP 2: 
%  _   _   _   _   _   _   _   _   _   _  
% / \ / \ / \ / \ / \ / \ / \ / \ / \ / \
%( f | u | n | c | t | i | o | n | a | l |
% \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ 

%  _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _
% / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \
%( m | o | t | i | o | n | - | c | o | r | r | e | c | t | i | o | n )
% \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/

% Motion-correction STEP 1


for sub = 1:size(subjects,1)
    
    for ses = 1:numel(sessions)
        
        cd(fullfile(inDir,subjects(sub).name,'func'))
        
        fprintf (['subject ' subjects(sub).name '--START'])
        
        mkdir(fullfile(outDir,subjects(sub).name,'func'))
        
        system(['cp *' sessions{ses} '*.nii.gz '  ...
            fullfile(outDir,subjects(sub).name,'func', [sessions{ses} '_raw.nii.gz']) ]);
        
        cd(fullfile(outDir,subjects(sub).name,'func'))
        
        system(['fslmaths ' sessions{ses} '_raw.nii.gz -Tmean '   sessions{ses} '_raw_mean.nii.gz' ]);
        
        system(['sct_propseg -i ' sessions{ses} '_raw_mean.nii.gz -c t2s' ])
        
        system(['fslmeants -i ' sessions{ses} '_raw_mean  -m ' sessions{ses} '_raw_mean_seg  --showall -o tmp.txt' ]);
        
        % Load data (remember that first 3 rows are coordinates)
        data     = load('tmp.txt');
        
        if length(unique(data(3,:)))== numSlices
            
            fprintf('Detected %d slices in the mask, no adjustment needed \n', numSlices)
        else
            warning ('Not enough slices [ number of slices in the mask =  %d ], the mask will be adjusted \n', length(unique(data(3,:))))
            % In the following part, if the automatic segmentation does not
            % propogate the options for sct_propseg function that affect the z-propogation
            % will be modified.
            
            % if the segmentation does not propogate, first modify the -max-area
            % that affects maximum cross-sectional area (see sct_propseg for details)
            unix(['sct_propseg -i ' sessions{ses} '_raw_mean.nii.gz  -c t2s -radius 5 -max-area 150']);
            system(['fslmeants -i ' sessions{ses} '_raw_mean  -m ' sessions{ses} '_raw_mean_seg  --showall -o tmp.txt' ]);
            % Load data (remember that first 3 rows are coordinates)
            data     = load('tmp.txt');
            
            if length(unique(data(3,:)))~= numSlices
                % if the segmentation still does not propogate,  modify the -max-area
                % that affects maximum cross-sectional area and min-contrast (see sct_propseg for
                % details)
                unix(['sct_propseg -i ' sessions{ses} '_raw_mean.nii.gz -c t2s -radius 5 -max-area 150 -min-contrast 30']);
                system(['fslmeants -i ' sessions{ses} '_raw_mean  -m ' sessions{ses} '_raw_mean_seg  --showall -o tmp.txt' ]);
                % Load data (remember that first 3 rows are coordinates)
                data     = load('tmp.txt');
                
                if length(unique(data(3,:)))~= numSlices
                    % if the segmentation still does not propogate,  modify the -max-area
                    % that affects maximum cross-sectional area, min-contrast and max-deformation (see sct_propseg for
                    % details)
                    [status, ~] = unix(['sct_propseg -i ' sessions{ses} '_raw_mean.nii.gz -c t2s -radius 5 -max-area 150 -min-contrast 30 -max-deformation 5']);
                    system(['fslmeants -i ' sessions{ses} '_raw_mean  -m ' sessions{ses} '_raw_mean_seg  --showall -o tmp.txt' ]);
                    % Load data (remember that first 3 rows are coordinates)
                    data     = load('tmp.txt');
                    if length(unique(data(3,:)))~= numSlices
                        error ('Not enough slices in the mask, please manually adjust your mask!!')
                        % If the image contrast is low or automatic segmentation
                        % does not work for some reason, the mask needs to be
                        % manually adjusted (fsleyes can be used).
                    end
                end
                
            end
            
        end
        
        if length(unique(data(3,:)))== numSlices
            
            system(['sct_create_mask -i ' sessions{ses} '_raw_mean.nii.gz -p centerline,' sessions{ses} '_raw_mean_seg.nii.gz -o ' sessions{ses} '_raw_mean_mask.nii.gz' ])
            system(['fslmerge -tr ' sessions{ses} '_raw_mean_merged '  sessions{ses} '_raw_mean '  sessions{ses}  '_raw 2.312 '  ])
            
            system(['sct_fmri_moco ' ...
                '-i '  sessions{ses} '_raw_mean_merged.nii.gz ' ...
                '-m ' sessions{ses} '_raw_mean_mask.nii.gz ' ...
                '-param iterAvg=0 -x spline']);
            
            system(['sct_image ' ...
                ' -i '  sessions{ses} '_raw_mean_merged_moco.nii.gz ' ...
                ' -remove-vol 0 ' ...
                ' -o ' sessions{ses} '_moco.nii.gz']);
            
            system(['fslmaths '  sessions{ses} '_moco.nii.gz' ...
                ' -Tmean ' ...
                sessions{ses} '_moco_mean '])
            
        else
            mocoProblemSubject{sub,ses} = subjects(sub).name;
            
            
        end
        
        fprintf (['subject ' subjects(sub).name '--DONE'])
        
        
    end
end


% Motion-correction: STEP 2

for sub = 1:size(subjects,1)
    
    for ses = 1:numel(sessions)
        
        
        fprintf (['subject ' subjects(sub).name '--START' newline])
        
        cd(fullfile(outDir,subjects(sub).name,'func'))
        
        system('rm *_T0*.nii.gz')
        
        % take the mean of moco images from Step 1
        system(['fslmaths ' sessions{ses} '_moco.nii.gz -Tmean '   sessions{ses} '_moco_mean.nii.gz' ]);
        
        %segment this mean image
        system(['sct_propseg -i ' sessions{ses} '_moco_mean.nii.gz -c t2s' ])
        
        %extract the signal using this segmentation (to check if segmentation propagated)
        system(['fslmeants -i ' sessions{ses} '_moco_mean ' ...
            '  -m ' sessions{ses} '_moco_mean_seg ' ...
            ' --showall -o tmp.txt' ]);
        
        % Load data (remember that first 3 rows are coordinates)
        data     = load('tmp.txt');
        
        if length(unique(data(3,:)))== numSlices
            
            fprintf('Detected %d slices in the mask, no adjustment needed \n', numSlices)
        else
            warning ('Not enough slices [ number of slices in the mask =  %d ], the mask will be adjusted \n', length(unique(data(3,:))))
            % In the following part, if the automatic segmentation does not
            % propogate the options for sct_propseg function that affect the z-propogation
            % will be modified.
            
            % if the segmentation does not propogate, first modify the -max-area
            % that affects maximum cross-sectional area (see sct_propseg for details)
            unix(['sct_propseg -i ' sessions{ses} '_moco_mean.nii.gz  -c t2s -radius 5 -max-area 150']);
            system(['fslmeants -i ' sessions{ses} '_moco_mean  -m ' sessions{ses} '_moco_mean_seg  --showall -o tmp.txt' ]);
            % Load data (remember that first 3 rows are coordinates)
            data     = load('tmp.txt');
            
            if length(unique(data(3,:)))~= numSlices
                % if the segmentation still does not propogate,  modify the -max-area
                % that affects maximum cross-sectional area and min-contrast (see sct_propseg for
                % details)
                unix(['sct_propseg -i ' sessions{ses} '_moco_mean.nii.gz -c t2s -radius 5 -max-area 150 -min-contrast 30']);
                system(['fslmeants -i ' sessions{ses} '_moco_mean  -m ' sessions{ses} '_moco_mean_seg  --showall -o tmp.txt' ]);
                % Load data (remember that first 3 rows are coordinates)
                data     = load('tmp.txt');
                
                if length(unique(data(3,:)))~= numSlices
                    % if the segmentation still does not propogate,  modify the -max-area
                    % that affects maximum cross-sectional area, min-contrast and max-deformation (see sct_propseg for
                    % details)
                    [status, ~] = unix(['sct_propseg -i ' sessions{ses} '_moco_mean.nii.gz -c t2s -radius 5 -max-area 150 -min-contrast 30 -max-deformation 5']);
                    system(['fslmeants -i ' sessions{ses} '_moco_mean  -m ' sessions{ses} '_moco_mean_seg  --showall -o tmp.txt' ]);
                    % Load data (remember that first 3 rows are coordinates)
                    data     = load('tmp.txt');
                    if length(unique(data(3,:)))~= numSlices
                        error ('Not enough slices in the mask, please manually adjust your mask!!')
                        % If the image contrast is low or automatic segmentation
                        % does not work for some reason, the mask needs to be
                        % manually adjusted (fsleyes can be used).
                    end
                end
                
            end
            
        end
        
        if length(unique(data(3,:)))== numSlices %if segmenation propagates
            
            %create a mask to be used in moco
            system(['sct_create_mask ' ...
                ' -i ' sessions{ses} '_moco_mean.nii.gz ' ...
                ' -p centerline,' sessions{ses} '_moco_mean_seg.nii.gz '...
                ' -o ' sessions{ses} '_moco_mean_mask.nii.gz' ])
            
            %append moco mean image (target of moco) to the RAW images
            system(['fslmerge -tr ' sessions{ses} '_moco_mean_merged ' ...
                sessions{ses} '_moco_mean '  sessions{ses}  '_raw 2.312 '  ])
            
            %do the moco
            system(['sct_fmri_moco ' ...
                '-i '  sessions{ses} '_moco_mean_merged.nii.gz ' ...
                '-m ' sessions{ses} '_moco_mean_mask.nii.gz ' ...
                '-param iterAvg=0 -x spline']);
            
            %remove the appended mean image
            system(['sct_image ' ...
                ' -i '  sessions{ses} '_moco_mean_merged_moco.nii.gz ' ...
                ' -remove-vol 0 ' ...
                ' -o ' sessions{ses} '_moco2.nii.gz']);
            
            %take the mean
            system(['fslmaths '  sessions{ses} '_moco2.nii.gz' ...
                ' -Tmean ' ...
                sessions{ses} '_moco2_mean '])
            
        else
            mocoProblemSubject{sub,ses} = subjects(sub).name;
            
            
        end
        
        fprintf (['subject ' subjects(sub).name '--DONE'])
        
    end
end
%% _   _   _   _     _   _   _   _   _   _   _   _   _   _     _   _   _
% / \ / \ / \ / \   / \ / \ / \ / \ / \ / \ / \ / \ / \ / \   / \ / \ / \
%( m | o | c | o ) ( p | a | r | a | m | a | t | e | r | s ) ( a | n | d )
% \_/ \_/ \_/ \_/   \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/   \_/ \_/ \_/
%  _   _   _   _   _   _   _   _
% / \ / \ / \ / \ / \ / \ / \ / \
%( o | u | t | l | i | e | r | s )
% \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/


% remove the first volume (target volume) from -x and -y parameter files obtained from SCT
for sub = 1:size(subjects,1)
    
    fprintf (['subject ' subjects(sub).name '--START' newline])
    
    for ses = 1:numel(sessions)
        
        cd(fullfile(outDir,subjects(sub).name,'func'))
        
        system(['sct_image ' ...
            ' -i ' sessions{ses} '_moco_mean_merged_moco_params_X.nii.gz'...
            ' -remove-vol 0 ' ...
            ' -o ' sessions{ses} '_moco2_paramX.nii.gz'])
        
        system(['sct_image ' ...
            ' -i ' sessions{ses} '_moco_mean_merged_moco_params_Y.nii.gz'...
            ' -remove-vol 0 ' ...
            ' -o ' sessions{ses} '_moco2_paramY.nii.gz'])
        
    end
    
end

%% NOTE

% Remember: for subject sub-ZS003 auto session the recording of triggers
% started a bit later than the functional acquisition. Remove the first 9
% volumes for which we do not have BrainAmp triggers. Also remember: in FSL,
% counting starts from 0 - not from 1.

cd(fullfile(outDir, 'sub-ZS003', 'func'))

system('fslroi auto_moco2.nii.gz auto_moco2.nii.gz 0 -1 0 -1 0 -1 9 -1')
system('fslroi auto_moco2_paramX.nii.gz  auto_moco2_paramX.nii.gz 0 -1 0 -1 0 -1 9 -1')
system('fslroi auto_moco2_paramY.nii.gz auto_moco2_paramY.nii.gz 0 -1 0 -1 0 -1 9 -1')

% calculate DVARS and REFRMS (note that REFRMS function is slightly
% modified by Dr. Falk Eippert)

for sub = 1:size(subjects,1)
    
    fprintf (['subject ' subjects(sub).name '--START' newline])
    
    
    for ses = 1:numel(sessions)
        
        cd(fullfile(outDir,subjects(sub).name,'func'))
        
        %calculate DVARS and REFRMS
        system(['fsl_motion_outliers ' ...
            ' -i ' sessions{ses} '_moco2.nii.gz ' ...
            ' -m ' sessions{ses} '_moco_mean_mask.nii.gz ' ...
            ' -s  dvars_' sessions{ses} '.txt ' ...
            ' -o dvars_' sessions{ses} ...
            ' --dvars --nomoco' ]);
        
        system(['fsl_motion_outliers_FALK ' ...
            ' -i ' sessions{ses} '_moco2.nii.gz ' ...
            ' -m ' sessions{ses} '_moco_mean_mask.nii.gz ' ...
            ' -s  refrms_' sessions{ses} '.txt ' ...
            ' -o  refrms_' sessions{ses} ...
            ' --refrms --nomoco' ])
        
        fprintf (['subject ' subjects(sub).name '--END' newline])
        
    end
end

% calculate outlier volumes >2 SD DVARs and REFRMS

for sub = 1:size(subjects,1)
    
    for ses = 1:numel(sessions)
        
        %where functional data is
        subPath = fullfile(outDir,subjects(sub).name, 'func');
        
        %create regressors
        REL_helper_CreateRegsOutliers(subPath,subjects(sub).name,sessions{ses},0,0);
        
    end
end

%% _   _   _   _   _   _   _   _   _   _   _   _
% / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \
%( s | e | g | m | e | n | t | a | t | i | o | n )
% \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/

%Segment the mean of images (after second step of moco) and manually
%correct them when necessary / output can be found under
%.../derivates/func/manual_moco2_mean_seg or auto_moco2_mean_seg .nii.gz

for sub = 1:size(subjects,1)
    
    cd(fullfile(outDir,subjects(sub).name,'func'))
    
    for ses = 1:numel(sessions)
        
        fprintf (['subject ' subjects(sub).name '--START' newline])
        
        system('rm *_T0*.nii.gz') % remove unnecessary files from second moco
        
        %segmentation
        system(['sct_propseg -i ' sessions{ses} '_moco2_mean.nii.gz -c t2s -radius 5 ' ...
            ' -qc ' funcsegQCdir ' -qc-subject ' subjects(sub).name])
        
        system(['fslmeants -i ' sessions{ses} '_moco2_mean ' ...
            '  -m ' sessions{ses} '_moco2_mean_seg ' ...
            ' --showall -o tmp.txt' ]);
        
        % Load data (remember that first 3 rows are coordinates)
        data     = load('tmp.txt');
        
        if length(unique(data(3,:)))== numSlices
            
            fprintf('Detected %d slices in the mask, no adjustment needed \n', numSlices)
        else
            warning ('Not enough slices [ number of slices in the mask =  %d ], the mask will be adjusted \n', length(unique(data(3,:))))
            % In the following part, if the automatic segmentation does not
            % propogate the options for sct_propseg function that affect the z-propogation
            % will be modified.
            
            % if the segmentation does not propogate, first modify the -max-area
            % that affects maximum cross-sectional area (see sct_propseg for details)
            unix(['sct_propseg -i ' sessions{ses} '_moco2_mean.nii.gz ' ...
                ' -c t2s -radius 5 -max-area 150 ' ...
                ' -qc ' funcsegQCdir ' -qc-subject ' subjects(sub).name]);
            system(['fslmeants -i ' sessions{ses} '_moco2_mean  -m ' sessions{ses} '_moco2_mean_seg  --showall -o tmp.txt' ]);
            % Load data (remember that first 3 rows are coordinates)
            data     = load('tmp.txt');
            
            if length(unique(data(3,:)))~= numSlices
                % if the segmentation still does not propogate,  modify the -max-area
                % that affects maximum cross-sectional area and min-contrast (see sct_propseg for
                % details)
                unix(['sct_propseg -i ' sessions{ses} '_moco2_mean.nii.gz ' ...
                    ' -c t2s -radius 5 -max-area 150 -min-contrast 30 '...
                    ' -qc ' funcsegQCdir ' -qc-subject ' subjects(sub).name]);
                system(['fslmeants -i ' sessions{ses} '_moco2_mean  -m ' sessions{ses} '_moco2_mean_seg  --showall -o tmp.txt' ]);
                % Load data (remember that first 3 rows are coordinates)
                data     = load('tmp.txt');
                
                if length(unique(data(3,:)))~= numSlices
                    % if the segmentation still does not propogate,  modify the -max-area
                    % that affects maximum cross-sectional area, min-contrast and max-deformation (see sct_propseg for
                    % details)
                    [status, ~] = unix(['sct_propseg -i ' sessions{ses} '_moco2_mean.nii.gz ' ...
                        ' -c t2s -radius 5 -max-area 150 -min-contrast 30 -max-deformation 5 ' ...
                        ' -qc ' funcsegQCdir ' -qc-subject ' subjects(sub).name]);
                    system(['fslmeants -i ' sessions{ses} '_moco2_mean  -m ' sessions{ses} '_moco2_mean_seg  --showall -o tmp.txt' ]);
                    % Load data (remember that first 3 rows are coordinates)
                    data     = load('tmp.txt');
                    
                end
            end
        end
    end
end

%%
%  _   _   _   _   _   _   _   _   _   _   _   _   _
% / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \ / \
%( n | o | r | m | a | l | i | z | a | t | i | o | n )
% \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/

% register PAM50 to native-space EPI images
% warp PAM50 templates to the native space

for sub = 1:size(subjects,1)
    
    cd(fullfile(outDir,subjects(sub).name,'func'))
    
    for ses = 1:numel(sessions)
        
        system(['sct_register_multimodal ' ...
            ' -i ' scttemplatepath   filesep 'PAM50' filesep 'template' filesep 'PAM50_t2.nii.gz ' ...
            ' -d ' sessions{ses}  '_moco2_mean.nii.gz ' ...
            ' -iseg ' scttemplatepath   filesep 'PAM50' filesep 'template' filesep 'PAM50_cord.nii.gz '...
            ' -dseg ' sessions{ses}  '_moco2_mean_seg.nii.gz ' ...
            ' -param step=1,type=seg,algo=centermass:step=2,type=seg,algo=bsplinesyn,metric=MeanSquares,smooth=1,slicewise=1,iter=3 ' ...
            ' -initwarp ' outDir subjects(sub).name filesep 'anat' filesep 'warps' filesep 'warp_PAM502T2.nii.gz ' ...
            ' -initwarpinv ' outDir subjects(sub).name filesep 'anat' filesep 'warps' filesep 'warp_T22PAM50.nii.gz '...
            ' -x spline '...
            ' -ofolder warps']);
        
        %make a new directory and put the warped masks etc there
        system(['sct_warp_template ' ...
            ' -d ' sessions{ses} '_moco2_mean.nii.gz ' ...
            ' -w warps' filesep 'warp_PAM50_t22' sessions{ses} '_moco2_mean.nii.gz' ...
            ' -ofolder  ' sessions{ses} '_pam50_templates' ])
        
        
    end
end
%% normalize segmental levels to the native space
cd(templateDir)
counter = 3

for lv = 1:7
    
    system(['fslmaths spinal_level_'  sprintf('%02d', counter) ...
        ' -thrp 0 ' ...
        sprintf('%02d', counter) '_thres_0prct'])
    
    system(['fslmaths '  sprintf('%02d', counter) '_thres_0prct.nii.gz ' ...
        ' -bin ' sprintf('%02d', counter) '_thres0_bin.nii.gz' ])
    
    counter = counter + 1
    
end



for sub = 1:size(subjects,1)
    counter = 3

    cd(fullfile(outDir,subjects(sub).name,'func'))
    
    for lv = 1:7
        
        for ses = 1:numel(sessions)
          
               system(['sct_apply_transfo ' ...
                ' -i ' fullfile(templateDir, ['spinal_level_' sprintf('%02d', counter) '.nii.gz ']) ...
                ' -d ' sessions{ses} '_moco2_mean.nii.gz ' ...
                ' -w warps' filesep 'warp_PAM50_t22' sessions{ses} '_moco2_mean.nii.gz' ...
                ' -o  ' sessions{ses} '_pam50_templates' filesep ['spinal_level_' sprintf('%02d', counter) '.nii.gz '] ...
                ' -x nn'])
            
            
        end
        
        counter = counter + 1
        
    end
end

%%
% threshold gray matter masks --> to ensure that there is no overlap!
masks = {'LD_unt', 'LV_unt', 'RD_unt', 'RV_unt'};

for sub = 1:size(subjects,1)  
    
    for ses = 1:numel(sessions)
        
        cd(fullfile(outDir,subjects(sub).name,'func', [sessions{ses} '_pam50_templates']))
        
        copyfile(['atlas' filesep 'PAM50_atlas_30.nii.gz'], ...
            'LV_unt.nii.gz')
        
        copyfile(['atlas' filesep 'PAM50_atlas_31.nii.gz'], ...
            'RV_unt.nii.gz')
        
        copyfile(['atlas' filesep 'PAM50_atlas_34.nii.gz'], ...
            'LD_unt.nii.gz')
        
        copyfile(['atlas' filesep 'PAM50_atlas_35.nii.gz'], ...
            'RD_unt.nii.gz')
        
        for m = 1:numel(masks)
            
            system('rm -rf tmp');
            
            mkdir('tmp')
            copyfile([masks{m} '.nii.gz'], ['tmp' filesep])
            
            cd('tmp')
            system(['fslsplit ' masks{m} ' -z'])
            
            volDir = dir('vol*.nii.gz');
            
            mkdir('tmp')
            
            
            for vx = 1:size(volDir,1)
                
                v = vx-1;
                
                [status,cmdout]   = system(['fslstats vol' sprintf('%04d',v)  ' -V']);
                cmdout = str2num(cmdout);
                numberofVoxels = cmdout(1);
                
                if isequal(numberofVoxels,0)
                    
                    subjectVoxels(sub,m,vx) = 0;
                    subjectThresholds(sub,m,vx) = NaN;
                    
                    copyfile(['vol' sprintf('%04d',v) '.nii.gz'], ['tmp' filesep ])
                    
                    
                elseif ~isequal(numberofVoxels,0)
                    
                    
                    system(['fslmaths vol' sprintf('%04d',v)  ' -thrp 70 vol' sprintf('%04d',v) '_70prct'])
                    [status,cmdout]   = system(['fslstats vol' sprintf('%04d',v)  '_70prct -V']);
                    cmdout = str2num(cmdout);
                    numberofVoxels = cmdout(1);
                    
                    subjectVoxels(sub,m,vx) = numberofVoxels;
                    subjectThresholds(sub,m,vx) = 70;
                    
                    copyfile(['vol' sprintf('%04d',v) '_70prct' '.nii.gz'], ['tmp' filesep ])
                    
                    
                    if isequal(numberofVoxels,0)
                        
                        %threshold 60 percent
                        
                        system(['fslmaths vol' sprintf('%04d',v)  ' -thrp 60 vol' sprintf('%04d',v) '_60prct'])
                        [status,cmdout]   = system(['fslstats vol' sprintf('%04d',v)  '_60prct -V']);
                        cmdout = str2num(cmdout);
                        numberofVoxels = cmdout(1);
                        
                        subjectVoxels(sub,m,vx) = numberofVoxels;
                        subjectThresholds(sub,m,vx) = 60;
                        
                        copyfile(['vol' sprintf('%04d',v) '_60prct' '.nii.gz'], ['tmp' filesep ])
                        
                        
                        if isequal(numberofVoxels,0)
                            
                            system(['fslmaths vol' sprintf('%04d',v)  ' -thrp 50 vol' sprintf('%04d',v) '_50prct'])
                            [status,cmdout]   = system(['fslstats vol' sprintf('%04d',v)  '_50prct -V']);
                            cmdout = str2num(cmdout);
                            numberofVoxels = cmdout(1);
                            
                            subjectVoxels(sub,m,vx) = numberofVoxels;
                            subjectThresholds(sub,m,vx) = 50;
                            
                            copyfile(['vol' sprintf('%04d',v) '_50prct' '.nii.gz'], ['tmp' filesep ])
                            
                            
                            
                            if isequal(numberofVoxels,0)
                                
                                system(['fslmaths vol' sprintf('%04d',v)  ' -thrp 40 vol' sprintf('%04d',v) '_40prct'])
                                [status,cmdout]   = system(['fslstats vol' sprintf('%04d',v)  '_40prct -V']);
                                cmdout = str2num(cmdout);
                                numberofVoxels = cmdout(1);
                                
                                subjectVoxels(sub,m,vx) = numberofVoxels;
                                subjectThresholds(sub,m,vx) = 40;
                                
                                copyfile(['vol' sprintf('%04d',v) '_40prct' '.nii.gz'], ['tmp' filesep ])
                                
                                
                                if isequal(numberofVoxels,0)
                                    
                                    system(['fslmaths vol' sprintf('%04d',v)  ' -thrp 30 vol' sprintf('%04d',v) '_30prct'])
                                    [status,cmdout]   = system(['fslstats vol' sprintf('%04d',v)  '_30prct -V']);
                                    cmdout = str2num(cmdout);
                                    numberofVoxels = cmdout(1);
                                    
                                    subjectVoxels(sub,m,vx) = numberofVoxels;
                                    subjectThresholds(sub,m,vx) = 30;
                                    
                                    copyfile(['vol' sprintf('%04d',v) '_30prct' '.nii.gz'], ['tmp' filesep ])
                                    
                                    
                                    if isequal(numberofVoxels,0)
                                        
                                        system(['fslmaths vol' sprintf('%04d',v)  ' -thrp 20 vol' sprintf('%04d',v) '_20prct'])
                                        [status,cmdout]   = system(['fslstats vol' sprintf('%04d',v)  '_20prct -V']);
                                        cmdout = str2num(cmdout);
                                        numberofVoxels = cmdout(1);
                                        
                                        subjectVoxels(sub,m,vx) = numberofVoxels;
                                        subjectThresholds(sub,m,vx) = 20;
                                        
                                        copyfile(['vol' sprintf('%04d',v) '_20prct' '.nii.gz'], ['tmp' filesep ])
                                        
                                        
                                        if isequal(numberofVoxels,0)
                                            
                                            system(['fslmaths vol' sprintf('%04d',v)  ' -thrp 10 vol' sprintf('%04d',v) '_10prct'])
                                            [status,cmdout]   = system(['fslstats vol' sprintf('%04d',v)  '_10prct -V']);
                                            cmdout = str2num(cmdout);
                                            numberofVoxels = cmdout(1);
                                            
                                            subjectVoxels(sub,m,vx) = numberofVoxels;
                                            subjectThresholds(sub,m,vx) = 10;
                                            
                                            copyfile(['vol' sprintf('%04d',v) '_10prct' '.nii.gz'], ['tmp' filesep ])
                                            
                                            
                                        end
                                    end
                                end
                                
                            end
                        end
                        
                    end
                end
            end
            
            cd tmp
            dirThrVol = dir('vol*');
            
            fileNames = {dirThrVol.name};
            strNames = string(fileNames);
            mergedstrNames = join(strNames," ");
            B = convertStringsToChars(mergedstrNames);
            
            system( ['fslmerge -z  ' masks{m}(1:end-3) 'thresholded ' B ]);
            
            copyfile([masks{m}(1:end-3) 'thresholded.nii.gz'], ...
                fullfile(outDir,subjects(sub).name,'func', [sessions{ses} '_pam50_templates']))
            
            cd(fullfile(outDir,subjects(sub).name, 'func', [sessions{ses} '_pam50_templates']));
            system('rm -rf tmp');
            
            system(['fslmaths ' masks{m}(1:end-3) 'thresholded.nii.gz -bin ' masks{m}(1:end-3) 'thresholded_binarized.nii.gz' ])
            
        end
    end
end

%%
% check if some slices have only ventral or dorsal masks and then remove them manually!

maskNames = {'LD_thresholded_binarized.nii.gz', 'LV_thresholded_binarized.nii.gz' ...
    'RD_thresholded_binarized.nii.gz', 'RV_thresholded_binarized.nii.gz'};

for sub = 1:numel(subjects)
    
    cd(fullfile(outDir,subjects(sub).name,'func', 'manual_pam50_templates'))
    
    LD = read_avw('LD_thresholded_binarized.nii.gz');
    LV = read_avw('LV_thresholded_binarized.nii.gz');
    RD = read_avw('RD_thresholded_binarized.nii.gz');
    RV = read_avw('RV_thresholded_binarized.nii.gz');
    
    
    for sli = 1:size(LD,3)
        
        if or(or(or(sum(sum(LD(:,:,sli))) == 0  && sum(sum(LV(:,:,sli))) ~= 0, ...
                sum(sum(RD(:,:,sli))) == 0  && sum(sum(RV(:,:,sli))) ~= 0), ...
                or(sum(sum(LD(:,:,sli))) == 0  && sum(sum(RV(:,:,sli))) ~= 0, ...
                sum(sum(RD(:,:,sli))) == 0  && sum(sum(LV(:,:,sli))) ~= 0)), ...
                or(sum(sum(LD(:,:,sli))) == 0  && sum(sum(RD(:,:,sli))) ~= 0, ...
                sum(sum(RV(:,:,sli))) == 0  && sum(sum(LV(:,:,sli))) ~= 0))
            
            maskProblem(sub,sli) = 1;
            
        end
        
    end
end
%%
% For each subject look which slice cprresponds to a certain segmental level
% and save these in a 2D matrix (subject x slices), Levels.

ses = 1; % same for each session, does not matter

for sub = 1 :size(subjects,1) 
    
    
    cd(fullfile(outDir,subjects(sub).name,'func'))
    
    counter = 3;
    
    
    for lv = 1:7
        
        
       tmp = read_avw([sessions{ses} '_pam50_templates' filesep  'spinal_level_' sprintf('%02d', counter) '.nii.gz ']);
        
        
        for sli = 1:size(tmp,3)
            
            
            if sum(sum(tmp(:,:,sli))) > 0
                
                Level_slices(lv,sli) = 1;
                voxNums(lv,sli) = sum(sum(tmp(:,:,sli)));
                
            end
            
        end
        
        
        counter = counter + 1;
        
        
    end
    
    
    for l = 1:size(Level_slices,1)
        
        C_max(l,1) = find(voxNums(l,:) == max(voxNums(l,:))) ;
        A(l,find(Level_slices(l,:))) = find(Level_slices(l,:)) ;

        
    end
    
    
    for l = 1:size(Level_slices,1)
        
        if numel(find(A(l,:))) == 1
            
            A_new(l,find(Level_slices(l,:))) = A(l,find(Level_slices(l,:)))
            
        elseif numel(find(A(l,:))) > 1
            
            A_new(l,C_max(l,1)) = C_max(l,1);
            A_new(l,(C_max(l,1)+1)) = C_max(l,1)+1;
            
            if ~isequal(C_max(l,1),1)
                
                A_new(l,(C_max(l,1)-1)) = C_max(l,1)-1;
                
            end
        end
        
            
     
    end

    overlapSlices = cell(size(Level_slices,2),1);
    nonoverlapSlices = cell(size(Level_slices,2),1);
    
    for sli = 1:size(A_new,2)
        
        [B,~] = find(A_new== sli);
        
        if numel(B)> 1
            
            overlapSlices{sli,1} = B;
            
        else
            nonoverlapSlices{sli,1} = B;
        end
        
    end
    
    for sli = 1:size(Level_slices,2)
        
        if ~isempty(overlapSlices{sli})
            
            [~,i] = max(voxNums(overlapSlices{sli}, sli));
            
            SubjectsSegLevels(sub,sli) =    (overlapSlices{sli}(i));
            
            
        elseif ~isempty((nonoverlapSlices{sli}))
            
            
            SubjectsSegLevels(sub,sli) =    (nonoverlapSlices{sli});
            
        else
            
            SubjectsSegLevels(sub,sli) = NaN;
        end
        
        
    end
    
    clear A B Level_slices voxNums C_max A_new
    
    
end
%%
[r c]=find(SubjectsSegLevels ==0);

SubjectsSegLevels_new = SubjectsSegLevels; 

for i = 1:numel(r)
    SubjectsSegLevels_new(r(i),c(i)) = NaN;
    
end

SubjectsSegLevels = SubjectsSegLevels_new;

save('SubjectsSegLevels', fullfile(outdDir,'subjectsSegmentalLevels.mat'));
%% create a 'whole' gray matter mask ---> by adding thresholded gray matter masks together for signal extraction

for sub = 1:numel(subjects)
    for z = 1:numel(sessions)
        
        cd(fullfile(outDir,subjects(sub).name,'func', [sessions{z} '_pam50_templates']))
       
        system('fslmaths LD_thresholded_binarized -add RD_thresholded_binarized Dorsal_GM ')
        system('fslmaths LV_thresholded_binarized -add RV_thresholded_binarized Ventral_GM ')
        system('fslmaths Dorsal_GM -add Ventral_GM Whole_GM')
        
    end
end