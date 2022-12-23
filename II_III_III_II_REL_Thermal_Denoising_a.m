% In this code, the thermal noise removal step which is listed for the following
% manuscript:  "" section: 2.3.3.2 Thermal noise 

% In order to run this code the Spinal Cord Toolbox (version 4.2.2 or higher -
% note that for the manuscript version 4.2.2 was used)
% and FSL should be added to the bash profile of the user to call their
% functions from MATLAB.


% Perform the thermal noise removal see:
% https://github.com/NYU-DiffusionMRI/mppca_denoise on RAW .nii.gz data
% and motion-correct the data afterwards.

%%
clc
clear all
close all
%% set the paths and directories

inDir  = '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/rawdata/'
outDir = '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/'
addpath(genpath('/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/Code/'))

% add fsl to the path
fslDir = '/afs/cbs.mpg.de/software/fsl/5.0.11/ubuntu-xenial-amd64/';
setenv('FSLDIR', fslDir);
pathFSL = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',pathFSL);
path(path, fsldirmpath);
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be

cd(inDir)
subjects = dir('sub-ZS*');

subjects([9,18,30]) = []; %exclude subjects with ECG problems
sessions  = {'manual', 'auto'};

%% thermal noise removal 
% see: https://github.com/NYU-DiffusionMRI/mppca_denoise
for sub = 1:numel(subjects)
    for ses = 1:numel(sessions)
        
        if doCalculations
            cd(fullfile(inDir,subjects(sub).name,'func'))
            
            system(['cp *' sessions{ses} '*.nii.gz '  ...
                fullfile(outDir,subjects(sub).name,'func', [sessions{ses} '_raw.nii.gz']) ]);
            thermalNoiseRemoval(sessions{ses},[sessions{ses} '_raw.nii.gz'])
        end
    end
end
%% Motion-correction of thermal denoised data
% Motion-correction STEP1
for sub = 1:size(subjects,1)
    
    for ses = 1:numel(sessions)
        
        cd(fullfile(inDir,subjects(sub).name,'func'))
        
        
        fprintf (['subject ' subjects(sub).name '--START'])
        
        mkdir(fullfile(outDir,subjects(sub).name,'func'))
        
        system(['cp ' sessions{ses} '_denoised.nii.gz '  ...
            fullfile(outDir,subjects(sub).name,'func', [sessions{ses} '_denoised.nii.gz']) ]);
        
        cd(fullfile(outDir,subjects(sub).name,'func'))
        
        system(['fslmaths ' sessions{ses} '_denoised.nii.gz -Tmean '   sessions{ses} '_denoised_mean.nii.gz' ]);
        
        system(['sct_propseg -i ' sessions{ses} '_denoised_mean.nii.gz -c t2s' ])
        
        system(['fslmeants -i ' sessions{ses} '_denoised_mean  -m ' sessions{ses} '_denoised_mean_seg  --showall -o tmp.txt' ]);
        
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
            unix(['sct_propseg -i ' sessions{ses} '_denoised_mean.nii.gz  -c t2s -radius 5 -max-area 150']);
            system(['fslmeants -i ' sessions{ses} '_denoised_mean  -m ' sessions{ses} '_denoised_mean_seg  --showall -o tmp.txt' ]);
            % Load data (remember that first 3 rows are coordinates)
            data     = load('tmp.txt');
            
            if length(unique(data(3,:)))~= numSlices
                % if the segmentation still does not propogate,  modify the -max-area
                % that affects maximum cross-sectional area and min-contrast (see sct_propseg for
                % details)
                unix(['sct_propseg -i ' sessions{ses} '_denoised_mean.nii.gz -c t2s -radius 5 -max-area 150 -min-contrast 30']);
                system(['fslmeants -i ' sessions{ses} '_denoised_mean  -m ' sessions{ses} '_denoised_mean_seg  --showall -o tmp.txt' ]);
                % Load data (remember that first 3 rows are coordinates)
                data     = load('tmp.txt');
                
                if length(unique(data(3,:)))~= numSlices
                    % if the segmentation still does not propogate,  modify the -max-area
                    % that affects maximum cross-sectional area, min-contrast and max-deformation (see sct_propseg for
                    % details)
                    [status, ~] = unix(['sct_propseg -i ' sessions{ses} '_denoised_mean.nii.gz -c t2s -radius 5 -max-area 150 -min-contrast 30 -max-deformation 5']);
                    system(['fslmeants -i ' sessions{ses} '_denoised_mean  -m ' sessions{ses} '_denoised_mean_seg  --showall -o tmp.txt' ]);
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
            
            system(['sct_create_mask -i ' sessions{ses} '_denoised_mean.nii.gz -p centerline,' sessions{ses} '_denoised_mean_seg.nii.gz -o ' sessions{ses} '_denoised_mean_mask.nii.gz' ])
            system(['fslmerge -tr ' sessions{ses} '_denoised_mean_merged '  sessions{ses} '_denoised_mean '  sessions{ses}  '_denoised 2.312 '  ])
            
            system(['sct_fmri_moco ' ...
                '-i '  sessions{ses} '_denoised_mean_merged.nii.gz ' ...
                '-m ' sessions{ses} '_denoised_mean_mask.nii.gz ' ...
                '-param iterAvg=0 -x spline']);
            
            system(['sct_image ' ...
                ' -i '  sessions{ses} '_denoised_mean_merged_moco.nii.gz ' ...
                ' -remove-vol 0 ' ...
                ' -o ' sessions{ses} '_denoised_moco.nii.gz']);
            
            
            %                 copyfile('moco_params_x.nii.gz', [ sessions{ses} '_moco_params_x.nii.gz'])
            %                 copyfile('moco_params_y.nii.gz', [ sessions{ses} '_moco_params_y.nii.gz'])
            %
            %                 copyfile('moco_params.tsv', [ sessions{ses} '_moco_params.tsv'])
            %
            system(['fslmaths '  sessions{ses} '_denoised_moco.nii.gz' ...
                ' -Tmean ' ...
                sessions{ses} '_denoised_moco_mean '])
            
        else
            mocoProblemSubject{sub,ses} = subjects(sub).name;
            
            
        end
        
        fprintf (['subject ' subjects(sub).name '--DONE'])
        
        
    end
    
end
%% MOTION-CORRECTION: Step 2

for sub = 1:size(subjects,1)
    
    for ses = 1:numel(sessions)
        
        
        fprintf (['subject ' subjects(sub).name '--START' newline])
        
        cd(fullfile(outDir,subjects(sub).name,'func'))
        
        system('rm *_T0*.nii.gz')
        
        system(['fslmaths ' sessions{ses} '_denoised_moco.nii.gz ' ...
            ' -Tmean '   sessions{ses} '_denoised_moco_mean.nii.gz' ]);
        
        system(['sct_propseg -i ' sessions{ses} '_denoised_moco_mean.nii.gz -c t2s' ])
        
        system(['fslmeants -i ' sessions{ses} '_denoised_moco_mean ' ...
            '  -m ' sessions{ses} '_denoised_moco_mean_seg ' ...
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
            unix(['sct_propseg -i ' sessions{ses} '_denoised_moco_mean.nii.gz  -c t2s -radius 5 -max-area 150']);
            system(['fslmeants -i ' sessions{ses} '_denoised_moco_mean  -m ' sessions{ses} '_denoised_moco_mean_seg  --showall -o tmp.txt' ]);
            % Load data (remember that first 3 rows are coordinates)
            data     = load('tmp.txt');
            
            if length(unique(data(3,:)))~= numSlices
                % if the segmentation still does not propogate,  modify the -max-area
                % that affects maximum cross-sectional area and min-contrast (see sct_propseg for
                % details)
                unix(['sct_propseg -i ' sessions{ses} '_denoised_moco_mean.nii.gz -c t2s -radius 5 -max-area 150 -min-contrast 30']);
                system(['fslmeants -i ' sessions{ses} '_denoised_moco_mean  -m ' sessions{ses} '_denoised_moco_mean_seg  --showall -o tmp.txt' ]);
                % Load data (remember that first 3 rows are coordinates)
                data     = load('tmp.txt');
                
                if length(unique(data(3,:)))~= numSlices
                    % if the segmentation still does not propogate,  modify the -max-area
                    % that affects maximum cross-sectional area, min-contrast and max-deformation (see sct_propseg for
                    % details)
                    [status, ~] = unix(['sct_propseg -i ' sessions{ses} '_denoised_moco_mean.nii.gz -c t2s -radius 5 -max-area 150 -min-contrast 30 -max-deformation 5']);
                    system(['fslmeants -i ' sessions{ses} '_denoised_moco_mean  -m ' sessions{ses} '_denoised_moco_mean_seg  --showall -o tmp.txt' ]);
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
            
            system(['sct_create_mask ' ...
                ' -i ' sessions{ses} '_denoised_moco_mean.nii.gz ' ...
                ' -p centerline,' sessions{ses} '_denoised_moco_mean_seg.nii.gz '...
                ' -o ' sessions{ses} '_denoised_moco_mean_mask.nii.gz' ])
            
            system(['fslmerge -tr ' sessions{ses} '_denoised_moco_mean_merged ' ...
                sessions{ses} '_denoised_moco_mean '  sessions{ses}  '_denoised 2.312 '  ])
            
            system(['sct_fmri_moco ' ...
                '-i '  sessions{ses} '_denoised_moco_mean_merged.nii.gz ' ...
                '-m ' sessions{ses} '_denoised_moco_mean_mask.nii.gz ' ...
                '-param iterAvg=0 -x spline']);
            
            system(['sct_image ' ...
                ' -i '  sessions{ses} '_denoised_moco_mean_merged_moco.nii.gz ' ...
                ' -remove-vol 0 ' ...
                ' -o ' sessions{ses} '_denoised_moco2.nii.gz']);
            
            
            %                 copyfile('moco_params_x.nii.gz', [ sessions{ses} '_moco_params_x.nii.gz'])
            %                 copyfile('moco_params_y.nii.gz', [ sessions{ses} '_moco_params_y.nii.gz'])
            %
            %                 copyfile('moco_params.tsv', [ sessions{ses} '_moco_params.tsv'])
            %
            system(['fslmaths '  sessions{ses} '_denoised_moco2.nii.gz' ...
                ' -Tmean ' ...
                sessions{ses} '_denoised_moco2_mean '])
            
        else
            mocoProblemSubject{sub,ses} = subjects(sub).name;
            
            
        end
        
        fprintf (['subject ' subjects(sub).name '--DONE'])
        
    end
end

%% NOTE
% Remember: for subject sub-ZS003 auto session the recording of triggers
% started a bit later than the functional acquisition. Remove the first 9
% volumes for which we do not have BrainAmp triggers. Also remember: in FSL,
% counting starts from 0 - not from 1.

cd(fullfile(outDir, 'sub-ZS003', 'func'))
system('fslroi auto_denoised_moco2_paramX.nii.gz auto_denoised_moco2_paramX.nii.gz 0 -1 0 -1 0 -1 9 -1')
system('fslroi auto_denoised_moco2_paramY.nii.gz auto_denoised_moco2_paramY.nii.gz 0 -1 0 -1 0 -1 9 -1')
system('fslroi auto_denoised_moco2.nii.gz auto_denoised_moco2.nii.gz 0 -1 0 -1 0 -1 9 -1')
