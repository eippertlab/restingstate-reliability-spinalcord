
function REL_helper_CalculateExpvar(outDir, subjects,fileNames, sessions)

for sub = 1:size(subjects,1)
         
    for ses = 1:numel(sessions)
        
        if ~exist(fullfile(outDir, subjects(sub).name, 'func', [sessions{ses} '_min.feat'], 'stats', 'res4d_sqr.nii.gz'))
            
            cd(fullfile(outDir, subjects(sub).name, 'func', [sessions{ses} '_min.feat'], 'stats'))
            
            system(['fslmaths res4d_sd -sqr  res4d_sqr'])
   
        end
          
        
            cd(fullfile(outDir, subjects(sub).name, 'func'))

        system(['fslmaths ' sessions{ses}  '_moco2_mean  -thr 10000000 thresholdedImg_0'])
        system('fslmaths thresholdedImg_0  -add 1 allOnesImg')
        
        
        
        for f = 1:numel(fileNames)
            
            cd(fullfile(outDir, subjects(sub).name, 'func', [sessions{ses} '_' fileNames{f} '.feat'], 'stats'))
            
            system(['fslmaths res4d_sd -sqr  res4d_sqr'])

            system(['fslmaths res4d_sqr -div ' ...
               fullfile(outDir, subjects(sub).name, 'func', [sessions{ses} '_min.feat'], 'stats', 'res4d_sqr') ...
                ' var_ratio' ])
            
            system(['fslmaths '   fullfile(outDir, subjects(sub).name, 'func','allOnesImg.nii.gz' ) ...
                ' -sub var_ratio ' ...
                ' expvar' ])
            
            
            
            
        end
        
        
        
        
    end
    
end

    
end