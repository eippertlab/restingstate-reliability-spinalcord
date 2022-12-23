function REL_helper_CalculatetSNR(outDir, subjects,fileNames, sessions)

parfor sub = 1:size(subjects,1)
    
    for ses = 1:numel(sessions)
        
        for f = 1:numel(fileNames)
            
            cd(fullfile(outDir, subjects(sub).name, 'func', [sessions{ses} '_' fileNames{f} '.feat']))

                system(['fslmaths stats' filesep 'res4d -Tstd stats' filesep 'res4d_sd'  ])
                system(['fslmaths mean_func -div stats' filesep 'res4d_sd stats' filesep 'tsnr.nii.gz' ])
                
                
        end
    end
    
end


