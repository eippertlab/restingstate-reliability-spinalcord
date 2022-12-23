function REL_helper_ExtractTimeSeries(outDir, subjects,niftiName,fileNames, masks, sessions, txtOutDir, showAll)

for sub = 1:size(subjects,1)
    
    if ~exist(fullfile(outDir, subjects(sub).name, 'func', txtOutDir))
        
        mkdir(fullfile(outDir, subjects(sub).name, 'func', txtOutDir))
        
    end
    
    for ses = 1:numel(sessions)
        
        
        if contains(sessions{ses},'auto')
            
            maskDir = 'auto_pam50_templates';
            
        elseif contains(sessions{ses}, 'manual')
            
            maskDir = 'manual_pam50_templates';
            
        end
        
        for f = 1:numel(fileNames)
            
            cd(fullfile(outDir, subjects(sub).name, 'func', [sessions{ses} '_' fileNames{f} '.feat'],'stats'))
            
            for m = 1:numel(masks)
                
%                 
                if showAll 
                    
                 system(['fslmeants -i ' niftiName ...
                        ' -m ' fullfile(outDir, subjects(sub).name, 'func', maskDir, masks{m}) ...
                        ' -o ' fullfile(outDir, subjects(sub).name, 'func', txtOutDir, [sessions{ses} '_' fileNames{f} '_' masks{m} '.txt'])  ...
                        ' --showall'])
                  
                else
%                     
                    system(['fslmeants -i ' niftiName ...
                        ' -m ' fullfile(outDir, subjects(sub).name, 'func', maskDir, masks{m}) ...
                        ' -o ' fullfile(outDir, subjects(sub).name, 'func', txtOutDir, [sessions{ses} '_' fileNames{f} '_' masks{m} '.txt'])])
                

                end
            end
        end
    end
    
end