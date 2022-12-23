function [eSmooth, eSmooth_classic] = thermalDenoise_3dFWMHx(dataDir,subjects,funcDir, fileNames,...
    thermDenoise, ...
    classic, mask, maskDir, maskName, detrend, outExt,makePlot, CalculateSmoothness)

%%
% AFNI 3dFWHMx estimates the spatial smoothness of your data. For more details, see
% the function's help, please. We use this function to see how much/if thermal
% denoising increases the intrinstic spatial 'smoothness' of our data; therefore it makes
% sense to compare this to your pipeline WITHOUT thermal denoising. There
% are multiple outputs but the most interesting parameter for us is the
% effective FWHMx computed with the new 'ACF' method provided in AFNI; this
% is why I output that here in a matrix as the main output. That value can
% be also found in the ...._stdout.txt file (the last column). 
% It would mostly make sense to calculate the smoothness of residuals of the data -->
% so, your default input is possibly residuals.


%% INPUTs

% dataDir = character; main data directory where your subject specific folders are
% located

% subjects = cell array; containing the name/list of subjects

% funcDir = character; where your functional data is under subject specific
% directory

% fileNames = cell array; file names of your functional data - such as
% {'func_pain_TE40_denoised.nii.gz', 'func_pain_TE30_denoised.nii.gz'} -->
% Remember to add file extension please!!

% thermDenoise = double or logical (1 or 0);do you wanna use thermal
% denoising on your data (if you have not done so) ?

% classic= double or logical (1 or 0); do you wanna output the ClassicFWHMX parameters? --> they are NOT
% recommended anymore, they are OUTDATED, but you may wanna do so for comparison purposes

% mask = double or logical (1 or 0); do you wanna use a mask for 3dFWMx? -->
% if yes, only nonzero voxels from that func time series will be used in the calculation


% maskDir = if mask == 1; please provide a CHARACTER array to indicate the
% mask directory if it is different than dataDir/subjects/funcDir

% maskName = if mask == 1; please provide a CHARACTER array to indicate the
% mask name such as 'myExtendedCord_mask.nii.gz' --> Remember to add file extension please!!


% detrend = double or logical (1 or 0); do you wanna detrend the func
% data? recommended ONLY for NON-denoised data

% outExt = character; extension for the AFNI outputs- makes sense to give a
% unique name acc to the unique combination of parameters that you use
% such as 'te40_Cordmask_pain_' 

% makePlot =  double or logical (1 or 0); do you want to have boxplot of
% the parameters and save it to outDir with the following name
% outExt_boxPlot.png


% CalculateSmoothness = double or logical (1 or 0); do you want to run the
% run the 3dFWHMx function or just want to load the values that you
% calculated previously


%% OUTPUTs

% eSmooth = number of subjects x number of fileNames matrix of the
% smootness estimate

% eSmooth_classic = number of subjects x number of fileNames matrix of the output
% parameter --> If you set classic to 1 and would like to see old style
% smothness estimate /Note: this is outdated, use  ONLY for comparison purposes!

% a .png file --> boxplots of smoothness estimates


% Merve Kaptan, mkaptan@cbs.mpg.de
%%


for sub = 1:numel(subjects)
    
    for f = 1:numel(fileNames)
        
        
        cd(fullfile(dataDir,subjects{sub},funcDir))
        
        
        if thermDenoise
            
            thermalNoiseRemoval(fileNames{f})
            
        end
        
        
        ExtraParams = " ";
        
        if classic
            
            tmp = " -ShowMeClassicFWHM "
            ExtraParams = ExtraParams + ' ' + tmp;
            
        end
        
        
        if mask
            
            tmp = convertCharsToStrings(fullfile(maskDir, maskName));
            ExtraParams = ExtraParams + ' -mask ' + tmp;
            
        end
        
        
        if detrend
            
            tmp = " -detrend ";
            ExtraParams = ExtraParams + ' ' + tmp;
            
        end
        
        ExtraParams = convertStringsToChars(ExtraParams);
        
        outDir = fullfile(dataDir,subjects{sub},funcDir);

        
  
        
       if CalculateSmoothness 
        system(['3dFWHMx ' ...
            ' -acf ' outDir filesep fileNames{f}(1:end-7) '_' outExt '_acf '...
            ' -input '  fileNames{f}  ' '  ...
            ExtraParams ' > ' ...
            outDir filesep fileNames{f}(1:end-7) '_' outExt '_stdout.txt'])
        
       end
        tmp = [];
        tmp = load([outDir filesep fileNames{f}(1:end-7) '_' outExt '_stdout.txt']);
        
        eSmooth(sub,f) = tmp(end,end);
        
        if classic
            
            eSmooth_classic(sub,f) = tmp(1,end);
            
        else 
            
             eSmooth_classic = [];
        
        end
  
    end
end


if makePlot
    
    fiGure = figure;
    fiGure.Position = [0 0 800 800];
    hold on
    box off
    
    h = boxplot(eSmooth, ...
        'Colors', [0 0 0] ,...
        'Symbol','.r', ...
        'Widths', 0.3);
    set(h,'linew',3)
    
    xlim([0.5 numel(fileNames)+0.5])
    ylabel('FWHMx')
    xticklabels(fileNames)
    title('AFNI 3dFWMHx ACF')
    
    saveas(gcf, [outExt '_boxPlot.png'])
    
    
    if classic
        
        fiGure = figure;
        fiGure.Position = [0 0 800 800];
        hold on
        box off
        
        h = boxplot(eSmooth_classic, ...
            'Colors', [0 0 0] ,...
            'Symbol','.r', ...
            'Widths', 0.3);
        set(h,'linew',3)
        
        xlim([0.5 numel(fileNames)+0.5])
        ylabel('FWHMx')
        xticklabels(fileNames)
        title('AFNI 3dFWMHx Classic')
        
        saveas(gcf, [outExt '_classic_boxPlot.png'])
        
    end  
end

end