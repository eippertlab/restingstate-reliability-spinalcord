function REL_helper_CreateRegsOutliers(dataPath,subjectName,sessionName,createFigures,saveFigures)

%%%% INPUTs
% dataPath     --> where functional files/runs are; give full path
% subjectName  --> name of the subject
% sessionName  --> name of the session
% saveFigures  --> save the output figure or not / saved as .fig
%%

cd(dataPath)

fprintf (['subject ' subjectName '--START' newline])

% load raw metric data
metricDVARS  = load(['dvars_' sessionName '.txt']);
metricREFRMS = load(['refrmsFALK_' sessionName '.txt']);
metricREFRMS = metricREFRMS(2:end);

% calculate own thresholds
tmpThreshDeDVARS  = prctile(metricDVARS ,75) + iqr(metricDVARS);
tmpThreshZ2DVARS  = mean(metricDVARS) + 2*std(metricDVARS);
tmpThreshZ3DVARS  = mean(metricDVARS) + 3*std(metricDVARS);

tmpThreshDeREFRMS = prctile(metricREFRMS,75) + iqr(metricREFRMS);
tmpThreshZ2REFRMS = mean(metricREFRMS) + 2*std(metricREFRMS);
tmpThreshZ3REFRMS = mean(metricREFRMS) + 3*std(metricREFRMS);


% get number of outliers according to own thresholds
threshZ2DVARS = find(metricDVARS >= tmpThreshZ2DVARS);
threshZ2REFRMS = find(metricREFRMS >= tmpThreshZ2REFRMS);

if createFigures
    % plot data with thresholds
    figure;
    
    subplot(2,1,1); title(['Subject: ' subjectName ', DVARS ' sessionName], 'Interpreter', 'none'); ...
        hold on; plot(metricDVARS, 'k');
    plot(1:numel(metricDVARS), tmpThreshZ2DVARS *ones(1,numel(metricDVARS)), 'b--'); ...
        plot(1:numel(metricDVARS), tmpThreshZ3DVARS *ones(1,numel(metricDVARS)), 'b');
    plot(1:numel(metricDVARS), tmpThreshDeDVARS *ones(1,numel(metricDVARS)), 'g');
    
    subplot(2,1,2); title(['Subject: ' subjectName ', REFRMS ' sessionName]); hold on; plot(metricREFRMS, 'k');
    plot(1:numel(metricREFRMS), tmpThreshZ2REFRMS*ones(1,numel(metricREFRMS)), 'b--'); ...
        plot(1:numel(metricREFRMS), tmpThreshZ3REFRMS*ones(1,numel(metricREFRMS)), 'b');
    plot(1:numel(metricREFRMS), tmpThreshDeREFRMS *ones(1,numel(metricREFRMS)), 'g');
    
    pause(3)
    
    close all
    
end


% get number of outliers according to default threshold of 2SD
numDVARSOutliers      = numel([threshZ2DVARS]);
indicesDVARSOutliers  = [threshZ2DVARS];
numREFRMSOutliers     = numel([threshZ2REFRMS]);
indicesREFRMSOutliers = [threshZ2REFRMS];
numTotalOutliers      = numel([threshZ2DVARS; threshZ2REFRMS]);
numUniqueOutliers     = numel(unique([threshZ2DVARS; threshZ2REFRMS]));
indicesOutliers       = unique([threshZ2DVARS; threshZ2REFRMS]);

%
if isempty(numUniqueOutliers) & isempty(indicesOutliers) ...
        & isempty(numTotalOutliers)
    %
    disp(['******* No outliers for this subject: ' name '&' sessionName ' *******']);
    %
else
    %
    % SAVE FILES
    
    save([sessionName '_moco_outliers.mat'], ...
        'numDVARSOutliers', 'indicesDVARSOutliers', ...
        'numREFRMSOutliers','indicesREFRMSOutliers', ....
        'indicesOutliers')
    
    if ~isempty(indicesOutliers)
        
        outlier = zeros(numel(metricDVARS),size(indicesOutliers,1));
        
        for c = 1:size(indicesOutliers,1)
            
            outlier(indicesOutliers(c),c) = 1;
            
        end
        
    elseif isempty(indicesOutliers)
        
        outlier = zeros(numel(metricDVARS),1);
        
    end
    
    writematrix(outlier,[sessionName '_moco_outlies.txt'],'Delimiter', 'space')
    
    
    clear indicesOutliers
    
    if saveFigures
        savefig(gcf, [sessionName '_moco_outliers'] )
    end
    
end


fprintf (['subject ' subjectName '--END' newline])


% %%%%%%%%%%

end