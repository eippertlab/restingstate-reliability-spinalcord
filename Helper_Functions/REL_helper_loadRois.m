function [tmpLD,tmpRD, tmpLV, tmpRV, tmpTSNR, CommonSlices] = REL_helper_loadRois(method,sessionId, fileName, fileType, maskType)

if isequal(fileType, 'tsnr')
    
    
    if ~isequal(maskType, 'whole') 
    
    tmpLD = load([sessionId '_' fileName '_LD_thresholded_binarized.txt']);
    tmpLD = tmpLD';
    
    
    tmpRD = load([sessionId '_' fileName '_RD_thresholded_binarized.txt']);
    tmpRD = tmpRD';
    
    
    tmpLV = load([sessionId '_' fileName  '_LV_thresholded_binarized.txt']);
    tmpLV = tmpLV';
    
    
    tmpRV = load([sessionId '_' fileName  '_RV_thresholded_binarized.txt']);
    tmpRV = tmpRV';
    tmpTSNR = [];
    
    else
        
        tmpTSNR = load([sessionId '_' fileName  '_Whole_GM.txt']);
        tmpLD = []
        tmpRD = []
        tmpLV = []
        tmpRV = [];
        
    end
    
    
    if isequal(method,'slicewise')
        
        
        CommonSlices = intersect(intersect(intersect(unique(tmpLD(:,3)),unique(tmpLV(:,3)), 'stable'),unique(tmpRD(:,3)), 'stable'),unique(tmpRV(:,3)), 'stable');
        CommonSlices2 = unique(tmpLD(:,3));
        
        if ~isequal(CommonSlices, CommonSlices2)
            
            
            error('Probably mask error!! check your mask!!')
            
            
        end
        
        
        
    end
    
    
elseif isequal(fileType, 'timeSeries')
    
    tmpLD = load([sessionId '_' fileName '_LD_thresholded_binarized.txt']);
    tmpLD = tmpLD';
    
    
    tmpRD = load([sessionId '_' fileName '_RD_thresholded_binarized.txt']);
    tmpRD = tmpRD';
    
    
    tmpLV = load([sessionId '_' fileName  '_LV_thresholded_binarized.txt']);
    tmpLV = tmpLV';
    
    
    tmpRV = load([sessionId '_' fileName  '_RV_thresholded_binarized.txt']);
    tmpRV = tmpRV';
    
    
    CommonSlices = intersect(intersect(intersect(unique(tmpLD(:,3)),unique(tmpLV(:,3)), 'stable'),unique(tmpRD(:,3)), 'stable'),unique(tmpRV(:,3)), 'stable');
    CommonSlices2 = unique(tmpLD(:,3));
    
    if ~isequal(CommonSlices, CommonSlices2)
        
        
        error('Probably mask error!! check your mask!!')
        
        
    end
    
elseif  isequal(fileType, 'expvar')
    
    
    if ~isequal(maskType, 'whole')
        
        tmpLD = load([sessionId '_' fileName '_LD_thresholded_binarized.txt']);
        tmpLD = tmpLD';
        
        
        tmpRD = load([sessionId '_' fileName '_RD_thresholded_binarized.txt']);
        tmpRD = tmpRD';
        
        
        tmpLV = load([sessionId '_' fileName  '_LV_thresholded_binarized.txt']);
        tmpLV = tmpLV';
        
        
        tmpRV = load([sessionId '_' fileName  '_RV_thresholded_binarized.txt']);
        tmpRV = tmpRV';
        
    else
        
        tmpTSNR = load([sessionId '_' fileName  '_Whole_GM.txt']);
        tmpLD = []
        tmpRD = []
        tmpLV = []
        tmpRV = [];
        
    end
    
    
    
    
end


if isequal(method, 'avgTS')
    
    tmpLD = mean(tmpLD(ismember(tmpLD(:,3),CommonSlices),4:end));
    tmpRD = mean(tmpRD(ismember(tmpRD(:,3),CommonSlices),4:end));
    
    
    tmpLV = mean(tmpLV(ismember(tmpLV(:,3),CommonSlices),4:end));
    tmpRV = mean(tmpRV(ismember(tmpRV(:,3),CommonSlices),4:end));
    
elseif isequal(method, 'slice-wise')
    
    
    tmpLD = tmpLD(ismember(tmpLD(:,3),CommonSlices),3:end);
    tmpRD = tmpRD(ismember(tmpRD(:,3),CommonSlices),3:end);
    
    tmpLV = tmpLV(ismember(tmpLV(:,3),CommonSlices),3:end);
    tmpRV = tmpRV(ismember(tmpRV(:,3),CommonSlices),3:end);
    
    
end


end

