
function [DD_session1, DD_session2, VV_session1, VV_session2,LDV_session1,LDV_session2,RDV_session1,RDV_session2, ...
    within_session1, within_session2, between_session1, between_session2, ICCvalue, ICCvalue_CI ] = REL_helper_CalculateConnandRel_Segmental(method,dataDir,subjects,analysisSpace, sessions, CorrType, roiTS,fileName, signalFolder ,CalculateICC, Levels)

levS = 7:-1:1;

eval(['DD_' sessions{1} '= nan(numel(subjects),24);']);
eval(['DD_' sessions{2} '= nan(numel(subjects),24);']);

eval(['VV_' sessions{1} '= nan(numel(subjects),24);']);
eval(['VV_' sessions{2} '= nan(numel(subjects),24);']);


eval(['LDV_' sessions{1} '= nan(numel(subjects),24);']);
eval(['LDV_' sessions{2} '= nan(numel(subjects),24);']);

eval(['RDV_' sessions{1} '= nan(numel(subjects),24);']);
eval(['RDV_' sessions{2} '= nan(numel(subjects),24);']);


eval(['LDRV_' sessions{1} '= nan(numel(subjects),24);']);
eval(['LDRV_' sessions{2} '= nan(numel(subjects),24);']);

eval(['RDLV_' sessions{1} '= nan(numel(subjects),24);']);
eval(['RDLV_' sessions{2} '= nan(numel(subjects),24);']);


for s = 1:numel(subjects)
    
    for ses = 1:numel(sessions)
        
        cd(fullfile(dataDir, subjects{s}, 'func',signalFolder ))
        
        sessionId = sessions{ses};
        
        if isequal(method, 'avgTS')
            
            clear tmpLD tmpRD tmpLV tmpRV CommonSlices LD RD LV RV
            
            [tmpLD,tmpRD, tmpLV, tmpRV, CommonSlices] = REL_helper_loadRois('slicewise',sessionId, fileName, 'timeSeries');
            
            CommonSlices = CommonSlices +1;
            
            subjectLevels = Levels(s,:);
            
            subjectLevels = subjectLevels(CommonSlices);
            
            tmpLD(:,3) = tmpLD(:,3)+1;
            tmpRD(:,3) = tmpRD(:,3)+1;
            
            tmpLV(:,3) = tmpLV(:,3)+1;
            tmpRV(:,3) = tmpRV(:,3)+1;
            
            for sli = 1:numel(CommonSlices)
                
                
                LD(sli,:) = mean(tmpLD(find(tmpLD(:,3)==CommonSlices(sli)),4:end),1);
                RD(sli,:) = mean(tmpRD(find(tmpRD(:,3)==CommonSlices(sli)),4:end),1);
                
                LV(sli,:) = mean(tmpLV(find(tmpLV(:,3)==CommonSlices(sli)),4:end),1);
                RV(sli,:) = mean(tmpRV(find(tmpRV(:,3)==CommonSlices(sli)),4:end),1);
                
                
            end
            
            
            for l = 1:numel(levS)
                clear LD_l RD_l LV_l RV_l
                
                LD_l = mean(LD(subjectLevels==levS(l),:),1);
                RD_l = mean(RD(subjectLevels==levS(l),:),1);
                
                LV_l = mean(LV(subjectLevels==levS(l),:),1);
                RV_l = mean(RV(subjectLevels==levS(l),:),1);
                
                if isequal(CorrType, 'pearson')
                    
                    x = [LD_l',RD_l'];
                    rho = corr(x);
                    eval(['DD_' sessions{ses}  '_lv(s,l) = rho(1,2);']);
                    
                    clear x rho
                    x = [LV_l',RV_l'];
                    rho = corr(x);
                    eval(['VV_' sessions{ses}  '_lv(s,l) = rho(1,2);']);
                    
                    clear x rho
                    x = [LD_l',LV_l'];
                    rho = corr(x);
                    eval(['LDV_' sessions{ses}  '_lv(s,l) = rho(1,2);']);
                    
                    clear x rho
                    x = [RD_l',RV_l'];
                    rho = corr(x);
                    eval(['RDV_' sessions{ses}  '_lv(s,l) = rho(1,2);']);
                    
                    clear x rho
                    x = [LD_l',RV_l'];
                    rho = corr(x);
                    eval(['LDRV_' sessions{ses}  '_lv(s,l) = rho(1,2);']);
                    
                    clear x rho
                    x = [RD_l',LV_l'];
                    rho = corr(x);
                    eval(['RDLV_' sessions{ses}  '_lv(s,l) = rho(1,2);']);
                    
                elseif isequal(CorrType,'partial')
                    x = [LD_l',RD_l'];
                    z = [LV_l', RV_l'];
                    rho = partialcorr(x,z);
                    eval(['DD_' sessions{ses}  '_lv(s,l) = rho(1,2);']);
                    clear x z rho
                    
                    x = [LV_l',RV_l'];
                    z = [LD_l', RD_l'];
                    rho = partialcorr(x,z);
                    eval(['VV_' sessions{ses}  '_lv(s,l) = rho(1,2);']);
                    
                    clear x rho
                    x = [LD_l',LV_l'];
                    z = [RD_l',RV_l'];
                    rho = partialcorr(x,z);
                    eval(['LDV_' sessions{ses}  '_lv(s,l) = rho(1,2);']);
                    
                    clear x rho
                    x = [RD_l',RV_l'];
                    z = [LD_l',LV_l'];
                    rho = partialcorr(x,z);
                    eval(['RDV_' sessions{ses}  '_lv(s,l) = rho(1,2);']);
                    
                    clear x rho
                    x = [LD_l',RV_l'];
                    z = [RD_l', LV_l'];
                    rho = partialcorr(x,z);
                    eval(['LDRV_' sessions{ses}  '_lv(s,l) = rho(1,2);']);
                    
                    clear x rho
                    x = [RD_l',LV_l'];
                    z = [LD_l', RV_l'];
                    rho = partialcorr(x,z);
                    eval(['RDLV_' sessions{ses}  '_lv(s,l) = rho(1,2);']);
                end
                
                
            end
            
            
            
        elseif isequal(method, 'slicewise')
            
            [tmpLD,tmpRD, tmpLV, tmpRV, CommonSlices] = REL_helper_loadRois(method,sessionId, fileName, 'timeSeries');
            CommonSlices = CommonSlices +1;
            
            subjectLevels = Levels(s,:);
            
            subjectLevels = subjectLevels(CommonSlices);
            
            tmpLD(:,3) = tmpLD(:,3)+1;
            tmpRD(:,3) = tmpRD(:,3)+1;
            
            tmpLV(:,3) = tmpLV(:,3)+1;
            tmpRV(:,3) = tmpRV(:,3)+1;
            
            
            for sli = 1:numel(CommonSlices)
                
                if isequal(roiTS,'mean')
                    
                    LD = mean(tmpLD(find(tmpLD(:,3)==CommonSlices(sli)),4:end),1);
                    RD = mean(tmpRD(find(tmpRD(:,3)==CommonSlices(sli)),4:end),1);
                    
                    LV = mean(tmpLV(find(tmpLV(:,3)==CommonSlices(sli)),4:end),1);
                    RV = mean(tmpRV(find(tmpRV(:,3)==CommonSlices(sli)),4:end),1);
                    
                elseif isequal(roiTS,'voxel')
                    
                    LD = tmpLD(find(tmpLD(:,3)==CommonSlices(sli)),4:end);
                    RD = tmpRD(find(tmpRD(:,3)==CommonSlices(sli)),4:end);
                    
                    LV = mean(tmpLV(find(tmpLV(:,3)==CommonSlices(sli)),4:end),1);
                    RV = mean(tmpRV(find(tmpRV(:,3)==CommonSlices(sli)),4:end),1);
                    
                end
                
                if isequal(CorrType, 'pearson')
                    x = [LD',RD'];
                    rho = corr(x);
                    
                    if isequal(roiTS,'mean')
                        eval(['DD_' sessions{ses}  '(s,sli) = rho(1,2);']);
                    elseif isequal(roiTS,'voxel')
                        tmp_rho = rho( 1:size(LD,1), end-(size(RD,1)-1):end);
                        eval(['DD_' sessions{ses}  '(s,sli) = prctile(tmp_rho(:),95);']);
                    end
                    
                elseif isequal(CorrType,'partial')
                    x = [LD',RD'];
                    z = [LV', RV'];
                    rho = partialcorr(x,z);
                    
                    if isequal(roiTS, 'mean')
                        eval(['DD_' sessions{ses}  '(s,sli) = rho(1,2);']);
                        
                    elseif isequal(roiTS,'voxel')
                        tmp_rho = rho( 1:size(LD,1), end-(size(RD,1)-1):end);
                        eval(['DD_' sessions{ses}  '(s,sli) = prctile(tmp_rho(:),95);']);
                    end
                    
                end
                clear x z rho tmp_rho
                clear LD RD LV RV
                
            end
            
            for l = 1:numel(levS)
                
                eval(['DD_' sessions{ses} '_lv(s,l) = '  'mean(DD_' sessions{ses} '(s,subjectLevels==levS(l)), ''omitnan'');'])
                
            end
            
            
            
            for sli = 1:numel(CommonSlices)
                
                if isequal(roiTS,'mean')
                    
                    LV = mean(tmpLV(find(tmpLV(:,3)==CommonSlices(sli)),4:end),1);
                    RV = mean(tmpRV(find(tmpRV(:,3)==CommonSlices(sli)),4:end),1);
                    
                    LD = mean(tmpLD(find(tmpLD(:,3)==CommonSlices(sli)),4:end),1);
                    RD = mean(tmpRD(find(tmpRD(:,3)==CommonSlices(sli)),4:end),1);
                    
                elseif isequal(roiTS,'voxel')
                    
                    LV = tmpLV(find(tmpLV(:,3)==CommonSlices(sli)),4:end);
                    RV = tmpRV(find(tmpRV(:,3)==CommonSlices(sli)),4:end);
                    
                    LD = mean(tmpLD(find(tmpLD(:,3)==CommonSlices(sli)),4:end),1);
                    RD = mean(tmpRD(find(tmpRD(:,3)==CommonSlices(sli)),4:end),1);
                    
                end
                
                
                if isequal(CorrType, 'pearson')
                    
                    x = [LV',RV'];
                    rho = corr(x);
                    
                    if isequal(roiTS,'mean')
                        eval(['VV_' sessions{ses}  '(s,sli) = rho(1,2);']);
                    elseif isequal(roiTS, 'voxel')
                        tmp_rho = rho( 1:size(LV,1), end-(size(RV,1)-1):end);
                        eval(['VV_' sessions{ses}  '(s,sli) = prctile(tmp_rho(:),95);']);
                    end
                    
                elseif isequal(CorrType,'partial')
                    x = [LV',RV'];
                    z = [LD', RD'];
                    rho = partialcorr(x,z);
                    
                    if isequal(roiTS, 'mean')
                        eval(['VV_' sessions{ses}  '(s,sli) = rho(1,2);']);
                    elseif isequal(roiTS,'voxel')
                        tmp_rho = rho( 1:size(LV,1), end-(size(RV,1)-1):end);
                        eval(['VV_' sessions{ses}  '(s,sli) =prctile(tmp_rho(:),95);']);
                    end
                end
                clear x z rho tmp_rho
                clear LD RD LV RV
                
            end
            
            
            
            for l = 1:numel(levS)
                
                eval(['VV_' sessions{ses} '_lv(s,l) = '  'mean(VV_' sessions{ses} '(s,subjectLevels==levS(l)), ''omitnan'');'])
                
            end
            
            for sli = 1:numel(CommonSlices)
                
                if isequal(roiTS,'mean')
                    
                    LD = mean(tmpLD(find(tmpLD(:,3)==CommonSlices(sli)),4:end),1);
                    LV = mean(tmpLV(find(tmpLV(:,3)==CommonSlices(sli)),4:end),1);
                    
                    
                    RD = mean(tmpRD(find(tmpRD(:,3)==CommonSlices(sli)),4:end),1);
                    RV = mean(tmpRV(find(tmpRV(:,3)==CommonSlices(sli)),4:end),1);
                    
                elseif isequal(roiTS,'voxel')
                    
                    LD = tmpLD(find(tmpLD(:,3)==CommonSlices(sli)),4:end);
                    LV = tmpLV(find(tmpLV(:,3)==CommonSlices(sli)),4:end);
                    
                    RD = mean(tmpRD(find(tmpRD(:,3)==CommonSlices(sli)),4:end),1);
                    RV = mean(tmpRV(find(tmpRV(:,3)==CommonSlices(sli)),4:end),1);
                    
                end
                
                if isequal(CorrType, 'pearson')
                    
                    x = [LD',LV'];
                    rho = corr(x);
                    
                    if isequal(roiTS,'mean')
                        eval(['LDV_' sessions{ses}  '(s,sli) = rho(1,2);']);
                    elseif isequal(roiTS, 'voxel')
                        tmp_rho = rho( 1:size(LD,1), end-(size(LV,1)-1):end);
                        eval(['LDV_' sessions{ses}  '(s,sli) = prctile(tmp_rho(:),95);']);
                    end
                    
                elseif isequal(CorrType,'partial')
                    x = [LD',LV'];
                    z = [RD', RV'];
                    rho = partialcorr(x,z);
                    
                    if isequal(roiTS, 'mean')
                        eval(['LDV_' sessions{ses}  '(s,sli) = rho(1,2);']);
                    elseif isequal(roiTS,'voxel')
                        tmp_rho = rho( 1:size(LD,1), end-(size(LV,1)-1):end);
                        eval(['LDV_' sessions{ses}  '(s,sli) =prctile(tmp_rho(:),95);']);
                    end
                end
                clear x z rho tmp_rho
                clear LD RD LV RV
                
            end
            
            
            for l = 1:numel(levS)
                
                eval(['LDV_' sessions{ses} '_lv(s,l) = '  'mean(LDV_' sessions{ses} '(s,subjectLevels==levS(l)), ''omitnan'');'])
                
            end
            
            
            for sli = 1:numel(CommonSlices)
                
                if isequal(roiTS,'mean')
                    
                    RD = mean(tmpRD(find(tmpRD(:,3)==CommonSlices(sli)),4:end),1);
                    RV = mean(tmpRV(find(tmpRV(:,3)==CommonSlices(sli)),4:end),1);
                    
                    LD = mean(tmpLD(find(tmpLD(:,3)==CommonSlices(sli)),4:end),1);
                    LV = mean(tmpLV(find(tmpLV(:,3)==CommonSlices(sli)),4:end),1);
                    
                elseif isequal(roiTS,'voxel')
                    
                    RD = tmpRD(find(tmpRD(:,3)==CommonSlices(sli)),4:end);
                    RV = tmpRV(find(tmpRV(:,3)==CommonSlices(sli)),4:end);
                    
                    LD = mean(tmpLD(find(tmpLD(:,3)==CommonSlices(sli)),4:end),1);
                    LV = mean(tmpLV(find(tmpLV(:,3)==CommonSlices(sli)),4:end),1);
                    
                end
                
                if isequal(CorrType, 'pearson')
                    
                    x = [RD',RV'];
                    rho = corr(x);
                    
                    if isequal(roiTS,'mean')
                        eval(['RDV_' sessions{ses}  '(s,sli) = rho(1,2);']);
                    elseif isequal(roiTS, 'voxel')
                        tmp_rho = rho( 1:size(RD,1), end-(size(RV,1)-1):end);
                        eval(['RDV_' sessions{ses}  '(s,sli) = prctile(tmp_rho(:),95);']);
                    end
                    
                elseif isequal(CorrType,'partial')
                    x = [RD',RV'];
                    z = [LD', LV'];
                    rho = partialcorr(x,z);
                    
                    if isequal(roiTS, 'mean')
                        eval(['RDV_' sessions{ses}  '(s,sli) = rho(1,2);']);
                    elseif isequal(roiTS,'voxel')
                        tmp_rho = rho( 1:size(RD,1), end-(size(RV,1)-1):end);
                        eval(['RDV_' sessions{ses}  '(s,sli) =prctile(tmp_rho(:),95);']);
                    end
                end
                clear x z rho tmp_rho
                clear LD RD LV RV
                
            end
            
            for l = 1:numel(levS)
                
                eval(['RDV_' sessions{ses} '_lv(s,l) = '  'mean(RDV_' sessions{ses} '(s,subjectLevels==levS(l)), ''omitnan'');'])
                
            end
            
            
            for sli = 1:numel(CommonSlices)
                
                if isequal(roiTS,'mean')
                    
                    LD = mean(tmpLD(find(tmpLD(:,3)==CommonSlices(sli)),4:end),1);
                    RV = mean(tmpRV(find(tmpRV(:,3)==CommonSlices(sli)),4:end),1);
                    
                    RD = mean(tmpRD(find(tmpRD(:,3)==CommonSlices(sli)),4:end),1);
                    LV = mean(tmpLV(find(tmpLV(:,3)==CommonSlices(sli)),4:end),1);
                    
                    
                elseif isequal(roiTS,'voxel')
                    
                    LD = tmpLD(find(tmpLD(:,3)==CommonSlices(sli)),4:end);
                    RV = tmpRV(find(tmpRV(:,3)==CommonSlices(sli)),4:end);
                    
                    RD = mean(tmpRD(find(tmpRD(:,3)==CommonSlices(sli)),4:end),1);
                    LV = mean(tmpLV(find(tmpLV(:,3)==CommonSlices(sli)),4:end),1);
                    
                end
                
                if isequal(CorrType, 'pearson')
                    
                    x = [LD',RV'];
                    rho = corr(x);
                    
                    if isequal(roiTS,'mean')
                        eval(['LDRV_' sessions{ses}  '(s,sli) = rho(1,2);']);
                    elseif isequal(roiTS, 'voxel')
                        tmp_rho = rho( 1:size(LD,1), end-(size(RV,1)-1):end);
                        eval(['LDRV_' sessions{ses}  '(s,sli) = prctile(tmp_rho(:),95);']);
                    end
                    
                elseif isequal(CorrType,'partial')
                    x = [LD',RV'];
                    z = [RD', LV'];
                    rho = partialcorr(x,z);
                    
                    if isequal(roiTS, 'mean')
                        eval(['LDRV_' sessions{ses}  '(s,sli) = rho(1,2);']);
                    elseif isequal(roiTS,'voxel')
                        tmp_rho = rho( 1:size(LD,1), end-(size(RV,1)-1):end);
                        eval(['LDRV_' sessions{ses}  '(s,sli) =prctile(tmp_rho(:),95);']);
                    end
                end
                clear x z rho tmp_rho
                clear LD RD LV RV
                
            end
            
            for l = 1:numel(levS)
                
                eval(['LDRV_' sessions{ses} '_lv(s,l) = '  'mean(LDRV_' sessions{ses} '(s,subjectLevels==levS(l)), ''omitnan'');'])
                
            end
            
            for sli = 1:numel(CommonSlices)
                
                if isequal(roiTS,'mean')
                    
                    RD = mean(tmpRD(find(tmpRD(:,3)==CommonSlices(sli)),4:end),1);
                    LV = mean(tmpLV(find(tmpLV(:,3)==CommonSlices(sli)),4:end),1);
                    
                    LD = mean(tmpLD(find(tmpLD(:,3)==CommonSlices(sli)),4:end),1);
                    RV = mean(tmpRV(find(tmpRV(:,3)==CommonSlices(sli)),4:end),1);
                    
                elseif isequal(roiTS,'voxel')
                    
                    RD = tmpRD(find(tmpRD(:,3)==CommonSlices(sli)),4:end);
                    LV = tmpLV(find(tmpLV(:,3)==CommonSlices(sli)),4:end);
                    
                    LD = mean(tmpLD(find(tmpLD(:,3)==CommonSlices(sli)),4:end),1);
                    RV = mean(tmpRV(find(tmpRV(:,3)==CommonSlices(sli)),4:end),1);
                    
                end
                
                if isequal(CorrType, 'pearson')
                    
                    x = [RD',LV'];
                    rho = corr(x);
                    
                    if isequal(roiTS,'mean')
                        eval(['RDLV_' sessions{ses}  '(s,sli) = rho(1,2);']);
                    elseif isequal(roiTS, 'voxel')
                        tmp_rho = rho( 1:size(RD,1), end-(size(LV,1)-1):end);
                        eval(['RDLV_' sessions{ses}  '(s,sli) = prctile(tmp_rho(:),95);']);
                    end
                    
                elseif isequal(CorrType,'partial')
                    x = [RD',LV'];
                    z = [LD', RV'];
                    rho = partialcorr(x,z);
                    
                    if isequal(roiTS, 'mean')
                        eval(['RDLV_' sessions{ses}  '(s,sli) = rho(1,2);']);
                    elseif isequal(roiTS,'voxel')
                        tmp_rho = rho( 1:size(RD,1), end-(size(LV,1)-1):end);
                        eval(['RDLV_' sessions{ses}  '(s,sli) =prctile(tmp_rho(:),95);']);
                    end
                end
                clear x z rho tmp_rho
                clear LD RD LV RV
                
            end
            
            for l = 1:numel(levS)
                
                eval(['RDLV_' sessions{ses} '_lv(s,l) = '  'mean(RDLV_' sessions{ses} '(s,subjectLevels==levS(l)), ''omitnan'');'])
                
            end
            
        end
        
    end
end
%%

for l = 1:numel(levS)
    
    
    eval(['within_' sessions{1}  '_lv(:,l) = (LDV_' sessions{1} '_lv(:,l)+RDV_' sessions{1} '_lv(:,l))./2;'])
    eval(['between_' sessions{1} '_lv(:,l) = (LDRV_' sessions{1} '_lv(:,l)+RDLV_' sessions{1} '_lv(:,l))./2;'])
    
    eval(['within_' sessions{2} '_lv(:,l) = (LDV_' sessions{2} '_lv(:,l)+RDV_' sessions{2} '_lv(:,l))./2;'])
    eval(['between_' sessions{2} '_lv(:,l) = (LDRV_' sessions{2} '_lv(:,l)+RDLV_' sessions{2} '_lv(:,l))./2;'])
    
    
    if CalculateICC
        
        clear tmp1 tmp2 X
        
        eval(['tmp1 = DD_' sessions{1} '_lv(:,l);']);
        eval(['tmp2 = DD_' sessions{2} '_lv(:,l);']);
        
        tmp1 = tmp1(~isnan(tmp1));
        tmp2 = tmp2(~isnan(tmp2));
        
        X = [tmp1, tmp2];
        ICCvalue(1,l,1) =  ICC(X, '2', '1');
        ICCvalue_CI(1,l,:) = bootci(100000,@myICC,X);
        
        clear tmp1 tmp2 X ci
        
        eval(['tmp1 = VV_' sessions{1} '_lv(:,l);']);
        eval(['tmp2 = VV_' sessions{2} '_lv(:,l);']);
        
        tmp1 = tmp1(~isnan(tmp1));
        tmp2 = tmp2(~isnan(tmp2));
        
        X = [tmp1, tmp2];
        ICCvalue(2,l,1) =  ICC(X, '2', '1');
        ICCvalue_CI(2,l,:) = bootci(100000,@myICC,X);
        
        %
        clear tmp1 tmp2 X ci
        
        eval(['tmp1 = within_' sessions{1} '_lv(:,l);']);
        eval(['tmp2 = within_' sessions{2} '_lv(:,l);']);
        
        tmp1 = tmp1(~isnan(tmp1));
        tmp2 = tmp2(~isnan(tmp2));
        
        
        X = [tmp1, tmp2];
        ICCvalue(3,l,1) =  ICC(X, '2', '1');
        ICCvalue_CI(3,l,:) = bootci(100000,@myICC,X);
        
        clear tmp1 tmp2 X ci
        
        eval(['tmp1 = between_' sessions{1} '_lv(:,l);']);
        eval(['tmp2 = between_' sessions{2} '_lv(:,l);']);
        
        tmp1 = tmp1(~isnan(tmp1));
        tmp2 = tmp2(~isnan(tmp2));
        
        
        X = [tmp1, tmp2];
        ICCvalue(4,l,1) =  ICC(X, '2', '1');
        ICCvalue_CI(4,l,:) = bootci(100000,@myICC,X);
             
        
    end
    
    
    DD_session1 = eval(['DD_' sessions{1} '_lv;']);
    DD_session2  = eval(['DD_' sessions{2} '_lv;']);
    VV_session1 = eval(['VV_' sessions{1} '_lv;']);
    VV_session2  = eval(['VV_' sessions{2} '_lv;']);
    LDV_session1 = eval(['LDV_' sessions{1} '_lv;']);
    LDV_session2 = eval(['LDV_' sessions{2} '_lv;']);
    RDV_session1 = eval(['RDV_' sessions{1} '_lv;']);
    RDV_session2 = eval(['RDV_' sessions{2} '_lv;']);
    within_session1 = eval(['within_' sessions{1} '_lv;']);
    within_session2 = eval(['within_' sessions{2} '_lv;']);
    between_session1  = eval(['between_' sessions{1} '_lv;']);
    between_session2  = eval(['between_' sessions{2} '_lv;']);
    
end
%%
