function  REL_helperPlotConnectivity(subjects, saveDir,yLims, fileNames, analysisSpace, ...
    order,sessions,tSNRCorrPlot,indLines,dataDir,plotType, printPath, printMode)

close all
%%
if isequal(analysisSpace, 'ts')
    aspace = 'templspace'  %'templspace'  %or subjspace
elseif isequal(analysisSpace, 'ns')
    aspace = 'nativespace'
end


for o = 1:numel(order)
    
    if isequal(plotType,'avgOverSes')
        
        for f = 1:numel(fileNames)
            
            load([saveDir filesep order{o} '_' aspace  '_' fileNames{f} 'corrs.mat']);
            
            DD = (DD_session1+DD_session2)./2;
            VV = (VV_session1+VV_session2)./2;
            wDV = (within_session1+within_session2)./2;
            bDV = (between_session1+between_session2)./2;
            
            corrDD{f} = DD;
            corrVV{f} = VV;
            corrwDV{f} = wDV;
            corrbDV{f} = bDV;
            
            clear DD_session1 DD_session2 VV_session1 VV_session2 within_session1 within_session2 ...
                between_session1 between_session2 DD VV wDV bD
            
            
        end
        
        connections = {'DD', 'VV', 'wDV', 'bDV'};
        %
        %         colorS(1,:) = [8,48,107]./255;
        %         colorS(2,:) = [33,113,181]./255 ;
        %         colorS(3,:) = [66,146,198]./255;
        %         colorS(4,:) = [107,174,214]./255;
        
        colorS(1,:) =  [102,194,165]./255;
        colorS(2,:) =  [252,141,98]./255;
        colorS(3,:) =  [141,160,203]./255;
        colorS(4,:) = [231,138,195]./255;
        
        
        for c = 1:numel(connections)
            
            fiGure = figure;
            fiGure.Position = [0 0 800 800];
            hold on
            box off
            
            ylabel('r');
            XLabels= fileNames;
            yLimVal = yLims;
            xlim([0 numel(fileNames)+1]);
            xticks([1:numel(fileNames)]);
            ylim(yLimVal)
            %yticks(-0.5:0.2:1)
            ax = gca;
            ax.FontSize = 20;
            %subtitle(order{o})
            yline(0)
            
            connName = connections{c};
            
            title([connName ]) %' N = ' num2str(numel(subjects))
            
            positions = 1:numel(fileNames);
            
            for f = 1:numel(fileNames)
                
                h =  boxplot(eval(['corr' connName '{f}']), 'positions', ...
                    ones(numel(eval(['corr' connName '{f}'])),1).*positions(f) ,...
                    'Colors', colorS(c,:) ,...
                    'Symbol','.r', ...
                    'Widths', 0.3);
                set(h,'linew',3)
                lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
                set(lines, 'Color', 'k');

                line([positions(f)-0.15 positions(f)+0.15], ...
                    [mean(eval(['corr' connName '{f}']) , 'omitnan')  mean(eval(['corr' connName '{f}']) , 'omitnan')],  ...
                    'LineWidth', 3, 'Color', colorS(c,:))
                
                
%                 plot(positions(f) , ...
%                     mean(eval(['corr' connName '{f}']), 'omitnan'), ...
%                     'r*','MarkerSize', 15, 'LineWidth',2)
%                 
                
            end
            
            
            if indLines
                
                for p = 1:size(eval(['corr' connName]),2)-1
                    
                    for j = 1:numel(eval(['corr' connName '{p}']))
                        
                        line([positions(p)+0.2 positions(p+1)-0.2], [eval(['corr' connName '{p}(j)']) eval(['corr' connName '{p+1}(j)'])], ...
                            'color',[colorS(c,:) 0.5] , 'linewidth', 0.8);
                    end
                end
                
                
                
                for m = 1:numel(fileNames)
                    
                    meanValues(m) = mean(eval(['corr' connName '{m}']), 'omitnan');
                    
                end
                
                for p = 1:size(eval(['corr' connName]),2)-1
                    
                    
                    line([positions(p)+0.2 positions(p+1)-0.2], [meanValues(p) meanValues(p+1)], ...
                        'color',[0 0 0] , 'linewidth',4);
                    
                end
                
            else
                
                tmp = eval(['corr' connName ]);
                tmp = cell2mat(tmp);
                
                for i = 1:size(tmp,2)
                    
                    scatter(ones(size(tmp,1),1).*(positions(i)+(rand(size(tmp,1),1)-0.5)/10),tmp(:,i),[] ,[colorS(c,:)] ,'filled', 'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.4)
                    
                end
                
                clear tmp
                
            end
            
            
            
            
            xticks(positions)
            xticklabels(fileNames)
            set(gca,'TickLabelInterpreter','none')
            
                  
            grid on
            box  off
            ylim(yLimVal)
            
            
            if printMode
                
                print([ printPath filesep ...
                    'Figure_' order{o} '_' aspace  '_' connName '_' plotType  '_Corr.tif'],'-dtiff')
                
                print([ printPath filesep ...
                    'Figure_' order{o} '_' aspace  '_' connName '_' plotType  '_Corr.svg'],'-dsvg')
                
                
            end
            
        end
    end
    
    if tSNRCorrPlot
        
        for sub = 1:numel(subjects)
            
            cd(fullfile(dataDir,subjects{sub},'func', 'tsnrFiles_censored'))
            
            for f = 1:numel(fileNames)
                
                for ses = 1:numel(sessions)
                    
                    [tmpLD,tmpRD, tmpLV, tmpRV] = REL_helper_loadRois('avgTS','manual', fileNames{f}, 'tsnr');
                    TSNR(sub,f,ses,1) = mean(tmpLD + tmpRD)./4;
                    TSNR(sub,f,ses,2) = mean(tmpLV + tmpRV)./4;
                    clear tmpLD tmpRD tmpLV tmpRV
                    
                end
                
            end
        end
        
        TSNR = squeeze(mean(TSNR, 3,'omitnan'));
        xLimVal_tSNR = [5 12];
        
        for j = 1:2
            
            fiGure = figure;
            fiGure.Position = [0 0 3000 3000];
            hold on
            
            TSNR_values = squeeze(TSNR(:,:,j));
            
            
            if j == 1
                
                sgtitle(['DD, ' order{o} ' N = ' num2str(size(TSNR_values,1))])
                connName = 'DD';
                
            elseif j == 2
                
                sgtitle(['VV, ' order{o}  ' N = ' num2str(size(TSNR_values,1))])
                connName = 'VV';
                
            end
            
            for i = 1:size(TSNR_values,2)
                
                subplot(2,ceil((size(TSNR_values,2)/2)), i); hold on;
                subtitle(fileNames{i}, 'Interpreter', 'none');
                
                if j == 1
                    
                    scatter(TSNR_values(:,i), corrDD{i}, 'filled')
                    R= corrcoef(TSNR_values(:,i), corrDD{i});
                    
                elseif j == 2
                    
                    scatter(TSNR_values(:,i), corrVV{i}, 'filled')
                    R= corrcoef(TSNR_values(:,i), corrVV{i});
                end
                
                
                h2 = lsline;
                h2.LineWidth = 2;
                h2.Color = 'k';
                
                ylim(yLimVal)
                xlim(xLimVal_tSNR)
                axis square
                yline(0)
                xlabel('tSNR')
                ylabel('r')
                ax = gca;
                ax.FontSize = 14;
                legend(['R^2 = ' num2str(round(R(1,2)^2,3)) ]);
                
                
                
                
            end
      
            
            if printMode
                
                print([ printPath filesep ...
                    'Figure_' sessions{:} '_' connName '_tSNRCorr_' num2str(numel(subjects)) '_' order{o} '.tif'],'-dtiff')
                
                print([ printPath filesep ...
                    'Figure_' sessions{:} '_' connName '_tSNRCorr_' num2str(numel(subjects)) '_' order{o} '.svg'],'-dsvg')
            end
            
            
        end
        
    end
    
end




end

