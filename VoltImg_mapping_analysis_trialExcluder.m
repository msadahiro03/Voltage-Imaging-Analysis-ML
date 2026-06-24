% stdOfImaging = std(holoSortedImagingAllTrials{1, 1}{5, 1}(:));
stdImagingAllTrials = cell(nConds, 1);
stdFiltImagingAllTrials = cell(nConds, 1);
for cc = 1:nConds 
    for hh = 1:nHolos(cc)
        stdImagingAllTrials{cc}{hh, 1} = std(holoSortedImagingAllTrials{cc}{hh}(:));
        stdFiltImagingAllTrials{cc}{hh, 1} = std(filtHoloSortedImagingAllTrials{cc}{hh}(:));
    end
end

exclHoloSortedImagingAllTrials = cell(nConds, 1);
exclFiltHoloSortedImagingAllTrials = cell(nConds, 1);
for cc = 1:nConds 
    for hh = 1:nHolos(cc)
        for tt = 1:size(holoSortedImagingAllTrials{cc}{hh}, 2)
            if any(holoSortedImagingAllTrials{cc}{hh}(:, tt) < -2*stdImagingAllTrials{cc}{hh, 1})
                exclHoloSortedImagingAllTrials{cc}{hh, 1}(:, tt) = nan(size(holoSortedImagingAllTrials{cc}{hh}(:, tt), 1), 1);
            else 
                exclHoloSortedImagingAllTrials{cc}{hh, 1}(:, tt) = holoSortedImagingAllTrials{cc}{hh}(:, tt);
            end

            if any(filtHoloSortedImagingAllTrials{cc}{hh}(:, tt) < -2*stdFiltImagingAllTrials{cc}{hh, 1})
                exclFiltHoloSortedImagingAllTrials{cc}{hh, 1}(:, tt) = nan(size(filtHoloSortedImagingAllTrials{cc}{hh}(:, tt), 1), 1);
            else 
                exclFiltHoloSortedImagingAllTrials{cc}{hh, 1}(:, tt) = filtHoloSortedImagingAllTrials{cc}{hh}(:, tt);
            end
        end
    end
end

for cc = 1:nConds 
    for hh = 1:nHolos(cc)
        for tt = 1:size(exclHoloSortedImagingAllTrials{cc}{hh}, 2)
            if any(isnan(exclFiltHoloSortedImagingAllTrials{cc}{hh}(:, tt)))
                continue
            end
            [~, maxImagingIndex] = max(exclHoloSortedImagingAllTrials{cc}{hh}(:, tt));
            [~, maxFiltImagingIndex] = max(exclFiltHoloSortedImagingAllTrials{cc}{hh}(:, tt));
            if maxImagingIndex > 40
                exclHoloSortedImagingAllTrials{cc}{hh}(:, tt) = nan(size(exclHoloSortedImagingAllTrials{cc}{hh}(:, tt), 1), 1);
            end
            if maxFiltImagingIndex > 40
                exclFiltHoloSortedImagingAllTrials{cc}{hh}(:, tt) = nan(size(exclFiltHoloSortedImagingAllTrials{cc}{hh}(:, tt), 1), 1);
            end
        end
    end
end

%% Calculate mean response (and CI) for each hologram across trials and per condition

exclHoloSortedImagingMean = cell(nConds, 1); % The mean response for each hologram across conditions
exclfiltHoloSortedImagingMean = cell(nConds, 1);
for cc = 1:nConds
    exclHoloSortedImagingMean{cc} = cell(nHolos(cc), 1);
    exclfiltHoloSortedImagingMean{cc} = cell(nHolos(cc), 1);
end
for cc = 1:nConds
    for hh = 1:nHolos(cc)
        exclHoloSortedImagingMean{cc}{hh} = nanmean(exclHoloSortedImagingAllTrials{cc}{hh}, 2);
        exclfiltHoloSortedImagingMean{cc}{hh} = nanmean(exclFiltHoloSortedImagingAllTrials{cc}{hh}, 2);        
        % Baseline the mean holo traces
%         holoSortedImagingMean{cc}{hh} = holoSortedImagingMean{cc}{hh} - mean(holoSortedImagingMean{cc}{hh}(1:preStimWindow/1000*imagingFreq));
%         filtHoloSortedImagingMean{cc}{hh} = filtHoloSortedImagingMean{cc}{hh} - mean(filtHoloSortedImagingMean{cc}{hh}(1:preStimWindow/1000*imagingFreq));
      
    %     figure(30);
    %     clf
    %     for nn = 1:length(nPulseCoordsImaging)
    %         xline(nPulseCoordsImaging(nn)/imagingFreq, '--', 'LineWidth', 1.5, 'color', [.8 .8 .8]);
    %     end
    %     hold on
    %     plot(linspace(0, length(holoSortedImagingMean{cc}{hh})/imagingFreq, length(holoSortedImagingMean{cc}{hh})), holoSortedImagingMean{cc}{hh}, 'LineWidth', 1.5);
    %     plot(linspace(0, length(filtHoloSortedImagingMean{cc}{hh})/imagingFreq, length(filtHoloSortedImagingMean{cc}{hh})), filtHoloSortedImagingMean{cc}{hh}, 'LineWidth', 1.5);
    % %     axis off
    %     hold off
    %     xlabel('Time(s)')
    %     ylabel('df/f')
        % pause
    end
end

% % Combine all hologram traces, calculate one grand mean. Not a useful step - wrote this for the hell of it.
% holoComboImagingAllTrials = [];
% for cc = 1
%     for hh = 1:nHolos(1)
%         holoComboImagingAllTrials = [holoComboImagingAllTrials, holoSortedImagingAllTrials_ALT{cc}{hh}(:, :)];
%     end
% end
% holoComboImagingGrandMean = nanmean(holoComboImagingAllTrials, 2);

exclCIDffAllConds = cell(nConds, 1);
exclFiltCIDffAllConds = cell(nConds, 1);
for cc = 1:nConds 
    for hh = 1:nHolos(cc)
        confidence_level = 0.95;
        means = nanmean(exclHoloSortedImagingAllTrials{cc}{hh, 1}, 2);
        filtMeans = nanmean(exclFiltHoloSortedImagingAllTrials{cc}{hh, 1}, 2);       
        std_errors = std(exclHoloSortedImagingAllTrials{cc}{hh, 1}, 0, 2, "omitnan") / sqrt(size(exclHoloSortedImagingAllTrials{cc}{hh, 1}, 2));
        filtStd_errors = std(exclFiltHoloSortedImagingAllTrials{cc}{hh, 1}, 0, 2, "omitnan") / sqrt(size(exclFiltHoloSortedImagingAllTrials{cc}{hh, 1}, 2));
  
        t_score = tinv((1 + confidence_level) / 2, size(exclHoloSortedImagingAllTrials{cc}{hh, 1}, 2) - 1);
        filtT_score = tinv((1 + confidence_level) / 2, size(exclFiltHoloSortedImagingAllTrials{cc}{hh, 1}, 2) - 1);
        margin_of_error = t_score * std_errors;
        filtMargin_of_error = filtT_score * filtStd_errors;
        lower_bounds = means - margin_of_error;
        filtLower_bounds = filtMeans - filtMargin_of_error;
        upper_bounds = means + margin_of_error;
        filtUpper_bounds = filtMeans + filtMargin_of_error;
        if UpOrDown == '2'
            exclCIDffAllConds{cc}{hh, 1} = [lower_bounds, upper_bounds];
            exclFiltCIDffAllConds{cc}{hh, 1} = [filtLower_bounds, filtUpper_bounds];
        elseif UpOrDown =='1'
            exclCIDffAllConds{cc}{hh, 1} = [-lower_bounds, -upper_bounds];
            exclFiltCIDffAllConds{cc}{hh, 1} = [-filtLower_bounds, -filtUpper_bounds];

        end
    end
end

voltMapping.holoSortedImagingAllTrials         = exclHoloSortedImagingAllTrials;
voltMapping.filtHoloSortedImagingAllTrials     = exclFiltHoloSortedImagingAllTrials;
voltMapping.exclHoloSortedImagingMean          = exclHoloSortedImagingMean;
voltMapping.exclfiltHoloSortedImagingMean      = exclfiltHoloSortedImagingMean;
voltMapping.exclCIDffAllConds                  = exclCIDffAllConds;
voltMapping.exclFiltCIDffAllConds              = exclFiltCIDffAllConds;

%% Align imaging with ephys traces
for cc = 1:nConds
    for hh = 1:nHolos(cc)
        figure(cc*1000+hh);
%         set(gcf, 'Position',  [100, 100, 600, 400])
        clf
        hold on
%         fill([linspace(0, size(CIDfAllConds{cc}{hh}, 1)/imagingFreq, size(CIDfAllConds{cc}{hh}, 1)), fliplr(linspace(0, size(CIDfAllConds{cc}{hh}, 1)/imagingFreq, size(CIDfAllConds{cc}{hh}, 1)))],...
%         [CIDfAllConds{cc}{hh}(:, 1)', fliplr(CIDfAllConds{cc}{hh}(:, 2)')], [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
%         % plot CI lowerbound
%         plot(linspace(0, size(CIDfAllConds{cc}{hh}, 1)/imagingFreq, size(CIDfAllConds{cc}{hh}, 1)), CIDfAllConds{cc}{hh}(:, 1), '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]); 
%         % plot CI upperbound
%         plot(linspace(0, size(CIDfAllConds{cc}{hh}, 1)/imagingFreq, size(CIDfAllConds{cc}{hh}, 1)), CIDfAllConds{cc}{hh}(:, 2), '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);   

        fill([linspace(0, size(exclFiltCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclFiltCIDffAllConds{cc}{hh}, 1)), fliplr(linspace(0, size(exclFiltCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclFiltCIDffAllConds{cc}{hh}, 1)))],...
        [exclFiltCIDffAllConds{cc}{hh}(:, 1)'*100, fliplr(exclFiltCIDffAllConds{cc}{hh}(:, 2)'*100)], [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
        % plot CI lowerbound
        plot(linspace(0, size(exclFiltCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclFiltCIDffAllConds{cc}{hh}, 1)), exclFiltCIDffAllConds{cc}{hh}(:, 1)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]); 
        % plot CI upperbound
        plot(linspace(0, size(exclFiltCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclFiltCIDffAllConds{cc}{hh}, 1)), exclFiltCIDffAllConds{cc}{hh}(:, 2)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);   
        
%         % plot ephys and voltage traces
%         [ax, hl, h2] = plotyy(linspace(0, size(holoSortedImagingMean{cc}(:, hh), 1)/imagingFreq, size(holoSortedImagingMean{cc}(:, hh), 1)), holoSortedImagingMean{cc}(:, hh), linspace(0, size(holoSortedDataMean{cc}(:, hh), 1)/Fs, size(holoSortedDataMean{cc}(:, hh), 1)), holoSortedDataMean{cc}(:, hh));
% %         set(ax, 'XAxisLocation', 'origin');
% %         ylim([min(holoSortedImagingMean{cc}(:, hh)), max(holoSortedImagingMean{cc}(:, hh))]);
% %         ylim([-1, 2])
%         ax(2).YLim = ax(1).YLim;
%         set(ax(1), 'ytick', floor(min(holoSortedImagingMean{cc}(:, hh))):1:ceil(max(holoSortedImagingMean{cc}(:, hh))));
%         set(ax, {'ycolor'},{'g';'k'});
%         set(ax, 'FontSize', 14);
%         ylabel(ax(1), 'dF/F');
%         ylabel(ax(2), 'dV');
%         xlabel(ax(2), 'Time (s)');
%         set(hl, 'Color', 'g');
%         set(h2, 'Color', [0.7 0.7 0.7]);
%         set(hl, 'LineWidth', 2);    

        % plot voltage imaging trace only
        % yyaxis right
        % plot(linspace(0, size(exclHoloSortedImagingMean{cc}{hh}, 1)/imagingFreq, size(exclHoloSortedImagingMean{cc}{hh}, 1)), exclHoloSortedImagingMean{cc}{hh}*100, '-', 'linewidth', 2, 'color', 'g');
        plot(linspace(0, size(exclfiltHoloSortedImagingMean{cc}{hh}, 1)/imagingFreq, size(exclfiltHoloSortedImagingMean{cc}{hh}, 1)), exclfiltHoloSortedImagingMean{cc}{hh}*100, '-', 'linewidth', 2, 'color', 'g');
        ylabel('dF/F (%)');
        xlabel('Time (s)');
%         ylim([-1 max(holoSortedImagingMean{cc}{hh})*100])
%         ylim([min(holoSortedImagingMean{cc}{hh})*100 max(holoSortedImagingMean{cc}{hh})*100])
        ylim([-1.5 3])
        ax.YColor = [0 1 0];
        % xticks([0:0.05:size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);

        % axis off
%         
        % plot ephys trace only
%         yyaxis left
%     %     axes('Position',[.70 .12 .2 .2]);
%     %     box on
%         plot(linspace(0, size(holoSortedDataMean{cc}(:, hh), 1)/Fs, size(holoSortedDataMean{cc}(:, hh), 1)), holoSortedDataMean{cc}(:, hh), 'linewidth', 1.5, 'color', [0 0 0]);
        % gca;
        % set(gca,'xtick',[], 'fontsize', 18);
% %         ylim([min(holoSortedDataMean{cc}(:, hh)), max(holoSortedDataMean{cc}(:, hh))]);
%         ylim([-0.5 2]);
%         ylabel('dV');
%         ax.YColor = [0, 0, 0];
        % axis off
%         xlim([0 size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);
%         xticks([0:0.05:size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);
%         
        % show line at dff = 0
%         yline(0, '-', 'LineWidth', 1.5, 'color', [0.9 0.9 0.9]);
        for nn = 1:length(nPulseCoords)
            xline(nPulseCoords(nn)/Fs, '-', 'LineWidth', 10, 'color', [1 0 0]);
        end
        
        % plot([0.01; 0.035], [0.5; 0.5], '-k', 'LineWidth', 2);
        % plot([0.01; 0.01], [0.5; 1], '-k', 'LineWidth', 2);
        hold off
     pause
    end
end