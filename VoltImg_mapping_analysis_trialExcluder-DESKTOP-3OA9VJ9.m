%%
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
            if any(holoSortedImagingAllTrials{cc}{hh}(:, tt) < -2*stdImagingAllTrials{cc}{hh, 1}) %| any(holoSortedImagingAllTrials{cc}{hh}(:, tt) > 2*stdImagingAllTrials{cc}{hh, 1})
                exclHoloSortedImagingAllTrials{cc}{hh, 1}(:, tt) = nan(size(holoSortedImagingAllTrials{cc}{hh}(:, tt), 1), 1);
            else 
                exclHoloSortedImagingAllTrials{cc}{hh, 1}(:, tt) = holoSortedImagingAllTrials{cc}{hh}(:, tt);
            end

            if any(filtHoloSortedImagingAllTrials{cc}{hh}(:, tt) < -2*stdFiltImagingAllTrials{cc}{hh, 1}) %| any(filtHoloSortedImagingAllTrials{cc}{hh}(:, tt) > 2*stdFiltImagingAllTrials{cc}{hh, 1})
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
            if maxImagingIndex > 45
                exclHoloSortedImagingAllTrials{cc}{hh}(:, tt) = nan(size(exclHoloSortedImagingAllTrials{cc}{hh}(:, tt), 1), 1);
            end
            if maxFiltImagingIndex > 45
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

voltMapping.exclHoloSortedImagingAllTrials     = exclHoloSortedImagingAllTrials;
voltMapping.exclFiltHoloSortedImagingAllTrials = exclFiltHoloSortedImagingAllTrials;
voltMapping.exclHoloSortedImagingMean          = exclHoloSortedImagingMean;
voltMapping.exclfiltHoloSortedImagingMean      = exclfiltHoloSortedImagingMean;
voltMapping.exclCIDffAllConds                  = exclCIDffAllConds;
voltMapping.exclFiltCIDffAllConds              = exclFiltCIDffAllConds;

