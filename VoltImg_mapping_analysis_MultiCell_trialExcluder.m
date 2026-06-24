%% Trial Excluder

% Create structs for all cells 
stdImagingAllTrialsCellNames = cell(nCells, 1);
stdFiltImagingAllTrialsCellNames = cell(nCells, 1);
for nn = 1:nCells
    stdImagingAllTrialsCellNames{nn} = ['stdImagingAllTrialsCellNames_', 'cell', num2str(nn)];
    stdFiltImagingAllTrialsCellNames{nn} = ['stdFiltImagingAllTrialsCellNames_', 'cell', num2str(nn)];
    analysisStruct.(stdImagingAllTrialsCellNames{nn}) = cell(nConds, 1);
    analysisStruct.(stdFiltImagingAllTrialsCellNames{nn}) = cell(nConds, 1);

    for cc = 1:nConds 
        for hh = 1:nHolos(cc)
            analysisStruct.(stdImagingAllTrialsCellNames{nn}){cc}{hh, 1} = std(analysisStruct.(holoSortedImagingCellNames{nn}){cc}{hh}(:));
            analysisStruct.(stdFiltImagingAllTrialsCellNames{nn}){cc}{hh, 1} = std(analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc}{hh}(:));
        end
    end
end

exclHoloSortedImagingAllTrialsCellNames = cell(nCells, 1);
exclFiltHoloSortedImagingAllTrialsCellNames = cell(nCells, 1);
for nn = 1:nCells
    exclHoloSortedImagingAllTrialsCellNames{nn} = ['exclHoloSortedImagingAllTrials_', 'cell', num2str(nn)];
    exclFiltHoloSortedImagingAllTrialsCellNames{nn} = ['exclFiltHoloSortedImagingAllTrials_', 'cell', num2str(nn)];
    analysisStruct.(exclHoloSortedImagingAllTrialsCellNames{nn}) = cell(nConds, 1);
    analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames{nn}) = cell(nConds, 1);

    for cc = 1:nConds
        analysisStruct.(exclHoloSortedImagingAllTrialsCellNames{nn}){cc} = cell(nHolos(cc), 1);
        analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames{nn}){cc} = cell(nHolos(cc), 1);
    end
end

% 
for nn = 1:nCells
    for cc = 1:nConds 
        for hh = 1:nHolos(cc)
            for tt = 1:size(analysisStruct.(holoSortedImagingCellNames{nn}){cc}{hh}, 2)
                % if any(analysisStruct.(holoSortedImagingCellNames{nn}){cc}{hh}(:, tt) > 2*analysisStruct.(stdImagingAllTrialsCellNames{nn}){cc}{hh, 1})
                if any(analysisStruct.(holoSortedImagingCellNames{nn}){cc}{hh}(:, tt) < -2.5*analysisStruct.(stdImagingAllTrialsCellNames{nn}){cc}{hh, 1}) %| any(analysisStruct.(holoSortedImagingCellNames{nn}){cc}{hh}(:, tt) > 3*analysisStruct.(stdImagingAllTrialsCellNames{nn}){cc}{hh, 1})
                    analysisStruct.(exclHoloSortedImagingAllTrialsCellNames{nn}){cc}{hh, 1}(:, tt) = nan(size(analysisStruct.(holoSortedImagingCellNames{nn}){cc}{hh}(:, tt), 1), 1);
                else 
                    analysisStruct.(exclHoloSortedImagingAllTrialsCellNames{nn}){cc}{hh, 1}(:, tt) = analysisStruct.(holoSortedImagingCellNames{nn}){cc}{hh}(:, tt);
                end
    
                % if any(analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc}{hh}(:, tt) > 2*analysisStruct.(stdFiltImagingAllTrialsCellNames{nn}){cc}{hh, 1})
                if any(analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc}{hh}(:, tt) < -2.5*analysisStruct.(stdFiltImagingAllTrialsCellNames{nn}){cc}{hh, 1}) %| any(analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc}{hh}(:, tt) > 3*analysisStruct.(stdFiltImagingAllTrialsCellNames{nn}){cc}{hh, 1}) 
                    analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames{nn}){cc}{hh, 1}(:, tt) = nan(size(analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc}{hh}(:, tt), 1), 1);
                else 
                    analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames{nn}){cc}{hh, 1}(:, tt) = analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc}{hh}(:, tt);
                end
            end
        end
    end
end

for nn = 1:nCells
    for cc = 1:nConds 
        for hh = 1:nHolos(cc)
            for tt = 1:size(analysisStruct.(exclHoloSortedImagingAllTrialsCellNames{nn}){cc}{hh}, 2)
                if any(isnan(analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames{nn}){cc}{hh}(:, tt)))
                    continue
                end
                [~, maxImagingIndex] = max(analysisStruct.(exclHoloSortedImagingAllTrialsCellNames{nn}){cc}{hh}(:, tt));
                [~, maxFiltImagingIndex] = max(analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames{nn}){cc}{hh}(:, tt));
                if maxImagingIndex > 45
                    analysisStruct.(exclHoloSortedImagingAllTrialsCellNames{nn}){cc}{hh}(:, tt) = nan(size(analysisStruct.(exclHoloSortedImagingAllTrialsCellNames{nn}){cc}{hh}(:, tt), 1), 1);
                end
                if maxFiltImagingIndex > 45
                    analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames{nn}){cc}{hh}(:, tt) = nan(size(analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames{nn}){cc}{hh}(:, tt), 1), 1);
                end
            end
        end
    end
end

%% Calculate mean response (and CI) for each hologram across trials and per condition
exclHoloSortedImagingMeanCellNames = cell(nCells, 1);
exclFiltHoloSortedImagingMeanCellNames = cell(nCells, 1);
for nn = 1:nCells
    exclHoloSortedImagingMeanCellNames{nn} = ['exclHoloSortedImagingMean_', 'cell', num2str(nn)];
    exclFiltHoloSortedImagingMeanCellNames{nn} = ['exclFiltHoloSortedImagingMean_', 'cell', num2str(nn)];
    analysisStruct.(exclHoloSortedImagingMeanCellNames{nn}) = cell(nConds, 1);
    analysisStruct.(exclFiltHoloSortedImagingMeanCellNames{nn}) = cell(nConds, 1);
    
    for cc = 1:nConds
        analysisStruct.(exclHoloSortedImagingMeanCellNames{nn}){cc} = cell(nHolos(cc), 1);
        analysisStruct.(exclFiltHoloSortedImagingMeanCellNames{nn}){cc} = cell(nHolos(cc), 1);
    end
end

for nn = 1:nCells
    for cc = 1:nConds
        for hh = 1:nHolos(cc)
            analysisStruct.(exclHoloSortedImagingMeanCellNames{nn}){cc}{hh} = nanmean(analysisStruct.(exclHoloSortedImagingAllTrialsCellNames{nn}){cc}{hh}, 2);
            analysisStruct.(exclFiltHoloSortedImagingMeanCellNames{nn}){cc}{hh} = nanmean(analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames{nn}){cc}{hh}, 2);        
        end
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

exclCIDffAllCondsCellNames = cell(nCells, 1);
exclFiltCIDffAllCondsCellNames = cell(nCells, 1);
for nn = 1:nCells
    exclCIDffAllCondsCellNames{nn} = ['exclCIDffAllConds_', 'cell', num2str(nn)];
    exclFiltCIDffAllCondsCellNames{nn} = ['exclFiltCIDffAllConds_', 'cell', num2str(nn)];
    analysisStruct.(exclCIDffAllCondsCellNames{nn}) = cell(nConds, 1);
    analysisStruct.(exclFiltCIDffAllCondsCellNames{nn}) = cell(nConds, 1);
end

for nn = 1:nCells
    for cc = 1:nConds 
        for hh = 1:nHolos(cc)
            confidence_level = 0.95;
            means = nanmean(analysisStruct.(exclHoloSortedImagingAllTrialsCellNames{nn}){cc}{hh, 1}, 2);
            filtMeans = nanmean(analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames{nn}){cc}{hh, 1}, 2);       
            std_errors = std(analysisStruct.(exclHoloSortedImagingAllTrialsCellNames{nn}){cc}{hh, 1}, 0, 2, "omitnan") / sqrt(size(analysisStruct.(exclHoloSortedImagingAllTrialsCellNames{nn}){cc}{hh, 1}, 2));
            filtStd_errors = std(analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames{nn}){cc}{hh, 1}, 0, 2, "omitnan") / sqrt(size(analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames{nn}){cc}{hh, 1}, 2));
      
            t_score = tinv((1 + confidence_level) / 2, size(analysisStruct.(exclHoloSortedImagingAllTrialsCellNames{nn}){cc}{hh, 1}, 2) - 1);
            filtT_score = tinv((1 + confidence_level) / 2, size(analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames{nn}){cc}{hh, 1}, 2) - 1);
            margin_of_error = t_score * std_errors;
            filtMargin_of_error = filtT_score * filtStd_errors;
            lower_bounds = means - margin_of_error;
            filtLower_bounds = filtMeans - filtMargin_of_error;
            upper_bounds = means + margin_of_error;
            filtUpper_bounds = filtMeans + filtMargin_of_error;
            if UpOrDown == '2'
                analysisStruct.(exclCIDffAllCondsCellNames{nn}){cc}{hh, 1} = [lower_bounds, upper_bounds];
                analysisStruct.(exclFiltCIDffAllCondsCellNames{nn}){cc}{hh, 1} = [filtLower_bounds, filtUpper_bounds];
            elseif UpOrDown =='1'
                analysisStruct.(exclCIDffAllCondsCellNames{nn}){cc}{hh, 1} = [-lower_bounds, -upper_bounds];
                analysisStruct.(exclFiltCIDffAllCondsCellNames{nn}){cc}{hh, 1} = [-filtLower_bounds, -filtUpper_bounds];
            end
        end
    end
end
