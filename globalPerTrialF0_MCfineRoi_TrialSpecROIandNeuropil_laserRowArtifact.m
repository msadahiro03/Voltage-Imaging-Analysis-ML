%% ========================================================================
%% Parallel pipeline: trial-common early baseline F0 + dF/F (_commonF0)
%% Requires: original "F, F0, dF, dF/F0" + ROI section already run so that
%%   fineRoiXAllCells, fineRoiYAllCells, bkgrndRoiXAllCells_trial,
%%   bkgrndRoiYAllCells_trial exist; same voltMapping / trials / mc paths.
%% ========================================================================
if ~exist('fineRoiXAllCells', 'var') || ~exist('fineRoiYAllCells', 'var')
    error('VoltImg:commonF0NeedsFineRois', ...
        'Run the original F/ROI/dF/F section first so fineRoiXAllCells / fineRoiYAllCells exist.');
end

input(['commonF0 pipeline: recomputes dF/F with early-trial common F0; ', ...
    'uses existing fine ROIs (ctrl+c to stop!)']);

f0WinMs_commonF0 = 50;

%% --- F, F0, dF, dF/F0 (_commonF0 field names; workspace *CellNames_commonF0) ---
holoSortedImagingCellNames_commonF0 = cell(nCells, 1);
filtHoloSortedImagingCellNames_commonF0 = cell(nCells, 1);
F0CellNames_commonF0 = cell(nCells, 1);
roiMeanFCellNames_commonF0 = cell(nCells, 1);
bkgrndMeanFCellNames_commonF0 = cell(nCells, 1);
subScalarCellNames_commonF0 = cell(nCells, 1);
roiMeanFCorrectedCellNames_commonF0 = cell(nCells, 1);
globalF0CellNames_commonF0 = cell(nCells, 1);
dFCellNames_commonF0 = cell(nCells, 1);
dFFCellNames_commonF0 = cell(nCells, 1);

for nn = 1:nCells
    F0CellNames_commonF0{nn}                    = ['F0AllTrials_', 'cell', num2str(nn), '_commonF0'];
    roiMeanFCellNames_commonF0{nn}              = ['roiMeanF_', 'cell', num2str(nn), '_commonF0'];
    bkgrndMeanFCellNames_commonF0{nn}           = ['bkgrndMeanF_', 'cell', num2str(nn), '_commonF0'];
    subScalarCellNames_commonF0{nn}             = ['subScalar_', 'cell', num2str(nn), '_commonF0'];
    roiMeanFCorrectedCellNames_commonF0{nn}     = ['roiMeanFCorrected_', 'cell', num2str(nn), '_commonF0'];
    globalF0CellNames_commonF0{nn}              = ['globalF0_', 'cell', num2str(nn), '_commonF0'];
    dFCellNames_commonF0{nn}                    = ['dF_', 'cell', num2str(nn), '_commonF0'];
    dFFCellNames_commonF0{nn}                   = ['dFF', 'cell', num2str(nn), '_commonF0'];
    holoSortedImagingCellNames_commonF0{nn}     = ['holoSortedImagingAllTrials_', 'cell', num2str(nn), '_commonF0'];
    filtHoloSortedImagingCellNames_commonF0{nn} = ['filtHoloSortedImagingAllTrials_', 'cell', num2str(nn), '_commonF0'];

    analysisStruct.(roiMeanFCellNames_commonF0{nn})              = [];
    analysisStruct.(bkgrndMeanFCellNames_commonF0{nn})           = [];
    analysisStruct.(subScalarCellNames_commonF0{nn})             = [];
    analysisStruct.(roiMeanFCorrectedCellNames_commonF0{nn})     = [];
    analysisStruct.(globalF0CellNames_commonF0{nn})              = [];
    analysisStruct.(dFCellNames_commonF0{nn})                    = [];
    analysisStruct.(dFFCellNames_commonF0{nn})                   = [];
    analysisStruct.(F0CellNames_commonF0{nn})                    = cell(nConds, 1);
    analysisStruct.(holoSortedImagingCellNames_commonF0{nn})     = cell(nConds, 1);
    analysisStruct.(filtHoloSortedImagingCellNames_commonF0{nn}) = cell(nConds, 1);

    for cc = 1:nConds
        analysisStruct.(F0CellNames_commonF0{nn}){cc}                    = cell(nHolos(cc), 1);
        analysisStruct.(holoSortedImagingCellNames_commonF0{nn}){cc}     = cell(nHolos(cc), 1);
        analysisStruct.(filtHoloSortedImagingCellNames_commonF0{nn}){cc} = cell(nHolos(cc), 1);
    end
end

startTimeImaging_commonF0 = floor(startTime * imagingFreq);

counter_commonF0 = 0;
for tt = 1:nTrials
    counter_commonF0 = counter_commonF0 + 1;
    disp(['commonF0 pipeline — trial ', num2str(counter_commonF0)]);

    rawName = ImgfolderContents(imagesIndex(tt)).name;
    [~, baseName, ~] = fileparts(rawName);
    mcName = [baseName, '_mc.tif'];
    currImgPath = fullfile(mcTiffFolder, mcName);

    t = Tiff(currImgPath, 'r');
    n = 1;
    while true
        try
            t.setDirectory(n);
            n = n + 1;
        catch
            n = n - 1;
            break
        end
    end
    numFrames = n;

    t.setDirectory(1);
    firstFrame = t.read();
    [H, W] = size(firstFrame);

    imageStack = zeros(H, W, numFrames, 'like', firstFrame);
    imageStack(:, :, 1) = firstFrame;
    for k = 2:numFrames
        t.setDirectory(k);
        imageStack(:, :, k) = t.read();
    end
    t.close();

    badRowMaskMc = false(H, numFrames);
    if laserArtifactMcSecondSweepForDff
        badRowsMatPath = fullfile(mcTiffFolder, [baseName, '_mc_badRows.mat']);
        if exist(badRowsMatPath, 'file') == 2
            S = load(badRowsMatPath, 'badRowMask');
            badRowMaskMc = logical(S.badRowMask);
        else
            badRowMaskMc = VoltImg_laserRowArtifact_badRowMaskStack(single(imageStack), ...
                laserArtifactGateColFirst, laserArtifactGateColLast, ...
                laserArtifactThreshMode, laserArtifactThreshParam);
        end
        if ~isequal(size(badRowMaskMc), [H, numFrames])
            warning('VoltImg:laserBadRowsSize', ...
                'commonF0: badRowMask size mismatch trial %d; ignoring mask.', tt);
            badRowMaskMc = false(H, numFrames);
        end
    end

    meanImgThisTrial = mean(single(imageStack), 3);
    meanImgThisTrialDouble = im2double(meanImgThisTrial);

    allTrialRoiMask = false(size(meanImgThisTrial));
    for nn = 1:nCells
        if ~isempty(fineRoiXAllCells{nn}{tt})
            trialInd = sub2ind(size(allTrialRoiMask), fineRoiXAllCells{nn}{tt}, fineRoiYAllCells{nn}{tt});
            allTrialRoiMask(trialInd) = true;
        end
    end

    for nn = 1:nCells
        if laserArtifactMcSecondSweepForDff
            roiMeanF = VoltImg_roiMeanPerFrameExcludeBadRows(imageStack, ...
                fineRoiXAllCells{nn}{tt}, fineRoiYAllCells{nn}{tt}, badRowMaskMc);
        else
            rawWholeRoiF = imageStack(fineRoiXAllCells{nn}{tt}, fineRoiYAllCells{nn}{tt}, :);
            roiMeanF = zeros(size(rawWholeRoiF, 3), 1);
            for ff = 1:size(rawWholeRoiF, 3)
                roiMeanF(ff, 1) = mean(mean(rawWholeRoiF(:, :, ff)));
            end
        end

        innerBuffer = 2;
        ringWidth   = 3;
        minArea     = 50;

        roiMaskThisCell = false(size(meanImgThisTrial));
        roiIndThisCell = sub2ind(size(roiMaskThisCell), fineRoiXAllCells{nn}{tt}, fineRoiYAllCells{nn}{tt});
        roiMaskThisCell(roiIndThisCell) = true;

        innerSelect = imdilate(roiMaskThisCell, strel('disk', innerBuffer));
        outerSelect = imdilate(roiMaskThisCell, strel('disk', innerBuffer + ringWidth));
        backgroundRing = outerSelect & ~innerSelect;
        backgroundRing = backgroundRing & ~allTrialRoiMask;

        valsBk = meanImgThisTrialDouble(backgroundRing);
        if ~isempty(valsBk)
            brightCut = prctile(valsBk, 95);
            ringClean = backgroundRing & (meanImgThisTrialDouble <= brightCut);
        else
            ringClean = backgroundRing;
        end

        ringClean = bwareaopen(ringClean, 7);
        if nnz(ringClean) < minArea
            ringClean = backgroundRing;
        end

        if nnz(ringClean) < 1
            ringGlobalMask = false(size(meanImgThisTrial));
            if ~isempty(bkgrndRoiXAllCells{nn})
                globalInd = sub2ind(size(ringGlobalMask), bkgrndRoiXAllCells{nn}, bkgrndRoiYAllCells{nn});
                ringGlobalMask(globalInd) = true;
            end
            ringClean = ringGlobalMask & ~allTrialRoiMask;
        end

        [bkgrndRoiXTrial, bkgrndRoiYTrial] = find(ringClean);

        if laserArtifactMcSecondSweepForDff
            bkgrndMeanF = VoltImg_roiMeanPerFrameExcludeBadRows(imageStack, ...
                bkgrndRoiXTrial, bkgrndRoiYTrial, badRowMaskMc);
        else
            rawWholeBkgrndF = imageStack(bkgrndRoiXTrial, bkgrndRoiYTrial, :);
            bkgrndMeanF = zeros(size(rawWholeBkgrndF, 3), 1);
            for ff = 1:size(rawWholeBkgrndF, 3)
                bkgrndMeanF(ff, 1) = mean(mean(rawWholeBkgrndF(:, :, ff)));
            end
        end

        baselineIndices = 1:startTimeImaging_commonF0;
        bFit = robustfit(bkgrndMeanF(baselineIndices), roiMeanF(baselineIndices));
        alphaScalar = bFit(2);
        alphaScalar = min(max(alphaScalar, 0), 1);
        if alphaScalar > 0.8
            alphaScalar = 0.8;
        end
        alphaScalar = 0.85;

        roiMeanFCorrected = roiMeanF - alphaScalar * bkgrndMeanF;

        if isempty((voltMapping.outParams.sequenceThisTrial{tt}))
            voltMapping.outParams.sequenceThisTrial{tt} = zeroDummySequence;
        end

        cutOffFreqIm = 40;
        [bIm, aIm] = butter(4, cutOffFreqIm / (imagingFreq / 2));
        if ~ismember(tt, excludeTrials)
            roiMeanFCorrectedFilt = filter(bIm, aIm, roiMeanFCorrected(:));
        else
            roiMeanFCorrectedFilt = [];
        end

        f0Trial = NaN;
        f0FiltTrial = NaN;
        if ~ismember(tt, excludeTrials)
            preEnd = min(startTimeImaging_commonF0, numel(roiMeanFCorrected));
            preEnd = max(preEnd, 1);
            preTrace = roiMeanFCorrected(1:preEnd);
            preTrace = preTrace(:);
            Lpre = numel(preTrace);
            winFramesF0 = max(2, round(f0WinMs_commonF0 / 1000 * imagingFreq));
            w = min(winFramesF0, Lpre);
            if Lpre >= 2 && w >= 2 && Lpre >= w
                nWin = Lpre - w + 1;
                winVars = zeros(nWin, 1);
                for sw = 1:nWin
                    winVars(sw) = var(preTrace(sw:sw+w-1), 1);
                end
                [~, iw] = min(winVars);
                f0Trial = mean(preTrace(iw:iw+w-1));
                preFilt = roiMeanFCorrectedFilt(1:preEnd);
                preFilt = preFilt(:);
                f0FiltTrial = mean(preFilt(iw:iw+w-1));
            elseif Lpre >= 1
                f0Trial = mean(preTrace);
                f0FiltTrial = mean(roiMeanFCorrectedFilt(1:preEnd));
            end
            if ~(f0Trial > 0) || ~isfinite(f0Trial)
                warning('VoltImg:f0TrialNonPositive_commonF0', ...
                    'Trial %d cell %d: commonF0 F0 non-positive or non-finite.', tt, nn);
            end
        end

        holoSeqThisTrial = (unique(voltMapping.outParams.sequenceThisTrial{tt}, 'stable') ...
            - min(unique(voltMapping.outParams.sequenceThisTrial{tt}, 'stable')) + 1)';

        for hh = 1:nHolos(voltMapping.trialCond(tt, 1))
            if isempty(voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)})
                voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)} = voltMapping.outParams.firstStimTimes{1, 2};
            end

            iHoloLo = floor((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh) - preStimWindow / 1000) * imagingFreq);
            iHoloHi = ceil((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh) - preStimWindow / 1000) * imagingFreq) + ...
                ceil((ipi * nPulses + (preStimWindow + postStimWindow)) / 1000 * imagingFreq);

            if ~ismember(tt, excludeTrials)
                roiFCorrectedThisHolo = roiMeanFCorrected(iHoloLo:iHoloHi);
                f0ThisHolo = f0Trial;
                dFThisHolo = roiFCorrectedThisHolo - f0ThisHolo;
                dFFThisHolo = dFThisHolo / f0ThisHolo;
                if UpOrDown == '2'
                    dFFThisHolo = -dFFThisHolo;
                elseif UpOrDown == '1'
                    dFFThisHolo = dFFThisHolo;
                end

                roiFCorrectedThisHoloFilt = roiMeanFCorrectedFilt(iHoloLo:iHoloHi);
                f0FiltThisHolo = f0FiltTrial;
                dFFiltThisHolo = (roiFCorrectedThisHoloFilt - f0FiltThisHolo) / f0FiltThisHolo;
                if UpOrDown == '2'
                    filtdffThisHolo = -dFFiltThisHolo;
                elseif UpOrDown == '1'
                    filtdffThisHolo = dFFiltThisHolo;
                end
            end

            if ismember(tt, excludeTrials)
                analysisStruct.(F0CellNames_commonF0{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = ...
                    [analysisStruct.(F0CellNames_commonF0{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, NaN];
                analysisStruct.(holoSortedImagingCellNames_commonF0{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = ...
                    [analysisStruct.(holoSortedImagingCellNames_commonF0{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, ...
                    NaN(ceil((ipi * nPulses + (preStimWindow + postStimWindow)) / 1000 * imagingFreq) + 2, 1)];
                analysisStruct.(filtHoloSortedImagingCellNames_commonF0{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = ...
                    [analysisStruct.(filtHoloSortedImagingCellNames_commonF0{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, ...
                    NaN(ceil((ipi * nPulses + (preStimWindow + postStimWindow)) / 1000 * imagingFreq) + 2, 1)];
            else
                analysisStruct.(F0CellNames_commonF0{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = ...
                    [analysisStruct.(F0CellNames_commonF0{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, f0ThisHolo];
                analysisStruct.(holoSortedImagingCellNames_commonF0{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = ...
                    [analysisStruct.(holoSortedImagingCellNames_commonF0{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, dFFThisHolo];
                analysisStruct.(filtHoloSortedImagingCellNames_commonF0{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = ...
                    [analysisStruct.(filtHoloSortedImagingCellNames_commonF0{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, filtdffThisHolo];
            end
        end

        if ismember(tt, excludeTrials)
            analysisStruct.(roiMeanFCellNames_commonF0{nn})(:, tt)          = NaN(numFrames, 1);
            analysisStruct.(bkgrndMeanFCellNames_commonF0{nn})(:, tt)       = NaN(numFrames, 1);
            analysisStruct.(subScalarCellNames_commonF0{nn})(tt, 1)         = NaN;
            analysisStruct.(roiMeanFCorrectedCellNames_commonF0{nn})(:, tt) = NaN(numFrames, 1);
        else
            analysisStruct.(roiMeanFCellNames_commonF0{nn})(:, tt)          = roiMeanF;
            analysisStruct.(bkgrndMeanFCellNames_commonF0{nn})(:, tt)       = bkgrndMeanF;
            analysisStruct.(subScalarCellNames_commonF0{nn})(tt, 1)         = alphaScalar;
            analysisStruct.(roiMeanFCorrectedCellNames_commonF0{nn})(:, tt) = roiMeanFCorrected;
        end
    end
end

%% Calculate mean response (and CI) for each hologram — _commonF0
holoSortedMeanCellNames_commonF0 = cell(nCells, 1);
filtHoloSortedMeanCellNames_commonF0 = cell(nCells, 1);
for nn = 1:nCells
    holoSortedMeanCellNames_commonF0{nn}     = ['holoSortedImagingMean_', 'cell', num2str(nn), '_commonF0'];
    filtHoloSortedMeanCellNames_commonF0{nn} = ['filtHoloSortedImagingMean_', 'cell', num2str(nn), '_commonF0'];
    analysisStruct.(holoSortedMeanCellNames_commonF0{nn})     = cell(nConds, 1);
    analysisStruct.(filtHoloSortedMeanCellNames_commonF0{nn}) = cell(nConds, 1);
    for cc = 1:nConds
        analysisStruct.(holoSortedMeanCellNames_commonF0{nn}){cc}     = cell(nHolos(cc), 1);
        analysisStruct.(filtHoloSortedMeanCellNames_commonF0{nn}){cc} = cell(nHolos(cc), 1);
    end
end

for nn = 1:nCells
    for cc = 1:nConds
        for hh = 1:nHolos(cc)
            analysisStruct.(holoSortedMeanCellNames_commonF0{nn}){cc}{hh} = ...
                nanmean(analysisStruct.(holoSortedImagingCellNames_commonF0{nn}){cc}{hh}, 2);
            analysisStruct.(filtHoloSortedMeanCellNames_commonF0{nn}){cc}{hh} = ...
                nanmean(analysisStruct.(filtHoloSortedImagingCellNames_commonF0{nn}){cc}{hh}, 2);
        end
    end
end

CIDffAllCondsCellNames_commonF0 = cell(nCells, 1);
filtCIDffAllCondsCellNames_commonF0 = cell(nCells, 1);
for nn = 1:nCells
    CIDffAllCondsCellNames_commonF0{nn}     = ['CIDffAllConds_', 'cell', num2str(nn), '_commonF0'];
    filtCIDffAllCondsCellNames_commonF0{nn} = ['filtCIDffAllConds_', 'cell', num2str(nn), '_commonF0'];
    analysisStruct.(CIDffAllCondsCellNames_commonF0{nn})     = cell(nConds, 1);
    analysisStruct.(filtCIDffAllCondsCellNames_commonF0{nn}) = cell(nConds, 1);
end

for nn = 1:nCells
    for cc = 1:nConds
        for hh = 1:nHolos(cc)
            confidence_level_commonF0 = 0.95;
            means_cf     = nanmean(analysisStruct.(holoSortedImagingCellNames_commonF0{nn}){cc}{hh, 1}, 2);
            filtMeans_cf = nanmean(analysisStruct.(filtHoloSortedImagingCellNames_commonF0{nn}){cc}{hh, 1}, 2);
            std_errors_cf     = std(analysisStruct.(holoSortedImagingCellNames_commonF0{nn}){cc}{hh, 1}, 0, 2, "omitnan") ...
                / sqrt(size(analysisStruct.(holoSortedImagingCellNames_commonF0{nn}){cc}{hh, 1}, 2));
            filtStd_errors_cf = std(analysisStruct.(filtHoloSortedImagingCellNames_commonF0{nn}){cc}{hh, 1}, 0, 2, "omitnan") ...
                / sqrt(size(analysisStruct.(filtHoloSortedImagingCellNames_commonF0{nn}){cc}{hh, 1}, 2));

            t_score_cf     = tinv((1 + confidence_level_commonF0) / 2, size(analysisStruct.(holoSortedImagingCellNames_commonF0{nn}){cc}{hh, 1}, 2) - 1);
            filtT_score_cf = tinv((1 + confidence_level_commonF0) / 2, size(analysisStruct.(filtHoloSortedImagingCellNames_commonF0{nn}){cc}{hh, 1}, 2) - 1);
            margin_of_error_cf     = t_score_cf * std_errors_cf;
            filtMargin_of_error_cf = filtT_score_cf * filtStd_errors_cf;
            lower_bounds_cf     = means_cf - margin_of_error_cf;
            filtLower_bounds_cf = filtMeans_cf - filtMargin_of_error_cf;
            upper_bounds_cf     = means_cf + margin_of_error_cf;
            filtUpper_bounds_cf = filtMeans_cf + filtMargin_of_error_cf;
            if UpOrDown == '2'
                analysisStruct.(CIDffAllCondsCellNames_commonF0{nn}){cc}{hh, 1}     = [lower_bounds_cf, upper_bounds_cf];
                analysisStruct.(filtCIDffAllCondsCellNames_commonF0{nn}){cc}{hh, 1} = [filtLower_bounds_cf, filtUpper_bounds_cf];
            elseif UpOrDown == '1'
                analysisStruct.(CIDffAllCondsCellNames_commonF0{nn}){cc}{hh, 1}     = [-lower_bounds_cf, -upper_bounds_cf];
                analysisStruct.(filtCIDffAllCondsCellNames_commonF0{nn}){cc}{hh, 1} = [-filtLower_bounds_cf, -filtUpper_bounds_cf];
            end
        end
    end
end

voltMapping.startTimeImaging_commonF0 = startTimeImaging_commonF0;
voltMapping.f0WinMs_commonF0 = f0WinMs_commonF0;

%% Trial excluder logic — _commonF0 (inlined; mirrors VoltImg_mapping_analysis_MultiCell_trialExcluder.m)
stdImagingAllTrialsCellNames_commonF0 = cell(nCells, 1);
stdFiltImagingAllTrialsCellNames_commonF0 = cell(nCells, 1);
for nn = 1:nCells
    stdImagingAllTrialsCellNames_commonF0{nn}     = ['stdImagingAllTrialsCellNames_', 'cell', num2str(nn), '_commonF0'];
    stdFiltImagingAllTrialsCellNames_commonF0{nn} = ['stdFiltImagingAllTrialsCellNames_', 'cell', num2str(nn), '_commonF0'];
    analysisStruct.(stdImagingAllTrialsCellNames_commonF0{nn})     = cell(nConds, 1);
    analysisStruct.(stdFiltImagingAllTrialsCellNames_commonF0{nn}) = cell(nConds, 1);
    for cc = 1:nConds
        for hh = 1:nHolos(cc)
            analysisStruct.(stdImagingAllTrialsCellNames_commonF0{nn}){cc}{hh, 1} = ...
                std(analysisStruct.(holoSortedImagingCellNames_commonF0{nn}){cc}{hh}(:));
            analysisStruct.(stdFiltImagingAllTrialsCellNames_commonF0{nn}){cc}{hh, 1} = ...
                std(analysisStruct.(filtHoloSortedImagingCellNames_commonF0{nn}){cc}{hh}(:));
        end
    end
end

exclHoloSortedImagingAllTrialsCellNames_commonF0 = cell(nCells, 1);
exclFiltHoloSortedImagingAllTrialsCellNames_commonF0 = cell(nCells, 1);
for nn = 1:nCells
    exclHoloSortedImagingAllTrialsCellNames_commonF0{nn}     = ['exclHoloSortedImagingAllTrials_', 'cell', num2str(nn), '_commonF0'];
    exclFiltHoloSortedImagingAllTrialsCellNames_commonF0{nn} = ['exclFiltHoloSortedImagingAllTrials_', 'cell', num2str(nn), '_commonF0'];
    analysisStruct.(exclHoloSortedImagingAllTrialsCellNames_commonF0{nn})     = cell(nConds, 1);
    analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames_commonF0{nn}) = cell(nConds, 1);
    for cc = 1:nConds
        analysisStruct.(exclHoloSortedImagingAllTrialsCellNames_commonF0{nn}){cc}     = cell(nHolos(cc), 1);
        analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames_commonF0{nn}){cc} = cell(nHolos(cc), 1);
    end
end

for nn = 1:nCells
    for cc = 1:nConds
        for hh = 1:nHolos(cc)
            for tt_ex = 1:size(analysisStruct.(holoSortedImagingCellNames_commonF0{nn}){cc}{hh}, 2)
                if any(analysisStruct.(holoSortedImagingCellNames_commonF0{nn}){cc}{hh}(:, tt_ex) ...
                        < -2.5 * analysisStruct.(stdImagingAllTrialsCellNames_commonF0{nn}){cc}{hh, 1})
                    analysisStruct.(exclHoloSortedImagingAllTrialsCellNames_commonF0{nn}){cc}{hh, 1}(:, tt_ex) = ...
                        nan(size(analysisStruct.(holoSortedImagingCellNames_commonF0{nn}){cc}{hh}(:, tt_ex), 1), 1);
                else
                    analysisStruct.(exclHoloSortedImagingAllTrialsCellNames_commonF0{nn}){cc}{hh, 1}(:, tt_ex) = ...
                        analysisStruct.(holoSortedImagingCellNames_commonF0{nn}){cc}{hh}(:, tt_ex);
                end

                if any(analysisStruct.(filtHoloSortedImagingCellNames_commonF0{nn}){cc}{hh}(:, tt_ex) ...
                        < -2.5 * analysisStruct.(stdFiltImagingAllTrialsCellNames_commonF0{nn}){cc}{hh, 1})
                    analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames_commonF0{nn}){cc}{hh, 1}(:, tt_ex) = ...
                        nan(size(analysisStruct.(filtHoloSortedImagingCellNames_commonF0{nn}){cc}{hh}(:, tt_ex), 1), 1);
                else
                    analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames_commonF0{nn}){cc}{hh, 1}(:, tt_ex) = ...
                        analysisStruct.(filtHoloSortedImagingCellNames_commonF0{nn}){cc}{hh}(:, tt_ex);
                end
            end
        end
    end
end

for nn = 1:nCells
    for cc = 1:nConds
        for hh = 1:nHolos(cc)
            for tt_ex = 1:size(analysisStruct.(exclHoloSortedImagingAllTrialsCellNames_commonF0{nn}){cc}{hh}, 2)
                if any(isnan(analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames_commonF0{nn}){cc}{hh}(:, tt_ex)))
                    continue
                end
                [~, maxImagingIndex] = max(analysisStruct.(exclHoloSortedImagingAllTrialsCellNames_commonF0{nn}){cc}{hh}(:, tt_ex));
                [~, maxFiltImagingIndex] = max(analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames_commonF0{nn}){cc}{hh}(:, tt_ex));
                if maxImagingIndex > 45
                    analysisStruct.(exclHoloSortedImagingAllTrialsCellNames_commonF0{nn}){cc}{hh}(:, tt_ex) = ...
                        nan(size(analysisStruct.(exclHoloSortedImagingAllTrialsCellNames_commonF0{nn}){cc}{hh}(:, tt_ex), 1), 1);
                end
                if maxFiltImagingIndex > 45
                    analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames_commonF0{nn}){cc}{hh}(:, tt_ex) = ...
                        nan(size(analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames_commonF0{nn}){cc}{hh}(:, tt_ex), 1), 1);
                end
            end
        end
    end
end

exclHoloSortedImagingMeanCellNames_commonF0 = cell(nCells, 1);
exclFiltHoloSortedImagingMeanCellNames_commonF0 = cell(nCells, 1);
for nn = 1:nCells
    exclHoloSortedImagingMeanCellNames_commonF0{nn}     = ['exclHoloSortedImagingMean_', 'cell', num2str(nn), '_commonF0'];
    exclFiltHoloSortedImagingMeanCellNames_commonF0{nn} = ['exclFiltHoloSortedImagingMean_', 'cell', num2str(nn), '_commonF0'];
    analysisStruct.(exclHoloSortedImagingMeanCellNames_commonF0{nn})     = cell(nConds, 1);
    analysisStruct.(exclFiltHoloSortedImagingMeanCellNames_commonF0{nn}) = cell(nConds, 1);
    for cc = 1:nConds
        analysisStruct.(exclHoloSortedImagingMeanCellNames_commonF0{nn}){cc}     = cell(nHolos(cc), 1);
        analysisStruct.(exclFiltHoloSortedImagingMeanCellNames_commonF0{nn}){cc} = cell(nHolos(cc), 1);
    end
end

for nn = 1:nCells
    for cc = 1:nConds
        for hh = 1:nHolos(cc)
            analysisStruct.(exclHoloSortedImagingMeanCellNames_commonF0{nn}){cc}{hh} = ...
                nanmean(analysisStruct.(exclHoloSortedImagingAllTrialsCellNames_commonF0{nn}){cc}{hh}, 2);
            analysisStruct.(exclFiltHoloSortedImagingMeanCellNames_commonF0{nn}){cc}{hh} = ...
                nanmean(analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames_commonF0{nn}){cc}{hh}, 2);
        end
    end
end

exclCIDffAllCondsCellNames_commonF0 = cell(nCells, 1);
exclFiltCIDffAllCondsCellNames_commonF0 = cell(nCells, 1);
for nn = 1:nCells
    exclCIDffAllCondsCellNames_commonF0{nn}     = ['exclCIDffAllConds_', 'cell', num2str(nn), '_commonF0'];
    exclFiltCIDffAllCondsCellNames_commonF0{nn} = ['exclFiltCIDffAllConds_', 'cell', num2str(nn), '_commonF0'];
    analysisStruct.(exclCIDffAllCondsCellNames_commonF0{nn})     = cell(nConds, 1);
    analysisStruct.(exclFiltCIDffAllCondsCellNames_commonF0{nn}) = cell(nConds, 1);
end

for nn = 1:nCells
    for cc = 1:nConds
        for hh = 1:nHolos(cc)
            confidence_level_excl_cf = 0.95;
            means_ex_cf     = nanmean(analysisStruct.(exclHoloSortedImagingAllTrialsCellNames_commonF0{nn}){cc}{hh, 1}, 2);
            filtMeans_ex_cf = nanmean(analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames_commonF0{nn}){cc}{hh, 1}, 2);
            std_errors_ex_cf     = std(analysisStruct.(exclHoloSortedImagingAllTrialsCellNames_commonF0{nn}){cc}{hh, 1}, 0, 2, "omitnan") ...
                / sqrt(size(analysisStruct.(exclHoloSortedImagingAllTrialsCellNames_commonF0{nn}){cc}{hh, 1}, 2));
            filtStd_errors_ex_cf = std(analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames_commonF0{nn}){cc}{hh, 1}, 0, 2, "omitnan") ...
                / sqrt(size(analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames_commonF0{nn}){cc}{hh, 1}, 2));

            t_score_ex_cf     = tinv((1 + confidence_level_excl_cf) / 2, size(analysisStruct.(exclHoloSortedImagingAllTrialsCellNames_commonF0{nn}){cc}{hh, 1}, 2) - 1);
            filtT_score_ex_cf = tinv((1 + confidence_level_excl_cf) / 2, size(analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames_commonF0{nn}){cc}{hh, 1}, 2) - 1);
            margin_of_error_ex_cf     = t_score_ex_cf * std_errors_ex_cf;
            filtMargin_of_error_ex_cf = filtT_score_ex_cf * filtStd_errors_ex_cf;
            lower_bounds_ex_cf     = means_ex_cf - margin_of_error_ex_cf;
            filtLower_bounds_ex_cf = filtMeans_ex_cf - filtMargin_of_error_ex_cf;
            upper_bounds_ex_cf     = means_ex_cf + margin_of_error_ex_cf;
            filtUpper_bounds_ex_cf = filtMeans_ex_cf + filtMargin_of_error_ex_cf;
            if UpOrDown == '2'
                analysisStruct.(exclCIDffAllCondsCellNames_commonF0{nn}){cc}{hh, 1}     = [lower_bounds_ex_cf, upper_bounds_ex_cf];
                analysisStruct.(exclFiltCIDffAllCondsCellNames_commonF0{nn}){cc}{hh, 1} = [filtLower_bounds_ex_cf, filtUpper_bounds_ex_cf];
            elseif UpOrDown == '1'
                analysisStruct.(exclCIDffAllCondsCellNames_commonF0{nn}){cc}{hh, 1}     = [-lower_bounds_ex_cf, -upper_bounds_ex_cf];
                analysisStruct.(exclFiltCIDffAllCondsCellNames_commonF0{nn}){cc}{hh, 1} = [-filtLower_bounds_ex_cf, -filtUpper_bounds_ex_cf];
            end
        end
    end
end

%% Reorganize voltMapping structs by cells — _commonF0 (adds fields; does not remove originals)
if ~exist('cellID', 'var')
    error('VoltImg:commonF0NeedsCellID', 'cellID not found; run original reorganize block first or define cellID.');
end

voltMapping.holoSortedImagingCellNames_commonF0                  = holoSortedImagingCellNames_commonF0;
voltMapping.filtHoloSortedImagingCellNames_commonF0              = filtHoloSortedImagingCellNames_commonF0;
voltMapping.holoSortedMeanCellNames_commonF0                     = holoSortedMeanCellNames_commonF0;
voltMapping.filtHoloSortedMeanCellNames_commonF0                 = filtHoloSortedMeanCellNames_commonF0;
voltMapping.CIDffAllCondsCellNames_commonF0                      = CIDffAllCondsCellNames_commonF0;
voltMapping.filtCIDffAllCondsCellNames_commonF0                  = filtCIDffAllCondsCellNames_commonF0;
voltMapping.stdImagingAllTrialsCellNames_commonF0                = stdImagingAllTrialsCellNames_commonF0;
voltMapping.stdFiltImagingAllTrialsCellNames_commonF0            = stdFiltImagingAllTrialsCellNames_commonF0;
voltMapping.exclHoloSortedImagingAllTrialsCellNames_commonF0     = exclHoloSortedImagingAllTrialsCellNames_commonF0;
voltMapping.exclFiltHoloSortedImagingAllTrialsCellNames_commonF0 = exclFiltHoloSortedImagingAllTrialsCellNames_commonF0;
voltMapping.exclHoloSortedImagingMeanCellNames_commonF0          = exclHoloSortedImagingMeanCellNames_commonF0;
voltMapping.exclFiltHoloSortedImagingMeanCellNames_commonF0      = exclFiltHoloSortedImagingMeanCellNames_commonF0;
voltMapping.exclCIDffAllCondsCellNames_commonF0                  = exclCIDffAllCondsCellNames_commonF0;
voltMapping.exclFiltCIDffAllCondsCellNames_commonF0             = exclFiltCIDffAllCondsCellNames_commonF0;
voltMapping.F0CellNames_commonF0                                 = F0CellNames_commonF0;
voltMapping.roiMeanFCellNames_commonF0                           = roiMeanFCellNames_commonF0;
voltMapping.bkgrndMeanFCellNames_commonF0                        = bkgrndMeanFCellNames_commonF0;
voltMapping.subScalarCellNames_commonF0                         = subScalarCellNames_commonF0;
voltMapping.roiMeanFCorrectedCellNames_commonF0                  = roiMeanFCorrectedCellNames_commonF0;

for nn = 1:nCells
    structname = ['voltMapping.', 'Cell', num2str(nn)];
    eval([structname '.holoSortedImagingAllTrials_commonF0 = analysisStruct.(holoSortedImagingCellNames_commonF0{nn});']);
    eval([structname '.filtHoloSortedImagingAllTrials_commonF0 = analysisStruct.(filtHoloSortedImagingCellNames_commonF0{nn});']);
    eval([structname '.holoSortedImagingMean_commonF0 = analysisStruct.(holoSortedMeanCellNames_commonF0{nn});']);
    eval([structname '.filtHoloSortedImagingMean_commonF0 = analysisStruct.(filtHoloSortedMeanCellNames_commonF0{nn});']);
    eval([structname '.CIDffAllConds_commonF0 = analysisStruct.(CIDffAllCondsCellNames_commonF0{nn});']);
    eval([structname '.filtCIDffAllConds_commonF0 = analysisStruct.(filtCIDffAllCondsCellNames_commonF0{nn});']);
    eval([structname '.stdImagingAllTrialsCellNames_commonF0 = analysisStruct.(stdImagingAllTrialsCellNames_commonF0{nn});']);
    eval([structname '.stdFiltImagingAllTrialsCellNames_commonF0 = analysisStruct.(stdFiltImagingAllTrialsCellNames_commonF0{nn});']);
    eval([structname '.exclHoloSortedImagingAllTrials_commonF0 = analysisStruct.(exclHoloSortedImagingAllTrialsCellNames_commonF0{nn});']);
    eval([structname '.exclFiltHoloSortedImagingAllTrials_commonF0 = analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames_commonF0{nn});']);
    eval([structname '.exclHoloSortedImagingMean_commonF0 = analysisStruct.(exclHoloSortedImagingMeanCellNames_commonF0{nn});']);
    eval([structname '.exclFiltHoloSortedImagingMean_commonF0 = analysisStruct.(exclFiltHoloSortedImagingMeanCellNames_commonF0{nn});']);
    eval([structname '.exclCIDffAllConds_commonF0 = analysisStruct.(exclCIDffAllCondsCellNames_commonF0{nn});']);
    eval([structname '.exclFiltCIDffAllConds_commonF0 = analysisStruct.(exclFiltCIDffAllCondsCellNames_commonF0{nn});']);
    eval([structname '.F0AllTrials_commonF0 = analysisStruct.(F0CellNames_commonF0{nn});']);
    eval([structname '.roiMeanF_commonF0 = analysisStruct.(roiMeanFCellNames_commonF0{nn});']);
    eval([structname '.bkgrndMeanF_commonF0 = analysisStruct.(bkgrndMeanFCellNames_commonF0{nn});']);
    eval([structname '.subScalar_commonF0 = analysisStruct.(subScalarCellNames_commonF0{nn});']);
    eval([structname '.roiMeanFCorrected_commonF0 = analysisStruct.(roiMeanFCorrectedCellNames_commonF0{nn});']);
end

%% Crash-safe checkpoint: dF/F + trial excluder — _commonF0
if exist('UpOrDown', 'var')
    voltMapping.UpOrDown_commonF0 = UpOrDown;
end
if exist('ePhysAvail', 'var')
    voltMapping.ePhysAvail_commonF0 = ePhysAvail;
end
if ~exist('confidence_level', 'var') || isempty(confidence_level)
    confidence_level = 0.95;
end
checkpointAfterDffFile_commonF0 = fullfile(saveDirectory, 'checkpoint_after_dff_calculation_commonF0.mat');
analysisCheckpointTimestamp_dff_commonF0 = datetime('now');
varsCheckpointAfterDff_commonF0 = { ...
    'voltMapping', 'analysisStruct', ...
    'F0CellNames_commonF0', 'roiMeanFCellNames_commonF0', 'bkgrndMeanFCellNames_commonF0', 'subScalarCellNames_commonF0', ...
    'roiMeanFCorrectedCellNames_commonF0', 'globalF0CellNames_commonF0', 'dFCellNames_commonF0', 'dFFCellNames_commonF0', ...
    'holoSortedImagingCellNames_commonF0', 'filtHoloSortedImagingCellNames_commonF0', ...
    'holoSortedMeanCellNames_commonF0', 'filtHoloSortedMeanCellNames_commonF0', ...
    'CIDffAllCondsCellNames_commonF0', 'filtCIDffAllCondsCellNames_commonF0', ...
    'stdImagingAllTrialsCellNames_commonF0', 'stdFiltImagingAllTrialsCellNames_commonF0', ...
    'exclHoloSortedImagingAllTrialsCellNames_commonF0', 'exclFiltHoloSortedImagingAllTrialsCellNames_commonF0', ...
    'exclHoloSortedImagingMeanCellNames_commonF0', 'exclFiltHoloSortedImagingMeanCellNames_commonF0', ...
    'exclCIDffAllCondsCellNames_commonF0', 'exclFiltCIDffAllCondsCellNames_commonF0', ...
    'mouseID', 'ePhysAvail', ...
    'nCells', 'nTrials', 'nConds', 'nHolos', 'imagingFreq', 'Fs', 'ipi', 'nPulses', ...
    'preStimWindow', 'postStimWindow', 'startTime', 'UpOrDown', 'excludeTrials', ...
    'powers', 'nextHoloDelay', 'pulseDurs', 'trialTime', ...
    'ImgfolderContents', 'savePath', 'saveDirectory', ...
    'normcorrePath', ...
    'cutOffFreq', 'blp', 'alp', 'cutOffFreqIm', 'bIm', 'aIm', ...
    'vThreshold', 'confidence_level', ...
    'mappingInputs', 'allBkgrndRois', 'nPlateauPowerPts', ...
    'useLaserRowArtifactFilter', 'laserArtifactGateColFirst', 'laserArtifactGateColLast', ...
    'laserArtifactThreshMode', 'laserArtifactThreshParam', 'laserArtifactMcMode', ...
    'mcUseGateColumnsOnly', 'laserArtifactMcSecondSweepForDff', ...
    'cellID', 'analysisCheckpointTimestamp_dff_commonF0', ...
    'startTimeImaging_commonF0', 'f0WinMs_commonF0'};
if exist('zeroDummySequence', 'var')
    varsCheckpointAfterDff_commonF0{end+1} = 'zeroDummySequence';
end
varsCheckpointAfterDff_commonF0 = varsCheckpointAfterDff_commonF0(cellfun(@(n) exist(n, 'var') == 1, varsCheckpointAfterDff_commonF0));
try
    save(checkpointAfterDffFile_commonF0, varsCheckpointAfterDff_commonF0{:}, '-v7.3');
catch
    try
        save(checkpointAfterDffFile_commonF0, varsCheckpointAfterDff_commonF0{:});
    catch ME
        warning('VoltImg:checkpointAfterDff_commonF0', 'Could not write commonF0 checkpoint: %s', ME.message);
    end
end
disp(['Checkpoint saved (commonF0 dF/F + trial excluder): ', checkpointAfterDffFile_commonF0]);

%% Four-panel comparison — _commonF0 pipeline (1×4 subplots; figure ids 40000+)
cellIdx_commonF0 = double(input('which cell number (4-panel mean + CI, commonF0 pipeline)? '));

holoSortedImagingMean_commonF0         = voltMapping.(cellID{cellIdx_commonF0}).holoSortedImagingMean_commonF0;
filtHoloSortedImagingMean_commonF0     = voltMapping.(cellID{cellIdx_commonF0}).filtHoloSortedImagingMean_commonF0;
exclHoloSortedImagingMean_commonF0     = voltMapping.(cellID{cellIdx_commonF0}).exclHoloSortedImagingMean_commonF0;
exclFiltHoloSortedImagingMean_commonF0 = voltMapping.(cellID{cellIdx_commonF0}).exclFiltHoloSortedImagingMean_commonF0;
CIDffAllConds_commonF0                 = voltMapping.(cellID{cellIdx_commonF0}).CIDffAllConds_commonF0;
filtCIDffAllConds_commonF0             = voltMapping.(cellID{cellIdx_commonF0}).filtCIDffAllConds_commonF0;
exclCIDffAllConds_commonF0             = voltMapping.(cellID{cellIdx_commonF0}).exclCIDffAllConds_commonF0;
exclFiltCIDffAllConds_commonF0         = voltMapping.(cellID{cellIdx_commonF0}).exclFiltCIDffAllConds_commonF0;

meanSeries_commonF0 = {holoSortedImagingMean_commonF0, filtHoloSortedImagingMean_commonF0, ...
    exclHoloSortedImagingMean_commonF0, exclFiltHoloSortedImagingMean_commonF0};
ciSeries_commonF0   = {CIDffAllConds_commonF0, filtCIDffAllConds_commonF0, ...
    exclCIDffAllConds_commonF0, exclFiltCIDffAllConds_commonF0};
yScales_commonF0    = [10, 10, 10, 10];
subTitles_commonF0 = {'Mean (common F0)', 'Filt mean (common F0)', 'Excl mean (common F0)', 'Excl filt mean (common F0)'};

FsEphys_cf = voltMapping.daqParams.Fs;
plotEphysWithImaging_commonF0 = exist('ePhysAvail', 'var') && ePhysAvail == 1 ...
    && isfield(voltMapping, 'ephys') && isfield(voltMapping.ephys, 'holoSortedDataMean');

for cc_cf = 1:nConds
    for hh_cf = 1:nHolos(cc_cf)
        figure(40000 + cc_cf * 1000 + hh_cf);
        set(gcf, 'Position', [0, 0, 2000, 420]);
        clf

        nT_cf = numel(holoSortedImagingMean_commonF0{cc_cf}{hh_cf});
        tAxis_cf = linspace(0, nT_cf / imagingFreq, nT_cf);
        pdMs_cf = voltMapping.outParams.pulseDur(:);
        pulseDurMsHere_cf = pdMs_cf(min(cc_cf, numel(pdMs_cf)));
        if pulseDurMsHere_cf <= 0 && ~isempty(pulseDurs)
            pulseDurMsHere_cf = pulseDurs(1);
        end
        pulseDurSecHere_cf = pulseDurMsHere_cf / 1000;

        ephysOnImgGrid_cf = [];
        if plotEphysWithImaging_commonF0 && hh_cf <= size(voltMapping.ephys.holoSortedDataMean{cc_cf}, 2)
            ephysMeanCol_cf = voltMapping.ephys.holoSortedDataMean{cc_cf}(:, hh_cf);
            nEp_cf = numel(ephysMeanCol_cf);
            if nEp_cf > 1
                tEp_cf = linspace(0, nEp_cf / FsEphys_cf, nEp_cf);
                ephysOnImgGrid_cf = interp1(tEp_cf(:), double(ephysMeanCol_cf(:)), tAxis_cf(:), 'linear', 'extrap');
            elseif nEp_cf == 1
                ephysOnImgGrid_cf = repmat(double(ephysMeanCol_cf), numel(tAxis_cf), 1);
            end
        end

        ephysRgb_cf = [0.9, 0.9, 0.9];
        for sp_cf = 1:4
            ax_cf = subplot(1, 4, sp_cf);
            m_cf = meanSeries_commonF0{sp_cf}{cc_cf}{hh_cf} * 10;
            ci_cf = ciSeries_commonF0{sp_cf}{cc_cf}{hh_cf} * 10;
            ys_cf = yScales_commonF0(sp_cf);

            yyaxis(ax_cf, 'left')
            hold(ax_cf, 'on')
            fill(ax_cf, [tAxis_cf, fliplr(tAxis_cf)], [ci_cf(:, 1)' * ys_cf, fliplr(ci_cf(:, 2)' * ys_cf)], ...
                [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
            ylStim_cf = ax_cf.YLim;
            for pulseIdx_cf = 1:length(nPulseCoords)
                tOn_cf = nPulseCoords(pulseIdx_cf) / FsEphys_cf;
                patch(ax_cf, [tOn_cf, tOn_cf + pulseDurSecHere_cf, tOn_cf + pulseDurSecHere_cf, tOn_cf], ...
                    [ylStim_cf(1), ylStim_cf(1), ylStim_cf(2), ylStim_cf(2)], [1, 0, 0], ...
                    'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HitTest', 'off');
            end
            plot(ax_cf, tAxis_cf, ci_cf(:, 1) * ys_cf, '--', 'linewidth', 1, 'color', [0.7, 0.7, 0.7]);
            plot(ax_cf, tAxis_cf, ci_cf(:, 2) * ys_cf, '--', 'linewidth', 1, 'color', [0.7, 0.7, 0.7]);
            imgMeanLine_cf = m_cf * ys_cf;
            plot(ax_cf, tAxis_cf, imgMeanLine_cf, '-', 'linewidth', 2, 'color', 'g');
            ylabel(ax_cf, 'dF/F (%)');

            if plotEphysWithImaging_commonF0 && ~isempty(ephysOnImgGrid_cf)
                imgPk_cf = max(abs(imgMeanLine_cf(:)), [], 'omitnan');
                epPk_cf = max(abs(ephysOnImgGrid_cf(:)), [], 'omitnan');
                yyaxis(ax_cf, 'right')
                hold(ax_cf, 'on')
                scaleEp_cf = nan;
                if isfinite(imgPk_cf) && isfinite(epPk_cf) && epPk_cf > 0 && imgPk_cf > 0
                    scaleEp_cf = imgPk_cf / epPk_cf;
                    ephysScaledRow_cf = ephysOnImgGrid_cf(:).' * scaleEp_cf;
                    plot(ax_cf, tAxis_cf, ephysScaledRow_cf, '-', 'linewidth', 2, 'color', ephysRgb_cf);
                else
                    plot(ax_cf, tAxis_cf, ephysOnImgGrid_cf(:).', '-', 'linewidth', 2, 'color', ephysRgb_cf);
                end
                if isfinite(scaleEp_cf) && scaleEp_cf > 0
                    yRightLbl_cf = sprintf('Rel. Vm (ephys); value / %.3g = mV', scaleEp_cf);
                else
                    yRightLbl_cf = 'Vm ephys (mV, rel. to baseline)';
                end
                ylabel(ax_cf, yRightLbl_cf, 'Color', ephysRgb_cf, 'Interpreter', 'none');
                ax_cf.YAxis(2).Color = ephysRgb_cf;
                ax_cf.YAxis(2).Label.Visible = 'on';
                ax_cf.YAxis(2).Label.FontWeight = 'normal';
            end

            yyaxis(ax_cf, 'left')
            hold(ax_cf, 'off')
            xlabel(ax_cf, 'Time (s)');
            title(ax_cf, subTitles_commonF0{sp_cf});
        end
        if plotEphysWithImaging_commonF0 && ~isempty(ephysOnImgGrid_cf)
            sgtitle(sprintf('Cell %d — cond %d, holo %d [commonF0] (ephys right axis)', cellIdx_commonF0, cc_cf, hh_cf));
        else
            sgtitle(sprintf('Cell %d — cond %d, holo %d [commonF0] (4-panel)', cellIdx_commonF0, cc_cf, hh_cf));
        end
        pause
    end
end