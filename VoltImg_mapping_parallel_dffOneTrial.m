function pack = VoltImg_mapping_parallel_dffOneTrial(tt, ctx)
%VOLTIMG_MAPPING_PARALLEL_DFFONETRIAL  Per-trial fine ROI + F traces + holo-sorted dF/F (one trial).

pack = struct();
pack.tt = tt;
pack.numFrames = 0;
pack.cc = ctx.voltMapping.trialCond(tt, 1);
pack.nCells = ctx.nCells;
pack.isExcluded = ismember(tt, ctx.excludeTrials);
cc = pack.cc;
nHolosThis = ctx.nHolos(cc);

rawName = ctx.ImgfolderContents(ctx.imagesIndex(tt)).name;
[~, baseName, ~] = fileparts(rawName);
mcName = [baseName, '_mc.tif'];
currImgPath = fullfile(ctx.mcTiffFolder, mcName);

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

pack.numFrames = numFrames;

badRowMaskMc = false(H, numFrames);
if ctx.laserArtifactMcSecondSweepForDff
    badRowsMatPath = fullfile(ctx.mcTiffFolder, [baseName, '_mc_badRows.mat']);
    if exist(badRowsMatPath, 'file') == 2
        S = load(badRowsMatPath, 'badRowMask');
        badRowMaskMc = logical(S.badRowMask);
    else
        badRowMaskMc = VoltImg_laserRowArtifact_badRowMaskStack(single(imageStack), ...
            ctx.laserArtifactGateColFirst, ctx.laserArtifactGateColLast, ...
            ctx.laserArtifactThreshMode, ctx.laserArtifactThreshParam);
    end
    if ~isequal(size(badRowMaskMc), [H, numFrames])
        warning('VoltImg:laserBadRowsSize', ...
            'badRowMask size mismatch trial %d; ignoring mask for this trial.', tt);
        badRowMaskMc = false(H, numFrames);
    end
end

meanImgThisTrial = mean(single(imageStack), 3);
startTimeImaging = floor(ctx.startTime * ctx.imagingFreq);
meanImgThisTrialDouble = im2double(meanImgThisTrial);

seqThis = ctx.voltMapping.outParams.sequenceThisTrial{tt};
if isempty(seqThis)
    seqThis = ctx.zeroDummySequence;
end
holoSeqThisTrial = (unique(seqThis, 'stable') - min(unique(seqThis, 'stable')) + 1)';

fstVecTrial = ctx.voltMapping.outParams.firstStimTimes{cc};
if isempty(fstVecTrial)
    fstVecTrial = ctx.voltMapping.outParams.firstStimTimes{1, 2};
end

fineRoiXAllCells = cell(ctx.nCells, 1);
fineRoiYAllCells = cell(ctx.nCells, 1);
for nn = 1:ctx.nCells
    roiMeanThisTrial = zeros(size(meanImgThisTrial));
    for rr = 1:length(ctx.roughRoiXAllCells{nn})
        roiMeanThisTrial(ctx.roughRoiXAllCells{nn}(rr), ctx.roughRoiYAllCells{nn}(rr)) = ...
            meanImgThisTrial(ctx.roughRoiXAllCells{nn}(rr), ctx.roughRoiYAllCells{nn}(rr));
    end

    roiMaxStackDouble = im2double(meanImgThisTrial);
    roiPixels = roiMaxStackDouble > 0;
    roiMaxStackFilt = imgaussfilt(roiMaxStackDouble, 0.7);
    vals = roiMaxStackFilt(roiPixels);
    if ~isempty(vals)
        lo = prctile(vals, 10);
        hi = prctile(vals, 99);
        roiMaxStackNorm = (roiMaxStackFilt - lo) / max(hi - lo, eps);
        roiMaxStackNorm = min(max(roiMaxStackNorm, 0), 1);
    else
        roiMaxStackNorm = roiMaxStackFilt;
    end

    roiMaxStackRidge = zeros(size(roiMaxStackNorm));
    roiMask = false(size(roiMaxStackNorm));
    for rr = 1:length(ctx.roughRoiXAllCells{nn})
        roiMask(ctx.roughRoiXAllCells{nn}(rr), ctx.roughRoiYAllCells{nn}(rr)) = true;
    end
    [rows, cols] = find(roiMask);
    rmin = max(min(rows) - 2, 1);
    rmax = min(max(rows) + 2, size(roiMask, 1));
    cmin = max(min(cols) - 2, 1);
    cmax = min(max(cols) + 2, size(roiMask, 2));
    subImg = roiMaxStackNorm(rmin:rmax, cmin:cmax);
    subMask = roiMask(rmin:rmax, cmin:cmax);
    subRidge = fibermetric(subImg, 'StructureSensitivity', 2);
    subRidge(~subMask) = 0;
    roiMaxStackRidge(rmin:rmax, cmin:cmax) = subRidge;

    valsR = nonzeros(roiMaxStackRidge);
    if ~isempty(valsR)
        thr = prctile(valsR, 60);
    else
        thr = 0;
    end
    roiMaxStackRidgeReduced = roiMaxStackRidge;
    roiMaxStackRidgeReduced(roiMaxStackRidgeReduced < thr) = 0;
    roiMaxStackRidgeReduced(roiMaxStackRidgeReduced > 0) = 1;

    [roiFineX, roiFineY] = find(roiMaxStackRidgeReduced);
    if isempty(roiFineX)
        roiFineX = ctx.roughRoiXAllCells{nn};
        roiFineY = ctx.roughRoiYAllCells{nn};
    end

    fineRoiXAllCells{nn} = roiFineX;
    fineRoiYAllCells{nn} = roiFineY;
end

pack.fineRoiXAllCells = fineRoiXAllCells;
pack.fineRoiYAllCells = fineRoiYAllCells;

allTrialRoiMask = false(size(meanImgThisTrial));
for nn = 1:ctx.nCells
    if ~isempty(fineRoiXAllCells{nn})
        trialInd = sub2ind(size(allTrialRoiMask), fineRoiXAllCells{nn}, fineRoiYAllCells{nn});
        allTrialRoiMask(trialInd) = true;
    end
end

roiMeanFcell = cell(ctx.nCells, 1);
bkgrndMeanFcell = cell(ctx.nCells, 1);
roiMeanFCorrectedCell = cell(ctx.nCells, 1);
alphaScalarVec = nan(ctx.nCells, 1);
bkgrndRoiXAllCells_trial = cell(ctx.nCells, 1);
bkgrndRoiYAllCells_trial = cell(ctx.nCells, 1);

f0Cell = cell(ctx.nCells, 1);
dFFCell = cell(ctx.nCells, 1);
filtDFFCell = cell(ctx.nCells, 1);

for nn = 1:ctx.nCells
    if ctx.laserArtifactMcSecondSweepForDff
        roiMeanF = VoltImg_roiMeanPerFrameExcludeBadRows(imageStack, ...
            fineRoiXAllCells{nn}, fineRoiYAllCells{nn}, badRowMaskMc);
    else
        rawWholeRoiF = imageStack(fineRoiXAllCells{nn}, fineRoiYAllCells{nn}, :);
        roiMeanF = zeros(size(rawWholeRoiF, 3), 1);
        for ff = 1:size(rawWholeRoiF, 3)
            roiMeanF(ff, 1) = mean(mean(rawWholeRoiF(:, :, ff)));
        end
    end

    innerBuffer = 2;
    ringWidth = 3;
    minArea = 50;

    roiMaskThisCell = false(size(meanImgThisTrial));
    roiIndThisCell = sub2ind(size(roiMaskThisCell), fineRoiXAllCells{nn}, fineRoiYAllCells{nn});
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
        if ~isempty(ctx.bkgrndRoiXAllCells{nn})
            globalInd = sub2ind(size(ringGlobalMask), ctx.bkgrndRoiXAllCells{nn}, ctx.bkgrndRoiYAllCells{nn});
            ringGlobalMask(globalInd) = true;
        end
        ringClean = ringGlobalMask & ~allTrialRoiMask;
    end

    [bkgrndRoiXTrial, bkgrndRoiYTrial] = find(ringClean);
    bkgrndRoiXAllCells_trial{nn} = bkgrndRoiXTrial;
    bkgrndRoiYAllCells_trial{nn} = bkgrndRoiYTrial;

    if ctx.laserArtifactMcSecondSweepForDff
        bkgrndMeanF = VoltImg_roiMeanPerFrameExcludeBadRows(imageStack, ...
            bkgrndRoiXTrial, bkgrndRoiYTrial, badRowMaskMc);
    else
        rawWholeBkgrndF = imageStack(bkgrndRoiXTrial, bkgrndRoiYTrial, :);
        bkgrndMeanF = zeros(size(rawWholeBkgrndF, 3), 1);
        for ff = 1:size(rawWholeBkgrndF, 3)
            bkgrndMeanF(ff, 1) = mean(mean(rawWholeBkgrndF(:, :, ff)));
        end
    end

    baselineIndices = 1:startTimeImaging;
    bFit = robustfit(bkgrndMeanF(baselineIndices), roiMeanF(baselineIndices));
    alphaScalar = bFit(2);
    alphaScalar = min(max(alphaScalar, 0), 1);
    if alphaScalar > 0.8
        alphaScalar = 0.8;
    end
    alphaScalar = 0.80;

    roiMeanFCorrected = roiMeanF - alphaScalar * bkgrndMeanF;

    cutOffFreqIm = 40;
    [bIm, aIm] = butter(4, cutOffFreqIm / (ctx.imagingFreq / 2));
    if ~pack.isExcluded
        roiMeanFCorrectedFilt = filter(bIm, aIm, roiMeanFCorrected(:));
    else
        roiMeanFCorrectedFilt = [];
    end

    f0H = nan(1, nHolosThis);
    dFFH = cell(1, nHolosThis);
    filtDFFH = cell(1, nHolosThis);

    for hh = 1:nHolosThis
        iHoloLo = floor((fstVecTrial(hh) - ctx.preStimWindow / 1000) * ctx.imagingFreq);
        iHoloHi = ceil((fstVecTrial(hh) - ctx.preStimWindow / 1000) * ctx.imagingFreq) + ...
            ceil((ctx.ipi * ctx.nPulses + (ctx.preStimWindow + ctx.postStimWindow)) / 1000 * ctx.imagingFreq);

        if ~pack.isExcluded
            roiFCorrectedThisHolo = roiMeanFCorrected(iHoloLo:iHoloHi);
            roiFCorrectedThisHoloPreStim = roiFCorrectedThisHolo(1:(ctx.preStimWindow / 1000 * ctx.imagingFreq) - 1);
            f0ThisHolo = mean(roiFCorrectedThisHoloPreStim);
            dFThisHolo = roiFCorrectedThisHolo - f0ThisHolo;
            dFFThisHolo = dFThisHolo / f0ThisHolo;
            if ctx.UpOrDown == '2'
                dFFThisHolo = -dFFThisHolo;
            end

            roiFCorrectedThisHoloFilt = roiMeanFCorrectedFilt(iHoloLo:iHoloHi);
            roiFCorrectedThisHoloPreStimFilt = roiFCorrectedThisHoloFilt(1:(ctx.preStimWindow / 1000 * ctx.imagingFreq) - 1);
            f0FiltThisHolo = mean(roiFCorrectedThisHoloPreStimFilt);
            dFFiltThisHolo = (roiFCorrectedThisHoloFilt - f0FiltThisHolo) / f0FiltThisHolo;
            if ctx.UpOrDown == '2'
                filtdffThisHolo = -dFFiltThisHolo;
            else
                filtdffThisHolo = dFFiltThisHolo;
            end
        else
            f0ThisHolo = NaN;
            dFFThisHolo = NaN(ceil((ctx.ipi * ctx.nPulses + (ctx.preStimWindow + ctx.postStimWindow)) / 1000 * ctx.imagingFreq) + 2, 1);
            filtdffThisHolo = dFFThisHolo;
        end

        f0H(hh) = f0ThisHolo;
        dFFH{hh} = dFFThisHolo;
        filtDFFH{hh} = filtdffThisHolo;
    end

    roiMeanFcell{nn} = roiMeanF;
    bkgrndMeanFcell{nn} = bkgrndMeanF;
    roiMeanFCorrectedCell{nn} = roiMeanFCorrected;
    alphaScalarVec(nn) = alphaScalar;
    f0Cell{nn} = f0H;
    dFFCell{nn} = dFFH;
    filtDFFCell{nn} = filtDFFH;
end

pack.bkgrndRoiXAllCells_trial = bkgrndRoiXAllCells_trial;
pack.bkgrndRoiYAllCells_trial = bkgrndRoiYAllCells_trial;
pack.roiMeanFcell = roiMeanFcell;
pack.bkgrndMeanFcell = bkgrndMeanFcell;
pack.roiMeanFCorrectedCell = roiMeanFCorrectedCell;
pack.alphaScalarVec = alphaScalarVec;
pack.holoSeqThisTrial = holoSeqThisTrial;
pack.nHolosThis = nHolosThis;
pack.f0Cell = f0Cell;
pack.dFFCell = dFFCell;
pack.filtDFFCell = filtDFFCell;

end
