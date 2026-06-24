%% Voltage Imaging Slice Test Analysis — BATCH / headless variant (042026_batch)
% Same pipeline as VoltImg_slice_test_analysis_MCfineROI_TrialSpecROIandNeuropil_042026.m
% but with fixed paths (no uigetdir), no blocking input(), and automatic rough ROI
% (disk around brightest pixel on motion-corrected mean) instead of drawfreehand.
% Edit the CONFIG block below for each experiment.
%
% New version (042026) based on VoltImg_slice_test_analysis_100124.m with:
% 1) Per-trial NoRMCorre rigid motion correction anchored to a global mean
%    template (same strategy as VoltImg_mapping_analysis_MultiCell_newDFF_
%    021226_MCfineROI_TrialSpecROIandNeuropil.m). Motion-corrected stacks are
%    saved as multi-page TIFFs (*_mc.tif).
% 2) One patched neuron: single hand-drawn rough ROI on motion-corrected data,
%    then per-trial fine ROI (gaussian + fibermetric inside rough ROI) and
%    trial-specific neuropil ring (excluding ROI pixels).
% 3) Neuropil subtraction and dF/F aligned with the mapping pipeline:
%    robustfit alpha on prestim frames (clamped), fixed alpha override 0.90,
%    F0 from prestim on neuropil-corrected trace, dF/F = (Fcorr - F0)/F0,
%    Butterworth filter on dF/F with padded filter() startup (40 Hz cutoff).
% 4) No holography: full-trial traces are kept; stimulation is current
%    injection; dvCondSequence still assigns trials to amplitude conditions.
%
% Requirements: NoRMCorre on MATLAB path (set normcorrePath below).
% Image Processing Toolbox: fibermetric, imgaussfilt, strel, imdilate, etc.

%% Initialize
clear all
% close all
set(0, 'DefaultFigureVisible', 'off');

%% CONFIG — fixed paths for non-interactive runs (edit per dataset)
ephysFilePath = '/Volumes/phoenixinthesky/Masato/Voltage Imaging Data_Phoenix/voltImgTest/DAQ_Ephys Data/250618/vGlutCre_DIOASAP7y_IC_Slice_1020nm_061825_cell1_dvTest_AH';
ImgsFilePath = '/Volumes/phoenixinthesky/Masato/Voltage Imaging Data_Phoenix/voltImgTest/Raw Imaging Data/061825_AH/vGlutCre_DIOASAP7y_IC_Slice_1020nm_061825_cell1_dvTest';
savePathRoot = '/Volumes/phoenixinthesky/Masato/Voltage Imaging Data_Phoenix/voltImgTest/MC Imaging Data';
% Disk radius (pixels) for automatic rough ROI around brightest pixel on mean MC image
batchRoughRoiRadiusPx = 25;
% GEVI polarity: '1' upward (ASAP6, FORCE), '2' downward (Jedi2P, ASAP5, ASAP7y, …)
UpOrDown = '2';

%% Load files and setup
% Step 1: Load ephys (.mat in folder; first file when sorted by name if multiple)
matFiles = dir(fullfile(ephysFilePath, '*.mat'));
if isempty(matFiles)
    error('No .mat files found in ephys folder: %s', ephysFilePath);
end
[~, sortIx] = sort({matFiles.name});
matFiles = matFiles(sortIx);
load(fullfile(matFiles(1).folder, matFiles(1).name));

% NoRMCorre (same as interactive script: empty = rely on MATLAB path)
normcorrePath = ''; % e.g. 'C:\path\to\NoRMCorre-master'
if ~isempty(normcorrePath) && exist(normcorrePath, 'dir')
    addpath(normcorrePath);
end

% Get dv parameters and condition sequence
dvCondSequence = ExpStruct.dvStepParams.dvCondSequence;
dvToTest = ExpStruct.dvStepParams.dvToTest;
nConds = length(unique(dvCondSequence));

pulseStart = ExpStruct.dvStepParams.pulseStart;
sweepDur = ExpStruct.dvStepParams.sweepDur;
nPulses = ExpStruct.dvStepParams.nPulses;
pulseFreq = ExpStruct.dvStepParams.pulseFreq;
imagingFreq = 330.22;
Fs = ExpStruct.Fs;

% Step 2: Imaging folder (fixed path from CONFIG)
ImgfolderContents = dir(ImgsFilePath);

% Step 2a: Avoid hidden files and non-image files
fileNames = [];
fileType = '.tif';
for ii = 1:length(ImgfolderContents)
    if ~ImgfolderContents(ii).isdir && ~startsWith(ImgfolderContents(ii).name, '.') && endsWith(ImgfolderContents(ii).name, fileType)
        fileNames{ii, 1} = ImgfolderContents(ii).name;
    end
end
imagesIndex = find(~cellfun(@isempty, fileNames));

% Step 2b: Revise cond sequence to match actual trials recorded
dvCondSequence = dvCondSequence(1:numel(imagesIndex));

% Step 3: Bad trials (baseline Vm)
if isfield(ExpStruct.dvStepParams, 'vsTest_inputs')
    vsTest_inputs = ExpStruct.dvStepParams.vsTest_inputs;
else
    vsTest_inputs = ExpStruct.dvStepParams.vsTest_inputs_Ch1;
end

baselineAllTrials = [];
excludeTrials = [];
for tt = 1:size(vsTest_inputs, 2)
    baseline = mean(vsTest_inputs(1:round(0.001 * pulseStart * Fs), tt));
    baselineAllTrials = [baselineAllTrials, baseline]; %#ok<AGROW>
    if baselineAllTrials(tt) > -55
        excludeTrials = [excludeTrials, tt]; %#ok<AGROW>
    end
end

% Step 4: GEVI orientation — set in CONFIG (UpOrDown)

% Step 5: Frame clock offsets (unchanged from original)
frameClock_inputs = [];
risingEdgeIndices = [];
risingEdgeTimes = [];
if isfield(ExpStruct.dvStepParams, 'frameClock_inputs')
    frameClock_inputs = ExpStruct.dvStepParams.frameClock_inputs;
    for tt = 1:size(frameClock_inputs, 2)
        trialTime = (0:size(frameClock_inputs, 1) - 1) * (1 / Fs);
        riseThreshold = 0;
        if any(frameClock_inputs(:, tt) == 1)
            risingEdgeIndices(tt) = find(frameClock_inputs(:, tt) > riseThreshold, 1);
            risingEdgeTimes(tt) = trialTime(risingEdgeIndices(tt));
        else
            risingEdgeIndices(tt) = 0;
            risingEdgeTimes(tt) = 0;
        end
    end
else
    risingEdgeTimes = repmat(0.0011, 1, size(vsTest_inputs, 2) - numel(excludeTrials));
end

% Auto-detect single vs interleaved two-color TIFF (same logic as mapping script)
if isempty(imagesIndex)
    error('No TIFF files found in selected imaging folder.');
end
testTiffPath = fullfile(ImgfolderContents(imagesIndex(1)).folder, ImgfolderContents(imagesIndex(1)).name);
tTest = Tiff(testTiffPath, 'r');
nDirs = 1;
while true
    try
        tTest.setDirectory(nDirs);
        nDirs = nDirs + 1;
    catch
        nDirs = nDirs - 1;
        break
    end
end
if nDirs < 4
    rawImgNChannels = 1;
    tTest.close();
else
    nSample = min(120, nDirs);
    frameMeans = zeros(nSample, 1);
    oddMeanImg = [];
    evenMeanImg = [];
    nOdd = 0;
    nEven = 0;
    for pp = 1:nSample
        tTest.setDirectory(pp);
        frame = single(tTest.read());
        frameMeans(pp) = mean(frame(:));
        if mod(pp, 2) == 1
            if isempty(oddMeanImg)
                oddMeanImg = zeros(size(frame), 'single');
            end
            oddMeanImg = oddMeanImg + frame;
            nOdd = nOdd + 1;
        else
            if isempty(evenMeanImg)
                evenMeanImg = zeros(size(frame), 'single');
            end
            evenMeanImg = evenMeanImg + frame;
            nEven = nEven + 1;
        end
    end
    tTest.close();
    oddMeanImg = oddMeanImg ./ max(1, nOdd);
    evenMeanImg = evenMeanImg ./ max(1, nEven);
    rImg = corrcoef(double(oddMeanImg(:)), double(evenMeanImg(:)));
    if numel(rImg) >= 4
        oddEvenImgCorr = rImg(1, 2);
    else
        oddEvenImgCorr = 1;
    end
    if numel(frameMeans) >= 3
        lag1 = corrcoef(frameMeans(1:end - 1), frameMeans(2:end));
        lag2 = corrcoef(frameMeans(1:end - 2), frameMeans(3:end));
        if numel(lag1) >= 4
            lag1Corr = lag1(1, 2);
        else
            lag1Corr = 0;
        end
        if numel(lag2) >= 4
            lag2Corr = lag2(1, 2);
        else
            lag2Corr = 0;
        end
        altStepDiff = mean(abs(diff(frameMeans)));
        sameChanDiff = mean(abs(frameMeans(3:end) - frameMeans(1:end - 2)));
    else
        lag1Corr = 0;
        lag2Corr = 0;
        altStepDiff = 0;
        sameChanDiff = 0;
    end
    isInterleaved = (oddEvenImgCorr < 0.90) && ...
        ((lag2Corr > lag1Corr + 0.10) || (altStepDiff > 1.15 * max(sameChanDiff, eps)));
    if isInterleaved
        rawImgNChannels = 2;
    else
        rawImgNChannels = 1;
    end
    disp(['Auto rawImgNChannels = ', num2str(rawImgNChannels), ...
        ' (odd-even mean image corr = ', num2str(oddEvenImgCorr, '%.3f'), ').']);
end

% Output root for motion-corrected TIFFs (fixed path from CONFIG)
expID = num2str(ExpStruct.mouseID);
mouseTag = ['voltImgSliceTest_Analysis_', expID, '_MCfineROI_batch'];
saveDirectory = fullfile(savePathRoot, mouseTag);
if ~exist(saveDirectory, 'dir')
    mkdir(saveDirectory);
end
mcTiffFolder = fullfile(saveDirectory, 'Motion_Corrected_Tiffs');
if ~exist(mcTiffFolder, 'dir')
    mkdir(mcTiffFolder);
end

voltImgTest_Analysis.ephysData.Fs = Fs;
voltImgTest_Analysis.pulseParams.pulseStart = pulseStart;
voltImgTest_Analysis.pulseParams.sweepDur = sweepDur;
voltImgTest_Analysis.pulseParams.nPulses = nPulses;
voltImgTest_Analysis.pulseParams.pulseFreq = pulseFreq;
voltImgTest_Analysis.imagingFreq = imagingFreq;
voltImgTest_Analysis.dvCondSequence = dvCondSequence;
voltImgTest_Analysis.dvToTest = dvToTest;
voltImgTest_Analysis.Rinput = ExpStruct.dvStepParams.Rinput;
voltImgTest_Analysis.fRinput = ExpStruct.dvStepParams.fRinput;
voltImgTest_Analysis.ephysData.vsTest_inputs = vsTest_inputs;
voltImgTest_Analysis.ephysData.baselineAllTrials = baselineAllTrials;
voltImgTest_Analysis.ephysData.excludeTrials = excludeTrials;
voltImgTest_Analysis.rawImgNChannels = rawImgNChannels;
voltImgTest_Analysis.mcTiffFolder = mcTiffFolder;
voltImgTest_Analysis.saveDirectory = saveDirectory;
voltImgTest_Analysis.batchConfig.ephysFilePath = ephysFilePath;
voltImgTest_Analysis.batchConfig.ImgsFilePath = ImgsFilePath;
voltImgTest_Analysis.batchConfig.savePathRoot = savePathRoot;
voltImgTest_Analysis.batchConfig.batchRoughRoiRadiusPx = batchRoughRoiRadiusPx;
voltImgTest_Analysis.batchConfig.UpOrDown = UpOrDown;

%% Motion correction: global template + per-trial NoRMCorre (anchored)
if ~exist('NoRMCorreSetParms', 'file') || ~exist('normcorre', 'file')
    error(['NoRMCorre not found. Add NoRMCorre to the MATLAB path (set normcorrePath ', ...
        'near the top of this script) and re-run.']);
end

firstImgPath = fullfile(ImgfolderContents(imagesIndex(1)).folder, ImgfolderContents(imagesIndex(1)).name);
infoFirst = imfinfo(firstImgPath);
globalTemplateAccum = zeros(infoFirst(1).Height, infoFirst(1).Width, 'single');
nTemplateTrials = 0;

for tt = 1:numel(imagesIndex)
    disp(['Global template: trial ', num2str(tt), ' / ', num2str(numel(imagesIndex))]);
    if ismember(tt, excludeTrials)
        continue
    end
    currImgPath = fullfile(ImgfolderContents(imagesIndex(tt)).folder, ImgfolderContents(imagesIndex(tt)).name);
    [imageStack, ~] = slice_loadRawTiffStack(currImgPath, rawImgNChannels);
    globalTemplateAccum = globalTemplateAccum + mean(single(imageStack), 3);
    nTemplateTrials = nTemplateTrials + 1;
end
if nTemplateTrials == 0
    error('No non-excluded trials available to build a global motion-correction template.');
end
globalTemplate = globalTemplateAccum ./ nTemplateTrials;

for tt = 1:numel(imagesIndex)
    disp(['Motion correcting trial ', num2str(tt), ' of ', num2str(numel(imagesIndex))]);
    rawName = ImgfolderContents(imagesIndex(tt)).name;
    [~, baseName, ~] = fileparts(rawName);
    mcName = [baseName, '_mc.tif'];
    mcPath = fullfile(mcTiffFolder, mcName);

    currImgPath = fullfile(ImgfolderContents(imagesIndex(tt)).folder, ImgfolderContents(imagesIndex(tt)).name);
    [imageStack, d1] = slice_loadRawTiffStack(currImgPath, rawImgNChannels);
    [~, d2, ~] = size(imageStack);
    Y = single(imageStack);
    options_rigid = NoRMCorreSetParms('d1', d1, 'd2', d2, 'bin_width', 20, 'max_shift', 4, 'us_fac', 50, 'init_batch', 1);
    [M_mc, ~] = normcorre(Y, options_rigid, globalTemplate);
    imageStack_mc = M_mc;

    if ~isa(imageStack_mc, 'uint16')
        mcMin = min(imageStack_mc(:));
        mcMax = max(imageStack_mc(:));
        if mcMax > mcMin
            imageStack_mc_uint16 = uint16((imageStack_mc - mcMin) ./ (mcMax - mcMin) * double(intmax('uint16')));
        else
            imageStack_mc_uint16 = uint16(zeros(size(imageStack_mc)));
        end
    else
        imageStack_mc_uint16 = imageStack_mc;
    end

    tOut = Tiff(mcPath, 'w');
    tagstruct.ImageLength = d1;
    tagstruct.ImageWidth = d2;
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 16;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Compression = Tiff.Compression.LZW;
    for k = 1:size(imageStack_mc_uint16, 3)
        tOut.setTag(tagstruct);
        tOut.write(imageStack_mc_uint16(:, :, k));
        if k < size(imageStack_mc_uint16, 3)
            tOut.writeDirectory();
        end
    end
    tOut.close();
end

voltImgTest_Analysis.globalMcTemplate = globalTemplate;

%% Mean stack for ROI drawing (motion-corrected), same logic as original max-dv subset
maxDvTrials = [];
if UpOrDown == '1'
    maxDvTrials = find(dvCondSequence == max(unique(dvCondSequence)));
elseif UpOrDown == '2'
    maxDvTrials = find(dvCondSequence == min(unique(dvCondSequence)));
end
if numel(maxDvTrials) > 100
    maxDvTrials = randsample(maxDvTrials, 100);
end

maxDvStackPieces = {};
counter = 0;
for tt = 1:numel(maxDvTrials)
    counter = counter + 1;
    disp(['Mean stack for ROI: ', num2str(counter), ' / ', num2str(numel(maxDvTrials))]);
    trialIdx = maxDvTrials(tt);
    if ismember(trialIdx, excludeTrials)
        continue
    end
    rawName = ImgfolderContents(imagesIndex(trialIdx)).name;
    [~, baseName, ~] = fileparts(rawName);
    mcPath = fullfile(mcTiffFolder, [baseName, '_mc.tif']);
    [imageStack_mc, ~] = slice_loadMcTiffStack(mcPath);
    imageStackMean = mean(single(imageStack_mc), 3);
    maxDvStackPieces{end + 1} = imageStackMean; %#ok<SAGROW>
end
if isempty(maxDvStackPieces)
    error('No max-dv trials available after exclusions to build mean image for ROI.');
end
maxDvStack = cat(3, maxDvStackPieces{:});
meanMaxDvStack = mean(maxDvStack, 3);
meanFluorMaxDvStack = meanMaxDvStack;

%%%%% Rough ROI: automatic disk around brightest pixel (batch; replaces drawfreehand)
if ~ismatrix(meanFluorMaxDvStack)
    error('meanFluorMaxDvStack must be 2-D.');
end
[Hmean, Wmean] = size(meanFluorMaxDvStack);
[~, linMax] = max(meanFluorMaxDvStack(:));
[r0, c0] = ind2sub([Hmean, Wmean], linMax);
[Ir, Ic] = ndgrid(1:Hmean, 1:Wmean);
distPx = hypot(double(Ir) - double(r0), double(Ic) - double(c0));
autoRoughMask = distPx <= double(batchRoughRoiRadiusPx);
[roughRoiX, roughRoiY] = find(autoRoughMask);
if isempty(roughRoiX)
    error('batchRoughRoiRadiusPx produced empty rough ROI; increase radius.');
end

roughRoiXAllCells = {roughRoiX};
roughRoiYAllCells = {roughRoiY};

%% Global fine ROI + fallback neuropil ring on mean image (for visualization / fallback)
roiMeanMaxDvStack = zeros(size(maxDvStack, 1), size(maxDvStack, 2));
for rr = 1:numel(roughRoiX)
    roiMeanMaxDvStack(roughRoiX(rr), roughRoiY(rr)) = mean(maxDvStack(roughRoiX(rr), roughRoiY(rr), :), 3);
end
figure(11); clf;
imagesc(roiMeanMaxDvStack); axis equal; axis image; axis off; colorbar off;

roiMaxStackDouble = im2double(meanFluorMaxDvStack);
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
roiMaxStackRidge = fibermetric(roiMaxStackNorm, 'StructureSensitivity', 2);
valsR = nonzeros(roiMaxStackRidge);
if ~isempty(valsR)
    thr = prctile(valsR, 50);
else
    thr = 0;
end
roiMaxStackRidgeReduced = roiMaxStackRidge;
roiMaxStackRidgeReduced(roiMaxStackRidgeReduced < thr) = 0;
roiMaxStackRidgeReduced(roiMaxStackRidgeReduced > 0) = 1;
[roiX_global, roiY_global] = find(roiMaxStackRidgeReduced);
if isempty(roiX_global)
    roiX_global = roughRoiX;
    roiY_global = roughRoiY;
end

innerBuffer = 2;
ringWidth = 3;
minArea = 50;
innerSelect = imdilate(roiMaxStackRidge, strel('disk', innerBuffer));
outerSelect = imdilate(roiMaxStackRidge, strel('disk', innerBuffer + ringWidth));
backgroundRing = outerSelect & ~innerSelect;
valsBk = roiMaxStackDouble(backgroundRing);
if ~isempty(valsBk)
    brightCut = prctile(valsBk, 95);
    ringClean = backgroundRing & (roiMaxStackDouble <= brightCut);
else
    ringClean = backgroundRing;
end
ringClean = bwareaopen(ringClean, 7);
if nnz(ringClean) < minArea
    ringClean = backgroundRing;
end
[bkgrndRoiX_global, bkgrndRoiY_global] = find(ringClean);

voltImgTest_Analysis.maxDvTrials = maxDvTrials;
voltImgTest_Analysis.meanFluorMaxDvStack = meanFluorMaxDvStack;
voltImgTest_Analysis.roughRoiXAllCells = roughRoiXAllCells;
voltImgTest_Analysis.roughRoiYAllCells = roughRoiYAllCells;
voltImgTest_Analysis.roughRoiX = roughRoiX;
voltImgTest_Analysis.roughRoiY = roughRoiY;
voltImgTest_Analysis.roiX_global = roiX_global;
voltImgTest_Analysis.roiY_global = roiY_global;
voltImgTest_Analysis.bkgrndRoiX_global = bkgrndRoiX_global;
voltImgTest_Analysis.bkgrndRoiY_global = bkgrndRoiY_global;
% Global fine ROI mask (same role as roiStack / roiX / roiY in the original slice script)
voltImgTest_Analysis.roiStack = roiMaxStackRidgeReduced;
voltImgTest_Analysis.roiX = roiX_global;
voltImgTest_Analysis.roiY = roiY_global;

%% Per-trial fine ROI, neuropil, corrected F, dF/F (mapping-style)
startTimeImaging = max(2, floor((pulseStart / 1000) * imagingFreq));

fineRoiXAll = cell(numel(imagesIndex), 1);
fineRoiYAll = cell(numel(imagesIndex), 1);
bkgrndRoiX_trial = cell(numel(imagesIndex), 1);
bkgrndRoiY_trial = cell(numel(imagesIndex), 1);
alphaScalarAll = nan(numel(imagesIndex), 1);

cutOffFreqIm = 40;
[bIm, aIm] = butter(4, cutOffFreqIm / (imagingFreq / 2));
padLenFilt = max(60, 6 * max(numel(bIm), numel(aIm)));

dfAllConds = cell(nConds, 1);
dffAllConds = cell(nConds, 1);
filtdffAllConds = cell(nConds, 1);
traceAllConds = cell(nConds, 1);

for tt = 1:size(vsTest_inputs, 2)
    disp(['Trace extraction trial ', num2str(tt), ' / ', num2str(size(vsTest_inputs, 2))]);
    rawName = ImgfolderContents(imagesIndex(tt)).name;
    [~, baseName, ~] = fileparts(rawName);
    mcPath = fullfile(mcTiffFolder, [baseName, '_mc.tif']);
    [imageStack, numFrames] = slice_loadMcTiffStack(mcPath);

    meanImgThisTrial = mean(single(imageStack), 3);
    meanImgThisTrialDouble = im2double(meanImgThisTrial);

    % Fine ROI (trial-specific)
    roiMask = false(size(meanImgThisTrial));
    for rr = 1:numel(roughRoiX)
        roiMask(roughRoiX(rr), roughRoiY(rr)) = true;
    end
    [rows, cols] = find(roiMask);
    rmin = max(min(rows) - 2, 1); rmax = min(max(rows) + 2, size(roiMask, 1));
    cmin = max(min(cols) - 2, 1); cmax = min(max(cols) + 2, size(roiMask, 2));

    roiMaxStackDouble = im2double(meanImgThisTrial);
    roiPixelsTrial = roiMaxStackDouble > 0;
    roiMaxStackFilt = imgaussfilt(roiMaxStackDouble, 0.7);
    valsTrial = roiMaxStackFilt(roiPixelsTrial);
    if ~isempty(valsTrial)
        lo = prctile(valsTrial, 10);
        hi = prctile(valsTrial, 99);
        roiMaxStackNorm = (roiMaxStackFilt - lo) / max(hi - lo, eps);
        roiMaxStackNorm = min(max(roiMaxStackNorm, 0), 1);
    else
        roiMaxStackNorm = roiMaxStackFilt;
    end
    roiMaxStackRidge = zeros(size(roiMaxStackNorm));
    subImg = roiMaxStackNorm(rmin:rmax, cmin:cmax);
    subMask = roiMask(rmin:rmax, cmin:cmax);
    subRidge = fibermetric(subImg, 'StructureSensitivity', 2);
    subRidge(~subMask) = 0;
    roiMaxStackRidge(rmin:rmax, cmin:cmax) = subRidge;
    valsR2 = nonzeros(roiMaxStackRidge);
    if ~isempty(valsR2)
        thr2 = prctile(valsR2, 60);
    else
        thr2 = 0;
    end
    roiMaxStackRidgeReduced = roiMaxStackRidge;
    roiMaxStackRidgeReduced(roiMaxStackRidgeReduced < thr2) = 0;
    roiMaxStackRidgeReduced(roiMaxStackRidgeReduced > 0) = 1;
    [roiFineX, roiFineY] = find(roiMaxStackRidgeReduced);
    if isempty(roiFineX)
        roiFineX = roughRoiX;
        roiFineY = roughRoiY;
    end
    fineRoiXAll{tt} = roiFineX;
    fineRoiYAll{tt} = roiFineY;

    rawWholeRoiF = imageStack(roiFineX, roiFineY, :);
    roiMeanF = zeros(numFrames, 1);
    for ff = 1:numFrames
        roiMeanF(ff) = mean(mean(rawWholeRoiF(:, :, ff)));
    end

    allTrialRoiMask = false(size(meanImgThisTrial));
    allTrialRoiMask(sub2ind(size(allTrialRoiMask), roiFineX, roiFineY)) = true;

    innerBuffer = 2;
    ringWidth = 3;
    minArea = 50;
    roiMaskThisCell = false(size(meanImgThisTrial));
    roiMaskThisCell(sub2ind(size(roiMaskThisCell), roiFineX, roiFineY)) = true;
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
        if ~isempty(bkgrndRoiX_global)
            ringGlobalMask(sub2ind(size(ringGlobalMask), bkgrndRoiX_global, bkgrndRoiY_global)) = true;
        end
        ringClean = ringGlobalMask & ~allTrialRoiMask;
    end
    [bkgrndRoiXTrial, bkgrndRoiYTrial] = find(ringClean);
    bkgrndRoiX_trial{tt} = bkgrndRoiXTrial;
    bkgrndRoiY_trial{tt} = bkgrndRoiYTrial;

    rawWholeBkgrndF = imageStack(bkgrndRoiXTrial, bkgrndRoiYTrial, :);
    bkgrndMeanF = zeros(numFrames, 1);
    for ff = 1:numFrames
        bkgrndMeanF(ff) = mean(mean(rawWholeBkgrndF(:, :, ff)));
    end

    biRobust = 1:min(startTimeImaging - 1, numFrames);
    if numel(biRobust) < 2
        biRobust = 1:min(2, numFrames);
    end
    bFit = robustfit(bkgrndMeanF(biRobust), roiMeanF(biRobust));
    alphaScalar = bFit(2);
    alphaScalar = min(max(alphaScalar, 0), 1);
    if alphaScalar > 0.8
        alphaScalar = 0.8;
    end
    alphaScalar = 0.90;
    alphaScalarAll(tt) = alphaScalar;

    roiMeanFCorrected = roiMeanF - alphaScalar * bkgrndMeanF;

    idxF0_start = max(1, floor((pulseStart - pulseStart / 2) / 1000 * imagingFreq));
    idxF0_end = min(numFrames, floor(pulseStart / 1000 * imagingFreq));
    if idxF0_start > idxF0_end
        idxF0_start = 1;
        idxF0_end = min(numFrames, max(1, startTimeImaging - 1));
    end
    f0 = mean(roiMeanFCorrected(idxF0_start:idxF0_end));

    df = roiMeanFCorrected - f0;
    dff = df / f0;

    if UpOrDown == '2'
        df = -df;
        dff = -dff;
    end

    dFFcol = dff(:);
    xPad = [repmat(dFFcol(1), padLenFilt, 1); dFFcol];
    yPad = filter(bIm, aIm, xPad);
    filtdff = yPad(padLenFilt + 1:end);

    condIdx = dvCondSequence(tt);
    dfAllConds{condIdx, 1}(:, end + 1) = df;
    dffAllConds{condIdx, 1}(:, end + 1) = dff;
    filtdffAllConds{condIdx, 1}(:, end + 1) = filtdff;
    traceAllConds{condIdx, 1}(:, end + 1) = vsTest_inputs(:, tt);
end

%% Trial means and CI (unchanged structure from original script)
meanDfAllConds = [];
meanDffAllConds = [];
meanFiltDffAllConds = [];
meanTraceAllConds = [];
for cc = 1:nConds
    meanDfAllConds = [meanDfAllConds, mean(dfAllConds{cc}, 2)]; %#ok<AGROW>
    meanDffAllConds = [meanDffAllConds, mean(dffAllConds{cc}, 2)]; %#ok<AGROW>
    meanFiltDffAllConds = [meanFiltDffAllConds, mean(filtdffAllConds{cc}, 2)]; %#ok<AGROW>
    meanTraceAllConds = [meanTraceAllConds, mean(traceAllConds{cc}, 2)]; %#ok<AGROW>
end

peakDff = [];
peakFiltDff = [];
nImRows = size(meanDffAllConds, 1);
for cc = 1:nConds
    for pp = 1:nPulses
        tMsStart = pulseStart + (1000 / pulseFreq) * (pp - 1);
        tMsEnd = pulseStart + (1000 / pulseFreq) * pp;
        % Frame indices from time in ms (full expression inside floor); clamp to valid rows
        idxStart = max(1, min(nImRows, floor((tMsStart / 1000) * imagingFreq)));
        idxEnd = max(idxStart, min(nImRows, floor((tMsEnd / 1000) * imagingFreq)));
        peakDff{cc}(:, pp) = max(meanDffAllConds(idxStart:idxEnd, cc));
        peakFiltDff{cc}(:, pp) = max(meanFiltDffAllConds(idxStart:idxEnd, cc));
    end
end

CIDfAllConds = [];
CIDffAllConds = [];
CIFiltDffAllConds = [];
for cc = 1:nConds
    confidence_level = 0.95;
    meansDff = mean(dffAllConds{cc}, 2);
    std_errorsDff = std(dffAllConds{cc}, 0, 2) / sqrt(size(dffAllConds{cc}, 2));
    t_scoreDff = tinv((1 + confidence_level) / 2, size(dffAllConds{cc}, 2) - 1);
    margin_of_errorDff = t_scoreDff .* std_errorsDff;
    lower_boundsDff = meansDff - margin_of_errorDff;
    upper_boundsDff = meansDff + margin_of_errorDff;

    meansDf = mean(dfAllConds{cc}, 2);
    std_errorsDf = std(dfAllConds{cc}, 0, 2) / sqrt(size(dfAllConds{cc}, 2));
    t_scoreDf = tinv((1 + confidence_level) / 2, size(dfAllConds{cc}, 2) - 1);
    margin_of_errorDf = t_scoreDf .* std_errorsDf;
    lower_boundsDf = meansDf - margin_of_errorDf;
    upper_boundsDf = meansDf + margin_of_errorDf;

    meansFiltDff = mean(filtdffAllConds{cc}, 2);
    std_errorsFiltDff = std(filtdffAllConds{cc}, 0, 2) / sqrt(size(filtdffAllConds{cc}, 2));
    t_scoreFiltDff = tinv((1 + confidence_level) / 2, size(filtdffAllConds{cc}, 2) - 1);
    margin_of_errorFiltDff = t_scoreFiltDff .* std_errorsFiltDff;
    lower_boundsFiltDff = meansFiltDff - margin_of_errorFiltDff;
    upper_boundsFiltDff = meansFiltDff + margin_of_errorFiltDff;

    CIDffAllConds{cc, 1} = [lower_boundsDff, upper_boundsDff];
    CIDfAllConds{cc, 1} = [lower_boundsDf, upper_boundsDf];
    CIFiltDffAllConds{cc, 1} = [lower_boundsFiltDff, upper_boundsFiltDff];
end

voltImgTest_Analysis.maxDvStack = maxDvStack;
voltImgTest_Analysis.fineRoiXAll = fineRoiXAll;
voltImgTest_Analysis.fineRoiYAll = fineRoiYAll;
voltImgTest_Analysis.bkgrndRoiX_trial = bkgrndRoiX_trial;
voltImgTest_Analysis.bkgrndRoiY_trial = bkgrndRoiY_trial;
voltImgTest_Analysis.alphaScalarAll = alphaScalarAll;
voltImgTest_Analysis.dfAllConds = dfAllConds;
voltImgTest_Analysis.dffAllConds = dffAllConds;
voltImgTest_Analysis.filtdffAllConds = filtdffAllConds;
voltImgTest_Analysis.traceAllConds = traceAllConds;
voltImgTest_Analysis.meanDffAllConds = meanDffAllConds;
voltImgTest_Analysis.meanDfAllConds = meanDfAllConds;
voltImgTest_Analysis.meanFiltDffAllConds = meanFiltDffAllConds;
voltImgTest_Analysis.meanTraceAllConds = meanTraceAllConds;
voltImgTest_Analysis.peakDff = peakDff;
voltImgTest_Analysis.peakFiltDff = peakFiltDff;
voltImgTest_Analysis.CIDffAllConds = CIDffAllConds;
voltImgTest_Analysis.CIDfAllConds = CIDfAllConds;
voltImgTest_Analysis.CIFiltDffAllConds = CIFiltDffAllConds;
voltImgTest_Analysis.mouseID = ['voltImgSliceTest_Analysis_', num2str(ExpStruct.mouseID), '_MCfineROI_batch'];

%% Align imaging with ephys traces
for cc = 1:nConds
    figure(30 + cc)
    clf
    hold on
    fill([linspace(0, size(CIDffAllConds{cc, 1}, 1) / imagingFreq, size(CIDffAllConds{cc, 1}, 1)), fliplr(linspace(0, size(CIDffAllConds{cc, 1}, 1) / imagingFreq, size(CIDffAllConds{cc, 1}, 1)))], ...
        [CIDffAllConds{cc, 1}(:, 1)', fliplr(CIDffAllConds{cc, 1}(:, 2)')], [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
    plot(linspace(0, size(CIDffAllConds{cc, 1}, 1) / imagingFreq, size(CIDffAllConds{cc, 1}, 1)), CIDffAllConds{cc, 1}(:, 1) * 100, '-', 'linewidth', 1, 'color', [0.9 0.9 0.9]);
    plot(linspace(0, size(CIDffAllConds{cc, 1}, 1) / imagingFreq, size(CIDffAllConds{cc, 1}, 1)), CIDffAllConds{cc, 1}(:, 2) * 100, '-', 'linewidth', 1, 'color', [0.9 0.9 0.9]);
    plot(linspace(0, size(meanFiltDffAllConds, 1) / imagingFreq, size(meanFiltDffAllConds, 1)), meanFiltDffAllConds(:, cc) * 100, '-', 'linewidth', 2, 'color', 'g');
    set(gca, 'fontsize', 16);
    nCi = size(CIDffAllConds{cc, 1}, 1);
    % ylim used fixed "time" windows via frame count; clamp to actual trace length (short stacks << 0.2 s in frames)
    iLo = max(1, min(nCi, ceil(0.2 * imagingFreq)));
    iMid = max(iLo, min(nCi, ceil(0.7 * imagingFreq)));
    iHi = max(iLo, min(nCi, ceil(0.8 * imagingFreq)));
    xlim([0.2, 0.9]);
    if nCi >= 1
        yMin = min(CIDffAllConds{cc, 1}(iLo:iMid, 1)) * 100;
        yMax = max(CIDffAllConds{cc, 1}(iLo:iHi, 2)) * 100;
        if yMin == yMax
            yMax = yMin + eps;
        end
        ylim([yMin, yMax]);
    end
    ylabel('dF/F (%)');
    title(['dV = ', num2str(abs(mean(meanTraceAllConds(1:(pulseStart / 1000) * Fs, cc))) + (max(meanTraceAllConds(((pulseStart / 1000) * Fs):((pulseStart / 1000) * Fs + ((1000 / pulseFreq) / 1000 * Fs)), cc)))), 'mV,', ' (1st pulse) ', ...
        num2str(abs(mean(meanTraceAllConds(1:(pulseStart / 1000) * Fs, cc))) + (max(meanTraceAllConds(((pulseStart / 1000) * Fs + ((1000 / pulseFreq) / 1000 * Fs)):end, cc)))), 'mV,', ' (2nd pulse)'], 'fontsize', 10);
    xlabel('Time (s)');
    yline(0, '-', 'LineWidth', 1.5, 'color', [0.9 0.9 0.9]);
    hold off

    yyaxis right
    nEp = size(meanTraceAllConds, 1);
    epRefLo = max(1, min(nEp, round((250 / 1000) * Fs - 1000)));
    epRefHi = max(epRefLo, min(nEp, round((250 / 1000) * Fs)));
    plot(linspace(0, nEp / Fs, nEp), meanTraceAllConds(:, cc) - mean(meanTraceAllConds(epRefLo:epRefHi, cc)), 'linewidth', 2, 'color', [0.8 0.8 0.8]);
    set(gca, 'xtick', [], 'fontsize', 18);
    ylabel('Vm');
    title(['dV = ', num2str(abs(mean(meanTraceAllConds(1:(pulseStart / 1000) * Fs, cc))) + (max(meanTraceAllConds(((pulseStart / 1000) * Fs):((pulseStart / 1000) * Fs + ((1000 / pulseFreq) / 1000 * Fs)), cc)))), 'mV,', ' (1st pulse) ', ...
        num2str(abs(mean(meanTraceAllConds(1:(pulseStart / 1000) * Fs, cc))) + (max(meanTraceAllConds(((pulseStart / 1000) * Fs + ((1000 / pulseFreq) / 1000 * Fs)):end, cc)))), 'mV,', ' (2nd pulse)'], 'fontsize', 14);
    xticks(0:0.1:size(meanTraceAllConds, 1) / Fs)
    hold off
end

%% Save Analysis Results
directory = saveDirectory;
fileName = [voltImgTest_Analysis.mouseID, '.mat'];
save(fullfile(directory, fileName), 'voltImgTest_Analysis', '-v7.3');

TimeNow = datetime('now', 'Format', 'dd-MMM-yyyy HH:mm:ss');
disp(['finished saving at: ', char(TimeNow)])

function [imageStack, d1] = slice_loadRawTiffStack(currImgPath, rawImgNChannels)
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
    numDirsTotal = n;
    if rawImgNChannels == 2
        dirList = 1:2:numDirsTotal;
    else
        dirList = 1:numDirsTotal;
    end
    nKeep = numel(dirList);
    t.setDirectory(dirList(1));
    firstFrame = t.read();
    [d1, d2] = size(firstFrame);
    imageStack = zeros(d1, d2, nKeep, 'like', firstFrame);
    imageStack(:, :, 1) = firstFrame;
    for ki = 2:nKeep
        t.setDirectory(dirList(ki));
        imageStack(:, :, ki) = t.read();
    end
    t.close();
end

function [imageStack, numFrames] = slice_loadMcTiffStack(mcPath)
    t = Tiff(mcPath, 'r');
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
end
