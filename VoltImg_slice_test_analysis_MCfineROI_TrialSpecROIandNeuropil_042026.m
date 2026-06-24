%% Voltage Imaging Slice Test Analysis (MC + trial ROI/neuropil + new dF/F)
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
% 5) SNR / SDT panels (Figs 90–91) plus publishable side-by-side summary (Fig 92).
% 6) Analysis .mat is written under .../voltImgTest/Analysis Results/ASAP7y (see Save section).
%
% Requirements: NoRMCorre on MATLAB path (set normcorrePath below).
% Image Processing Toolbox: fibermetric, imgaussfilt, strel, imdilate, etc.

%% Initialize
clear all
close all

%% Load files and setup
% Step 1: Read the ephys file
ephysFilePath = char(uigetdir('/Volumes/phoenixinthesky/Masato/Voltage Imaging Data_Phoenix/voltImgTest/DAQ_Ephys Data'));
ephysFileDir = dir(ephysFilePath);
load([ephysFileDir(3).folder, '/', ephysFileDir(3).name]);

% NoRMCorre (edit if needed)
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

% Step 2: Imaging folder
ImgsFilePath = char(uigetdir('/Volumes/phoenixinthesky/Masato/Voltage Imaging Data_Phoenix/voltImgTest/Raw Imaging Data'));
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

% Step 4: GEVI orientation
UpOrDown = input('1 for upward GEVI (ASAP6, FORCE, etc.), 2 for downward GEVI (Jedi2Pd, ASAP5, 7y etc.) ', 's');

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

% Single-channel imaging: every TIFF IFD is one frame (no interleaved channel handling).
if isempty(imagesIndex)
    error('No TIFF files found in selected imaging folder.');
end

% Output root for motion-corrected TIFFs (created under this folder)
savePathRoot = char(uigetdir('/Volumes/phoenixinthesky/Masato/Voltage Imaging Data_Phoenix/voltImgTest/MC Imaging Data'));
expID = num2str(ExpStruct.mouseID);
mouseTag = ['voltImgSliceTest_Analysis_', expID, '_MCfineROI'];
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
voltImgTest_Analysis.mcTiffFolder = mcTiffFolder;
voltImgTest_Analysis.saveDirectory = saveDirectory;

%% Motion correction: global template + per-trial NoRMCorre (anchored)
if ~exist('NoRMCorreSetParms', 'file') || ~exist('normcorre', 'file')
    error(['NoRMCorre not found. Add NoRMCorre to the MATLAB path (set normcorrePath ', ...
        'near the top of this script) and re-run.']);
end
input('Run NoRMCorre on all trials (Ctrl+C to cancel). Requires normcorre on path.');

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
    [imageStack, ~] = slice_loadRawTiffStack(currImgPath);
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
    [imageStack, d1] = slice_loadRawTiffStack(currImgPath);
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
    mcPath = fullfile(mcTiffFolder, [baseName, '.tif']); % change to '_mc.tif' when running motion-corrected data again
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

greenScale = [zeros(256, 1), (0:255)' / 255, zeros(256, 1)];
figure(10); clf;
imagesc(meanFluorMaxDvStack); colormap(greenScale); axis equal; axis image; axis off; colorbar off;

%%%%% Rough ROI: hand select (single patched neuron)
roughRoiX = [];
roughRoiY = [];
roiHandSelect = drawfreehand;
roiHandSelectMask = createMask(roiHandSelect);
[roughRoiX, roughRoiY] = find(roiHandSelectMask);

roughRoiXAllCells = {roughRoiX};
roughRoiYAllCells = {roughRoiY};

% Global fine ROI + fallback neuropil ring on mean image (for visualization / fallback)
roiMeanMaxDvStack = zeros(size(maxDvStack, 1), size(maxDvStack, 2));
for rr = 1:numel(roughRoiX)
    roiMeanMaxDvStack(roughRoiX(rr), roughRoiY(rr)) = mean(maxDvStack(roughRoiX(rr), roughRoiY(rr), :), 3);
end
figure(11); clf;
imagesc(roiMeanMaxDvStack); axis equal; axis image; axis off; colorbar off;

roiMaxStackDouble = im2double(roiMeanMaxDvStack);
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

% Summary ROI maps (single-cell analog to mapping script figs 20–21)
globalRoiMap = zeros(size(meanFluorMaxDvStack, 1), size(meanFluorMaxDvStack, 2));
for rr = 1:numel(roiX_global)
    globalRoiMap(roiX_global(rr), roiY_global(rr)) = meanFluorMaxDvStack(roiX_global(rr), roiY_global(rr));
end
centerXY_global(1) = (min(roiX_global) + max(roiX_global)) / 2;
centerXY_global(2) = (min(roiY_global) + max(roiY_global)) / 2;

globalBkgrndMap = zeros(size(meanFluorMaxDvStack, 1), size(meanFluorMaxDvStack, 2));
for rr = 1:numel(bkgrndRoiX_global)
    globalBkgrndMap(bkgrndRoiX_global(rr), bkgrndRoiY_global(rr)) = ...
        meanFluorMaxDvStack(bkgrndRoiX_global(rr), bkgrndRoiY_global(rr));
end

figure(12); set(gcf, 'Position', [100, 100, 1800, 800]); clf;
subplot(2, 1, 1);
colormap('winter'); imagesc(meanFluorMaxDvStack); axis equal; axis image; colorbar; set(gca, 'fontsize', 12);
title('Mean fluorescence (max-dv trials)', 'FontSize', 12);
subplot(2, 1, 2);
colormap('winter'); imagesc(globalRoiMap); axis equal; axis image; colorbar; set(gca, 'fontsize', 12);
title('Global fine ROI (values from mean image)', 'FontSize', 12);
hold on;
plot(centerXY_global(2), centerXY_global(1), 'r+', 'LineWidth', 2, 'MarkerSize', 8);
dx = 5; dy = 5;
text(centerXY_global(2) + dx, centerXY_global(1) + dy, '1', ...
    'Color', 'w', 'FontSize', 18, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
hold off;

figure(13); set(gcf, 'Position', [100, 100, 1800, 300]); clf;
colormap('winter'); imagesc(globalRoiMap); axis equal; axis image; colorbar; set(gca, 'fontsize', 12);
title('Global ROI with neuropil ring (red)', 'FontSize', 12);
hold on;
[yBk2, xBk2] = find(globalBkgrndMap);
plot(xBk2, yBk2, 'r.', 'MarkerSize', 6);
hold off;

greenScale = [zeros(256, 1), (0:255)' / 255, zeros(256, 1)];
figure(14); clf;
imagesc(meanFluorMaxDvStack); colormap(greenScale); axis equal; axis image; axis off; colorbar off;
caxis(prctile(meanFluorMaxDvStack(:), [0.5 99]));

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
voltImgTest_Analysis.globalRoiMap = globalRoiMap;
voltImgTest_Analysis.globalBkgrndMap = globalBkgrndMap;
voltImgTest_Analysis.centerXY_global = centerXY_global;
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
    mcPath = fullfile(mcTiffFolder, [baseName, '.tif']); % change to '_mc.tif' when running motion-corrected data again
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
voltImgTest_Analysis.mouseID = ['voltImgSliceTest_Analysis_', num2str(ExpStruct.mouseID), '_MCfineROI'];

%% Align imaging with ephys traces
% Usage: (1) Normal pipeline — variables already exist; voltImgTest_Analysis is filled above.
%        (2) Standalone — load saved analysis, then run this section only, e.g.
%            load('/path/to/voltImgSliceTest_Analysis_*.mat', 'voltImgTest_Analysis');
%            % then select this section and run, or copy into a small script.
S = [];
if exist('voltImgTest_Analysis', 'var') && isstruct(voltImgTest_Analysis)
    S = voltImgTest_Analysis;
end
fieldsFromStruct = all(isfield(S, {'meanTraceAllConds', 'meanFiltDffAllConds', 'CIDffAllConds', 'imagingFreq'})) ...
    && isfield(S, 'ephysData') && isfield(S.ephysData, 'Fs') ...
    && isfield(S, 'pulseParams') && all(isfield(S.pulseParams, {'pulseStart', 'pulseFreq'}));
if fieldsFromStruct
    meanTraceAllConds = S.meanTraceAllConds;
    meanFiltDffAllConds = S.meanFiltDffAllConds;
    CIDffAllConds = S.CIDffAllConds;
    imagingFreq = S.imagingFreq;
    Fs = S.ephysData.Fs;
    pulseStart = S.pulseParams.pulseStart;
    pulseFreq = S.pulseParams.pulseFreq;
    nConds = size(meanFiltDffAllConds, 2);
elseif exist('meanTraceAllConds', 'var') && exist('meanFiltDffAllConds', 'var') && exist('CIDffAllConds', 'var') ...
        && exist('imagingFreq', 'var') && exist('Fs', 'var') && exist('pulseStart', 'var') && exist('pulseFreq', 'var')
    if ~exist('nConds', 'var')
        nConds = size(meanFiltDffAllConds, 2);
    end
else
    error(['Align imaging: need voltImgTest_Analysis from a saved .mat (with meanTraceAllConds, meanFiltDffAllConds, ', ...
        'CIDffAllConds, imagingFreq, ephysData.Fs, pulseParams.pulseStart/Freq), or run the pipeline through CI computation first.']);
end

% yyaxis limits shared by all exploratory + publication panels (Vm + dF/F matched to ref condition)
ccRefYUnify = min(5, nConds); % 5th condition when nConds >= 5; otherwise condition nConds
refCi = CIDffAllConds{ccRefYUnify, 1};
refDffPct = meanFiltDffAllConds(:, ccRefYUnify) * 100;
ciUpperRef = refCi(:, 2) * 100;
ciLowerRef = refCi(:, 1) * 100;
yHiRef = max(ciUpperRef, [], 'omitnan');
if isempty(yHiRef) || ~isfinite(yHiRef) || yHiRef <= 0
    yHiRef = max(max(refDffPct, [], 'omitnan'), eps);
end
yLoRef = min(ciLowerRef, [], 'omitnan');
if isempty(yLoRef) || ~isfinite(yLoRef)
    yLoRef = 0;
end
yLoRef = min(yLoRef, 0);
if yLoRef >= yHiRef
    yHiRef = yLoRef + eps;
end
nEpRef = size(meanTraceAllConds, 1);
epLoRef = max(1, min(nEpRef, round((250 / 1000) * Fs - 1000)));
epHiRef = max(epLoRef, min(nEpRef, round((250 / 1000) * Fs)));
vmRef = meanTraceAllConds(:, ccRefYUnify) - mean(meanTraceAllConds(epLoRef:epHiRef, ccRefYUnify));
vr = vmRef(isfinite(vmRef));
if isempty(vr)
    yRref = [0, 1];
else
    vmLo = min(vr);
    vmHi = max(vr);
    if vmLo > 0
        yRref = [0, vmHi];
    elseif vmHi < 0
        yRref = [vmLo, 0];
    else
        yRref = [vmLo, vmHi];
    end
end
denRref = yRref(2) - yRref(1);
if denRref > eps
    tZref = (0 - yRref(1)) / denRref;
    tZref = min(max(tZref, eps), 1 - eps);
    bRef = max(yHiRef, eps);
    aRef = -tZref * bRef / (1 - tZref);
    ephysAlignYlimLeft = [aRef, bRef];
    ephysAlignYlimRight = yRref;
else
    ephysAlignYlimLeft = [yLoRef, yHiRef];
    ephysAlignYlimRight = [0, 1];
end

ciFillCol = [0.80 0.80 0.80];
ciEdgeCol = [0.40 0.40 0.40];
ciLineCol = [0.40 0.40 0.40];
for cc = 1:nConds
    figure(30 + cc)
    clf
    hold on
    fill([linspace(0, size(CIDffAllConds{cc, 1}, 1) / imagingFreq, size(CIDffAllConds{cc, 1}, 1)), fliplr(linspace(0, size(CIDffAllConds{cc, 1}, 1) / imagingFreq, size(CIDffAllConds{cc, 1}, 1)))], ...
        [CIDffAllConds{cc, 1}(:, 1)', fliplr(CIDffAllConds{cc, 1}(:, 2)')], ciFillCol, 'EdgeColor', ciEdgeCol);
    plot(linspace(0, size(CIDffAllConds{cc, 1}, 1) / imagingFreq, size(CIDffAllConds{cc, 1}, 1)), CIDffAllConds{cc, 1}(:, 1) * 100, '-', 'linewidth', 1, 'color', ciLineCol);
    plot(linspace(0, size(CIDffAllConds{cc, 1}, 1) / imagingFreq, size(CIDffAllConds{cc, 1}, 1)), CIDffAllConds{cc, 1}(:, 2) * 100, '-', 'linewidth', 1, 'color', ciLineCol);
    dffTracePct = meanFiltDffAllConds(:, cc) * 100;
    plot(linspace(0, size(meanFiltDffAllConds, 1) / imagingFreq, size(meanFiltDffAllConds, 1)), dffTracePct, '-', 'linewidth', 2, 'color', 'g');
    set(gca, 'fontsize', 16);
    nCi = size(CIDffAllConds{cc, 1}, 1);
    ciMat = CIDffAllConds{cc, 1};
    ciUpperPct = ciMat(:, 2) * 100;
    ciLowerPct = ciMat(:, 1) * 100;
    % Upper dF/F bound follows peak upper CI (%); preliminary lower from lower CI (refined after Vm for y=0 match)
    yHiTarget = max(ciUpperPct, [], 'omitnan');
    if isempty(yHiTarget) || ~isfinite(yHiTarget) || yHiTarget <= 0
        yHiTarget = max(max(dffTracePct, [], 'omitnan'), eps);
    end
    yLoCandidate = min(ciLowerPct, [], 'omitnan');
    if isempty(yLoCandidate) || ~isfinite(yLoCandidate)
        yLoCandidate = 0;
    end
    yLoCandidate = min(yLoCandidate, 0);
    if yLoCandidate >= yHiTarget
        yHiTarget = yLoCandidate + eps;
    end
    xlim([0.2, 0.9]);
    if nCi >= 1
        ylim([yLoCandidate, yHiTarget]);
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
    vmTrace = meanTraceAllConds(:, cc) - mean(meanTraceAllConds(epRefLo:epRefHi, cc));
    plot(linspace(0, nEp / Fs, nEp), vmTrace, 'linewidth', 2, 'color', [0.8 0.8 0.8]);
    set(gca, 'xtick', [], 'fontsize', 18);
    ylabel('Vm');
    title(['dV = ', num2str(abs(mean(meanTraceAllConds(1:(pulseStart / 1000) * Fs, cc))) + (max(meanTraceAllConds(((pulseStart / 1000) * Fs):((pulseStart / 1000) * Fs + ((1000 / pulseFreq) / 1000 * Fs)), cc)))), 'mV,', ' (1st pulse) ', ...
        num2str(abs(mean(meanTraceAllConds(1:(pulseStart / 1000) * Fs, cc))) + (max(meanTraceAllConds(((pulseStart / 1000) * Fs + ((1000 / pulseFreq) / 1000 * Fs)):end, cc)))), 'mV,', ' (2nd pulse)'], 'fontsize', 14);
    xticks(0:0.1:size(meanTraceAllConds, 1) / Fs)
    % Same yyaxis limits for every condition (set from ccRefYUnify above).
    yyaxis right
    ylim(ephysAlignYlimRight);
    yyaxis left
    ylim(ephysAlignYlimLeft);
    hold off
end

%% Publication-style ephys + dF/F (square panel, light mode, scale bars)
% Same variables as "Align imaging with ephys traces" (run that section's unpack header first if
% you only loaded voltImgTest_Analysis). yyaxis limits: ephysAlignYlimLeft/Right from ref condition
% min(5, nConds), computed in the block above. Figures 130+cc (separate from exploratory Figs 31+).
% Scale bars: time + dF/F (bottom-left), Vm (far-left, right y-axis); axis off after drawing.
pubFigInch = 5;
pubLineW = 1.25;
pubDffLineW = pubLineW + 0.5;
pubCiFill = [0.90 0.90 0.90];
pubCiEdge = [0.86 0.86 0.86];
pubCiLine = [0.78 0.78 0.78];
pubDffCol = 'g';
pubVmCol = [0 0 0];
pubBarCol = [0.08 0.08 0.08];
pubFont = 'Arial';
for cc = 1:nConds
    fhPub = figure(130 + cc);
    clf(fhPub)
    set(fhPub, 'Color', 'w', 'InvertHardcopy', 'off', 'Units', 'inches', ...
        'Position', [0.5 0.5 pubFigInch pubFigInch], 'PaperPositionMode', 'auto');
    axPub = axes(fhPub, 'Color', 'w', 'FontName', pubFont, 'LineWidth', 0.75);
    hold(axPub, 'on')
    yyaxis(axPub, 'left')
    nCiPub = size(CIDffAllConds{cc, 1}, 1);
    tIm = linspace(0, nCiPub / imagingFreq, nCiPub);
    fill(axPub, [tIm, fliplr(tIm)], [CIDffAllConds{cc, 1}(:, 1)', fliplr(CIDffAllConds{cc, 1}(:, 2)')] * 100, ...
        pubCiFill, 'EdgeColor', pubCiEdge, 'LineWidth', 0.5);
    plot(axPub, tIm, CIDffAllConds{cc, 1}(:, 1) * 100, '-', 'LineWidth', 0.65, 'Color', pubCiLine);
    plot(axPub, tIm, CIDffAllConds{cc, 1}(:, 2) * 100, '-', 'LineWidth', 0.65, 'Color', pubCiLine);
    dffTracePctPub = meanFiltDffAllConds(:, cc) * 100;
    tDff = linspace(0, size(meanFiltDffAllConds, 1) / imagingFreq, size(meanFiltDffAllConds, 1));
    plot(axPub, tDff, dffTracePctPub, '-', 'LineWidth', pubDffLineW, 'Color', pubDffCol);
    ciMatPub = CIDffAllConds{cc, 1};
    ciUpperPctPub = ciMatPub(:, 2) * 100;
    ciLowerPctPub = ciMatPub(:, 1) * 100;
    yHiTargetPub = max(ciUpperPctPub, [], 'omitnan');
    if isempty(yHiTargetPub) || ~isfinite(yHiTargetPub) || yHiTargetPub <= 0
        yHiTargetPub = max(max(dffTracePctPub, [], 'omitnan'), eps);
    end
    yLoCandidatePub = min(ciLowerPctPub, [], 'omitnan');
    if isempty(yLoCandidatePub) || ~isfinite(yLoCandidatePub)
        yLoCandidatePub = 0;
    end
    yLoCandidatePub = min(yLoCandidatePub, 0);
    if yLoCandidatePub >= yHiTargetPub
        yHiTargetPub = yLoCandidatePub + eps;
    end
    xlim(axPub, [0.2, 0.9]);
    if nCiPub >= 1
        ylim(axPub, [yLoCandidatePub, yHiTargetPub]);
    end
    yline(axPub, 0, '-', 'LineWidth', 0.75, 'Color', [0.82 0.82 0.82]);

    yyaxis(axPub, 'right')
    nEpPub = size(meanTraceAllConds, 1);
    epRefLoPub = max(1, min(nEpPub, round((250 / 1000) * Fs - 1000)));
    epRefHiPub = max(epRefLoPub, min(nEpPub, round((250 / 1000) * Fs)));
    vmTracePub = meanTraceAllConds(:, cc) - mean(meanTraceAllConds(epRefLoPub:epRefHiPub, cc));
    tEp = linspace(0, nEpPub / Fs, nEpPub);
    plot(axPub, tEp, vmTracePub, 'LineWidth', pubLineW, 'Color', pubVmCol);

    yyaxis(axPub, 'right')
    ylim(axPub, ephysAlignYlimRight);
    yyaxis(axPub, 'left')
    ylim(axPub, ephysAlignYlimLeft);

    % Square plot box; light-mode tick/label colors suppressed next via axis off
    pbaspect(axPub, [1 1 1])
    set(axPub, 'Box', 'off', 'Clipping', 'off', 'XColor', 'none', 'YColor', 'none');

    xlPub = xlim(axPub);
    yyaxis(axPub, 'left')
    ylPub = ylim(axPub);
    xSpan = xlPub(2) - xlPub(1);
    ySpanL = ylPub(2) - ylPub(1);
    % Time scale (nice 100 / 200 / 500 ms)
    candDt = [0.05, 0.1, 0.2, 0.5, 1];
    dtPub = candDt(find(candDt <= 0.35 * xSpan, 1, 'last'));
    if isempty(dtPub)
        dtPub = 0.05;
    end
    padX = 0.05 * xSpan;
    padY = 0.052 * ySpanL;
    xT0 = xlPub(1) + padX;
    yT0 = ylPub(1) + padY;
    plot(axPub, [xT0, xT0 + dtPub], [yT0, yT0], '-', 'Color', pubBarCol, 'LineWidth', 1.4, 'Clipping', 'off');
    if dtPub >= 0.1
        tLab = sprintf('%g s', dtPub);
    else
        tLab = sprintf('%g ms', dtPub * 1000);
    end
    text(axPub, mean([xT0, xT0 + dtPub]), yT0 - 0.04 * ySpanL, tLab, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', 'FontSize', 8, 'Color', pubBarCol, 'FontName', pubFont);

    % dF/F scale bar (left; %), magnitude ~1/4 of positive span from 0 if possible
    yPosTop = min(ylPub(2), max(0, ylPub(1)) + 0.25 * ySpanL);
    y0D = ylPub(1) + padY;
    dDffRaw = max(0.15 * (yPosTop - max(y0D, ylPub(1))), 0.5);
    mag10 = 10^floor(log10(dDffRaw));
    nrm = dDffRaw / mag10;
    if nrm <= 1
        dDffPub = mag10;
    elseif nrm <= 2
        dDffPub = 2 * mag10;
    elseif nrm <= 5
        dDffPub = 5 * mag10;
    else
        dDffPub = 10 * mag10;
    end
    if dDffPub > ySpanL * 0.45
        dDffPub = max(mag10, ySpanL * 0.2);
    end
    xD0 = xlPub(1) + padX;
    plot(axPub, [xD0, xD0], [y0D, y0D + dDffPub], '-', 'Color', pubBarCol, 'LineWidth', 1.4, 'Clipping', 'off');
    text(axPub, xD0 + 0.025 * xSpan, y0D + 0.5 * dDffPub, sprintf('%g%%', dDffPub), 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle', 'FontSize', 8, 'Color', pubBarCol, 'FontName', pubFont);

    % Vm scale bar (right y-axis): far left, slightly lower than the dF/F % bar column
    yyaxis(axPub, 'right')
    ylRPub = ylim(axPub);
    ySpanR = ylRPub(2) - ylRPub(1);
    dVmRaw = 0.2 * ySpanR;
    mag10v = 10^floor(log10(max(dVmRaw, eps)));
    nrmv = dVmRaw / mag10v;
    if nrmv <= 1
        dVmPub = mag10v;
    elseif nrmv <= 2
        dVmPub = 2 * mag10v;
    elseif nrmv <= 5
        dVmPub = 5 * mag10v;
    else
        dVmPub = 10 * mag10v;
    end
    if dVmPub > ySpanR * 0.45
        dVmPub = max(mag10v, ySpanR * 0.2);
    end
    xVmBar = xlPub(1) + 0.028 * xSpan;
    padYVm = 0.038 * ySpanR;
    yBotVm = ylRPub(1) + padYVm;
    yTopVm = yBotVm + dVmPub;
    if yTopVm > ylRPub(2) - 0.04 * ySpanR
        yTopVm = ylRPub(2) - 0.04 * ySpanR;
        yBotVm = yTopVm - dVmPub;
        yBotVm = max(yBotVm, ylRPub(1) + 0.02 * ySpanR);
        yTopVm = yBotVm + dVmPub;
    end
    plot(axPub, [xVmBar, xVmBar], [yBotVm, yTopVm], '-', 'Color', pubBarCol, 'LineWidth', 1.4, 'Clipping', 'off');
    text(axPub, xVmBar + 0.022 * xSpan, mean([yBotVm, yTopVm]), sprintf('%g mV', dVmPub), 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle', 'FontSize', 8, 'Color', pubBarCol, 'FontName', pubFont);

    axis(axPub, 'off')
    hold(axPub, 'off')
end

%% SNR (trial-averaged trace) + detection thresholds vs ΔVm
% Baseline σ: std of filtered dF/F on prestim frames (no current injection),
% pooled across all trials and conditions (single defensible noise scale).
snrCfg.baselineFrameStart = 1;
snrCfg.baselineFrameEnd = max(2, floor((pulseStart / 1000) * imagingFreq) - 1);
snrCfg.snrThresholdList = [2, 3];
snrCfg.dPrimeTarget = 1;
baselineSeg = snrCfg.baselineFrameStart:snrCfg.baselineFrameEnd;
allBaselineFiltdff = [];
for cc = 1:nConds
    nTrHere = size(filtdffAllConds{cc}, 2);
    for tr = 1:nTrHere
        allBaselineFiltdff = [allBaselineFiltdff; filtdffAllConds{cc}(baselineSeg, tr)]; %#ok<AGROW>
    end
end
sigmaDffBaseline = std(allBaselineFiltdff, 0, 'omitnan');
nBaselineSamples = numel(allBaselineFiltdff);
if sigmaDffBaseline < eps
    sigmaDffBaseline = eps;
end

nImRows = size(meanFiltDffAllConds, 1);
nEp = size(meanTraceAllConds, 1);
iPreEnd = max(1, min(nEp, floor((pulseStart / 1000) * Fs)));

if numel(dvToTest) >= nConds
    dVmNominal = abs(dvToTest(1:nConds).');
else
    dVmNominal = nan(nConds, 1);
end

NtrialsPerCond = zeros(nConds, 1);
SNR_meanFiltDff = nan(nConds, nPulses);
peakMeanFiltDff = nan(nConds, nPulses);
dVmMeasMeanTrace = nan(nConds, nPulses);
for cc = 1:nConds
    NtrialsPerCond(cc) = size(filtdffAllConds{cc}, 2);
    vBase = mean(meanTraceAllConds(1:iPreEnd, cc), 'omitnan');
    for pp = 1:nPulses
        tMsStart = pulseStart + (1000 / pulseFreq) * (pp - 1);
        tMsEnd = pulseStart + (1000 / pulseFreq) * pp;
        idxStart = max(1, min(nImRows, floor((tMsStart / 1000) * imagingFreq)));
        idxEnd = max(idxStart, min(nImRows, floor((tMsEnd / 1000) * imagingFreq)));
        pk = max(meanFiltDffAllConds(idxStart:idxEnd, cc), [], 'omitnan');
        peakMeanFiltDff(cc, pp) = pk;
        SNR_meanFiltDff(cc, pp) = pk / sigmaDffBaseline;
        ep0 = max(1, min(nEp, floor((tMsStart / 1000) * Fs) + 1));
        ep1 = max(ep0, min(nEp, floor((tMsEnd / 1000) * Fs)));
        dVmMeasMeanTrace(cc, pp) = max(meanTraceAllConds(ep0:ep1, cc) - vBase, [], 'omitnan');
    end
end

% Collapse pulses for calibration curves (mean SNR / ΔVm across pulses in a sweep)
dVmForCalib = dVmNominal;
for cc = 1:nConds
    if isnan(dVmForCalib(cc))
        dVmForCalib(cc) = mean(dVmMeasMeanTrace(cc, :), 'omitnan');
    end
end
snrCalibVec = mean(SNR_meanFiltDff, 2, 'omitnan');
dVmCalibVec = dVmForCalib(:);
[snrOrd, ordBySnr] = sort(snrCalibVec, 'ascend', 'MissingPlacement', 'last');
dvBySnr = dVmCalibVec(ordBySnr);
validSnr = ~isnan(snrOrd) & ~isnan(dvBySnr);
snrOrd = snrOrd(validSnr);
dvBySnr = dvBySnr(validSnr);

dVmAtSnr = struct();
for kk = 1:numel(snrCfg.snrThresholdList)
    thr = snrCfg.snrThresholdList(kk);
    tag = ['snr', strrep(num2str(thr, '%g'), '.', 'p')];
    if numel(snrOrd) >= 2 && numel(unique(snrOrd)) >= 2
        dVmAtSnr.(tag) = interp1(snrOrd, dvBySnr, thr, 'linear', 'extrap');
    else
        dVmAtSnr.(tag) = nan;
    end
end

snrTrialAverage.N = NtrialsPerCond;
snrTrialAverage.nBaselineSamples = nBaselineSamples;
snrTrialAverage.sigmaDffBaseline = sigmaDffBaseline;
snrTrialAverage.baselineFrameRange = [snrCfg.baselineFrameStart, snrCfg.baselineFrameEnd];
snrTrialAverage.SNR_meanFiltDff = SNR_meanFiltDff;
snrTrialAverage.peakMeanFiltDff = peakMeanFiltDff;
snrTrialAverage.dVmNominal = dVmNominal;
snrTrialAverage.dVmMeasMeanTrace = dVmMeasMeanTrace;
snrTrialAverage.dVmAtSnrThreshold = dVmAtSnr;
snrTrialAverage.snrThresholdList = snrCfg.snrThresholdList;
snrTrialAverage.note = ['SNR = peak on trial-mean filtered dF/F divided by prestim σ (single-trial σ); ', ...
    'N = trials per condition in snrTrialAverage.N. Trial-mean noise is ~σ/√N for i.i.d. samples (upper bound).'];
snrTrialAverage.snrCalib_meanAcrossPulses = snrCalibVec;
snrTrialAverage.dVmCalib = dVmCalibVec;

figure(90); clf;
set(gcf, 'Position', [120, 120, 700, 520]);
hold on;
for cc = 1:nConds
    scatter(dVmCalibVec(cc), snrCalibVec(cc), 70, 'filled', ...
        'DisplayName', sprintf('cond %d, N=%d', cc, NtrialsPerCond(cc)));
end
for kk = 1:numel(snrCfg.snrThresholdList)
    thr = snrCfg.snrThresholdList(kk);
    fld = ['snr', strrep(num2str(thr, '%g'), '.', 'p')];
    if isfield(dVmAtSnr, fld)
        xv = dVmAtSnr.(fld);
        if ~isnan(xv)
            xline(xv, '--', sprintf('SNR=%g (~%.2f mV)', thr, xv), 'HandleVisibility', 'off');
        end
    end
end
yline(snrCfg.snrThresholdList, ':', 'LineWidth', 1, 'HandleVisibility', 'off');
grid on;
xlabel('|\DeltaV| nominal or measured (mV)');
ylabel('SNR (trial-mean peak / \sigma_{baseline})');
title(sprintf('Trial-averaged SNR; \\sigma from %d pooled baseline samples', nBaselineSamples));
legend('Location', 'northwest', 'FontSize', 9);
hold off;

axesInset = axes('Position', [0.55, 0.18, 0.32, 0.28]);
hold(axesInset, 'on');
plot(axesInset, dVmCalibVec, snrCalibVec, 'o-', 'LineWidth', 1.2, 'MarkerSize', 5);
grid(axesInset, 'on');
xlabel(axesInset, '|\DeltaV| (mV)', 'FontSize', 9);
ylabel(axesInset, 'SNR', 'FontSize', 9);
title(axesInset, 'Inset: SNR vs |\DeltaV|', 'FontSize', 10);
hold(axesInset, 'off');

%% Signal detection (single-trial): matched filter, d'', ROC AUC vs |\DeltaV|
sdtCfg.nNullDrawsPerTrial = 80;
sdtCfg.rngSeed = 42;
rng(sdtCfg.rngSeed);

dPrime_MF = nan(nConds, nPulses);
auc_MF = nan(nConds, nPulses);
dPrime_peak = nan(nConds, nPulses);
auc_peak = nan(nConds, nPulses);

for cc = 1:nConds
    nTr = size(filtdffAllConds{cc}, 2);
    if nTr < 2
        continue
    end
    for pp = 1:nPulses
        tMsStart = pulseStart + (1000 / pulseFreq) * (pp - 1);
        tMsEnd = pulseStart + (1000 / pulseFreq) * pp;
        idxStart = max(1, min(nImRows, floor((tMsStart / 1000) * imagingFreq)));
        idxEnd = max(idxStart, min(nImRows, floor((tMsEnd / 1000) * imagingFreq)));
        winLen = idxEnd - idxStart + 1;
        tmpl = meanFiltDffAllConds(idxStart:idxEnd, cc);
        tmpl = tmpl - mean(meanFiltDffAllConds(baselineSeg, cc), 'omitnan');
        nrm = norm(tmpl(:));
        if nrm < eps
            continue
        end
        u = tmpl(:) / nrm;

        sigScoresMF = zeros(nTr, 1);
        sigScoresPk = zeros(nTr, 1);
        nullScoresMF = [];
        nullScoresPk = [];
        maxStartBase = snrCfg.baselineFrameEnd - winLen + 1;
        for tr = 1:nTr
            seg = filtdffAllConds{cc}(idxStart:idxEnd, tr);
            muB = mean(filtdffAllConds{cc}(baselineSeg, tr), 'omitnan');
            segC = seg - muB;
            sigScoresMF(tr) = dot(segC(:), u);
            sigScoresPk(tr) = max(segC, [], 'omitnan');
            if maxStartBase < snrCfg.baselineFrameStart
                continue
            end
            for nn = 1:sdtCfg.nNullDrawsPerTrial
                s0 = randi([snrCfg.baselineFrameStart, maxStartBase]);
                bseg = filtdffAllConds{cc}(s0:s0 + winLen - 1, tr);
                bsegC = bseg - muB;
                nullScoresMF(end + 1, 1) = dot(bsegC(:), u); %#ok<AGROW>
                nullScoresPk(end + 1, 1) = max(bsegC, [], 'omitnan'); %#ok<AGROW>
            end
        end

        if isempty(nullScoresMF)
            continue
        end

        muS = mean(sigScoresMF, 'omitnan');
        muN = mean(nullScoresMF, 'omitnan');
        vS = var(sigScoresMF, 0, 'omitnan');
        vN = var(nullScoresMF, 0, 'omitnan');
        sigPooled = sqrt(0.5 * (vS + vN));
        if sigPooled < eps
            dPrime_MF(cc, pp) = nan;
        else
            dPrime_MF(cc, pp) = (muS - muN) / sigPooled;
        end

        muSp = mean(sigScoresPk, 'omitnan');
        muNp = mean(nullScoresPk, 'omitnan');
        vSp = var(sigScoresPk, 0, 'omitnan');
        vNp = var(nullScoresPk, 0, 'omitnan');
        sigPooledP = sqrt(0.5 * (vSp + vNp));
        if sigPooledP < eps
            dPrime_peak(cc, pp) = nan;
        else
            dPrime_peak(cc, pp) = (muSp - muNp) / sigPooledP;
        end

        labels = [zeros(numel(nullScoresMF), 1); ones(numel(sigScoresMF), 1)];
        scores = [nullScoresMF(:); sigScoresMF(:)];
        if exist('perfcurve', 'file') == 2
            [~, ~, ~, aucV] = perfcurve(labels, scores, 1);
            auc_MF(cc, pp) = aucV(1);
            scoresP = [nullScoresPk(:); sigScoresPk(:)];
            [~, ~, ~, aucVp] = perfcurve(labels, scoresP, 1);
            auc_peak(cc, pp) = aucVp(1);
        else
            auc_MF(cc, pp) = nan;
            auc_peak(cc, pp) = nan;
        end
    end
end

dPrimeCalib_MF = mean(dPrime_MF, 2, 'omitnan');
dPrimeCalib_peak = mean(dPrime_peak, 2, 'omitnan');
aucCalib_MF = mean(auc_MF, 2, 'omitnan');
dVmAtDprime = nan;
[dpOrd, ordByDp] = sort(dPrimeCalib_MF, 'ascend', 'MissingPlacement', 'last');
dvForDp = dVmCalibVec(ordByDp);
maskDp = ~isnan(dpOrd) & ~isnan(dvForDp);
dpOrd = dpOrd(maskDp);
dvForDp = dvForDp(maskDp);
if numel(dpOrd) >= 2 && numel(unique(dpOrd)) >= 2
    dVmAtDprime = interp1(dpOrd, dvForDp, snrCfg.dPrimeTarget, 'linear', 'extrap');
end

sdtSingleTrial.nNullDrawsPerTrial = sdtCfg.nNullDrawsPerTrial;
sdtSingleTrial.rngSeed = sdtCfg.rngSeed;
sdtSingleTrial.dPrime_matchedFilter = dPrime_MF;
sdtSingleTrial.auc_matchedFilter = auc_MF;
sdtSingleTrial.dPrime_peakWindow = dPrime_peak;
sdtSingleTrial.auc_peakWindow = auc_peak;
sdtSingleTrial.dPrimeCalib_matchedFilter_meanPulses = dPrimeCalib_MF;
sdtSingleTrial.dPrimeCalib_peak_meanPulses = dPrimeCalib_peak;
sdtSingleTrial.aucCalib_matchedFilter_meanPulses = aucCalib_MF;
sdtSingleTrial.dVmCalib = dVmCalibVec;
sdtSingleTrial.dVmNominal = dVmNominal;
sdtSingleTrial.dVmAtDprimeTarget = dVmAtDprime;
sdtSingleTrial.dPrimeTarget = snrCfg.dPrimeTarget;
sdtSingleTrial.note = ['d'' uses (mu_signal - mu_null)/sqrt(0.5*(var_s+var_n)) on MF scores; ', ...
    'null = random prestim windows matched to pulse length.'];

figure(91); clf;
set(gcf, 'Position', [140, 140, 720, 520]);
yyaxis left;
scatter(dVmCalibVec, dPrimeCalib_MF, 72, 'filled');
ylabel('d'' (matched filter)');
yyaxis right;
scatter(dVmCalibVec, aucCalib_MF, 52, '^', 'filled', 'MarkerFaceAlpha', 0.55);
ylabel('AUC (MF)');
grid on;
xlabel('|\DeltaV| nominal or measured (mV)');
title('Single-trial detectability vs |\DeltaV| (mean across pulses)');
hold on;
if ~isnan(dVmAtDprime)
    xline(dVmAtDprime, '--', sprintf('d''=%g (~%.2f mV)', snrCfg.dPrimeTarget, dVmAtDprime), 'HandleVisibility', 'off');
end
hold off;

axesInset2 = axes('Position', [0.52, 0.2, 0.35, 0.28]);
hold(axesInset2, 'on');
plot(axesInset2, dVmCalibVec, dPrimeCalib_MF, 'o-', 'LineWidth', 1.1);
yline(axesInset2, snrCfg.dPrimeTarget, ':');
grid(axesInset2, 'on');
xlabel(axesInset2, '|\DeltaV| (mV)', 'FontSize', 9);
ylabel(axesInset2, 'd'' (MF)', 'FontSize', 9);
title(axesInset2, 'Inset: d'' vs |\DeltaV|', 'FontSize', 10);
hold(axesInset2, 'off');

%% Publishable side-by-side: trial-averaged SNR vs single-trial d' (Fig. 92)
% Simplified panels for print: shared |ΔV| axis range, square plot boxes, readable fonts.
% Explicit light-theme colors (overrides MATLAB dark desktop defaults on axes/legend).
pubFig.fontLabel = 14;
pubFig.fontTitle = 14;
pubFig.fontTick = 12;
pubFig.fontPanel = 17;
pubFig.markerSize = 86;
pubFig.figureNumber = 92;
pubFig.bgFigure = [1 1 1];
pubFig.bgAxes = [1 1 1];
pubFig.fgText = [0 0 0];
% Okabe–Ito palette (Wong): sky blue vs vermillion — strong contrast for protan/deutan vision.
pubFig.colorDprime = [86, 180, 233] / 255;
pubFig.colorAuc = [213, 94, 0] / 255;
pubFig.colorDprimeEdge = min(1, max(0, pubFig.colorDprime * 0.45 + [0.05 0.08 0.12]));
pubFig.colorAucEdge = min(1, max(0, pubFig.colorAuc * 0.55));

figPub = figure(pubFig.figureNumber);
clf(figPub);
set(figPub, 'Color', pubFig.bgFigure, 'Position', [80, 80, 1000, 460], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
tlPub = tiledlayout(figPub, 1, 2, 'Padding', 'tight', 'TileSpacing', 'compact');
if isprop(tlPub, 'BackgroundColor')
    tlPub.BackgroundColor = pubFig.bgFigure;
end

xvCal = dVmCalibVec(:);
xvCal = xvCal(~isnan(xvCal));
if isempty(xvCal)
    xLimPub = [0, 1];
else
    dSort = diff(sort(xvCal));
    if isempty(dSort)
        dxv = max(max(xvCal) - min(xvCal), 0.5);
    else
        dxv = max([max(dSort), max(xvCal) - min(xvCal), 0.5]);
    end
    padX = max(0.07 * (max(xvCal) - min(xvCal)), 0.15 * dxv);
    xLimPub = [min(xvCal) - padX, max(xvCal) + padX];
end

markerBlack = [0 0 0];

axPubA = nexttile(tlPub, 1);
hold(axPubA, 'on');
set(axPubA, 'Color', pubFig.bgAxes, 'XColor', pubFig.fgText, 'YColor', pubFig.fgText, ...
    'TickLabelInterpreter', 'tex', 'XMinorTick', 'off', 'YMinorTick', 'off', ...
    'DefaultTextColor', pubFig.fgText, 'XGrid', 'off', 'YGrid', 'off');
for cc = 1:nConds
    scatter(axPubA, dVmCalibVec(cc), snrCalibVec(cc), pubFig.markerSize, markerBlack, 'filled', ...
        'MarkerEdgeColor', markerBlack, 'LineWidth', 0.35);
end
yline(axPubA, snrCfg.snrThresholdList(1), '--', 'Color', [0 0 0], 'LineWidth', 1.4, ...
    'Label', sprintf('SNR = %g', snrCfg.snrThresholdList(1)), 'LabelHorizontalAlignment', 'left', ...
    'FontSize', pubFig.fontTick - 1);
if numel(snrCfg.snrThresholdList) > 1
    yline(axPubA, snrCfg.snrThresholdList(2), ':', 'Color', [0 0 0], 'LineWidth', 1.2, ...
        'Label', sprintf('SNR = %g', snrCfg.snrThresholdList(2)), 'LabelHorizontalAlignment', 'left', ...
        'FontSize', pubFig.fontTick - 1);
end
xlim(axPubA, xLimPub);
ylA = ylim(axPubA);
if ylA(2) > ylA(1)
    ylim(axPubA, [max(0, ylA(1) - 0.05 * diff(ylA)), ylA(2) + 0.08 * diff(ylA)]);
end
xlabel(axPubA, '|\DeltaV| (mV)', 'FontSize', pubFig.fontLabel, 'Color', pubFig.fgText);
ylabel(axPubA, 'SNR (trial-mean peak / \sigma_{pre})', 'FontSize', pubFig.fontLabel, 'Color', pubFig.fgText);
title(axPubA, 'Trial-averaged optical SNR', 'FontSize', pubFig.fontTitle, 'FontWeight', 'normal', 'Color', pubFig.fgText);
set(axPubA, 'FontSize', pubFig.fontTick, 'LineWidth', 1.1, 'TickDir', 'out', 'Box', 'off');
axis(axPubA, 'square');
hold(axPubA, 'off');
% text(axPubA, 0.03, 0.97, 'A', 'Units', 'normalized', 'FontSize', pubFig.fontPanel, ...
    % 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Color', pubFig.fgText);

dpB = dPrimeCalib_MF(:);
dpB = dpB(~isnan(dpB));
if isempty(dpB)
    yLimPubB = [0, 1];
else
    padY = max(0.12 * (max(dpB) - min(dpB)), 0.15);
    yLimPubB = [min(dpB) - padY, max(dpB) + padY];
    if snrCfg.dPrimeTarget < yLimPubB(1) || snrCfg.dPrimeTarget > yLimPubB(2)
        yLimPubB(1) = min(yLimPubB(1), snrCfg.dPrimeTarget - 0.2);
        yLimPubB(2) = max(yLimPubB(2), snrCfg.dPrimeTarget + 0.2);
    end
end
aucPub = aucCalib_MF(:);
aucPub = aucPub(~isnan(aucPub));
if isempty(aucPub)
    yLimPubBR = [0.42, 1];
else
    loA = min(aucPub);
    hiA = max(aucPub);
    padA = max(0.04 * (hiA - loA), 0.025);
    yLimPubBR = [max(0.4, loA - padA), min(1.0, hiA + padA)];
    if diff(yLimPubBR) < 0.12
        yLimPubBR(2) = min(1, yLimPubBR(1) + 0.15);
    end
end

axPubB = nexttile(tlPub, 2);
hold(axPubB, 'on');
yyaxis(axPubB, 'left');
set(axPubB, 'Color', pubFig.bgAxes, 'XColor', pubFig.fgText, 'YColor', pubFig.colorDprime, ...
    'TickLabelInterpreter', 'tex', 'XMinorTick', 'off', 'YMinorTick', 'off', ...
    'DefaultTextColor', pubFig.fgText, 'XGrid', 'off', 'YGrid', 'off');
for cc = 1:nConds
    scatter(axPubB, dVmCalibVec(cc), dPrimeCalib_MF(cc), pubFig.markerSize, pubFig.colorDprime, 'filled', ...
        'MarkerEdgeColor', pubFig.colorDprimeEdge, 'LineWidth', 0.45);
end
yline(axPubB, snrCfg.dPrimeTarget, '--', 'Color', pubFig.colorDprimeEdge, 'LineWidth', 1.4, ...
    'Label', sprintf('d'' = %g', snrCfg.dPrimeTarget), 'LabelHorizontalAlignment', 'left', ...
    'FontSize', pubFig.fontTick - 1);
xlim(axPubB, xLimPub);
ylim(axPubB, yLimPubB);
ylabel(axPubB, 'd'' (matched filter, mean pulses)', 'FontSize', pubFig.fontLabel, 'Color', pubFig.colorDprime);

yyaxis(axPubB, 'right');
set(axPubB, 'YColor', pubFig.colorAuc);
for cc = 1:nConds
    % Same scatter pattern as Fig. 91: marker char then 'filled'; color via name-value (not numeric between '^' and 'filled').
    scatter(axPubB, dVmCalibVec(cc), aucCalib_MF(cc), round(pubFig.markerSize * 0.62), '^', 'filled', ...
        'MarkerEdgeColor', pubFig.colorAucEdge, 'MarkerFaceColor', pubFig.colorAuc, ...
        'MarkerFaceAlpha', 0.92, 'LineWidth', 0.45);
end
ylim(axPubB, yLimPubBR);
ylabel(axPubB, 'AUC (matched filter)', 'FontSize', pubFig.fontLabel, 'Color', pubFig.colorAuc);

yyaxis(axPubB, 'left');
xlabel(axPubB, '|\DeltaV| (mV)', 'FontSize', pubFig.fontLabel, 'Color', pubFig.fgText);
title(axPubB, 'Single-trial detectability', 'FontSize', pubFig.fontTitle, 'FontWeight', 'normal', 'Color', pubFig.fgText);
set(axPubB, 'FontSize', pubFig.fontTick, 'LineWidth', 1.1, 'TickDir', 'out', 'Box', 'off');
axis(axPubB, 'square');
hold(axPubB, 'off');
% text(axPubB, 0.03, 0.97, 'B', 'Units', 'normalized', 'FontSize', pubFig.fontPanel, ...
    % 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Color', pubFig.fgText);

sgtitle(figPub, 'Optical reporting vs stimulus amplitude', 'FontSize', pubFig.fontTitle + 1.5, ...
    'FontWeight', 'bold', 'Color', pubFig.fgText);

pubFig.xLimShared = xLimPub;
pubFig.yLimPanelB_dPrime = yLimPubB;
pubFig.yLimPanelB_auc = yLimPubBR;
pubFig.note = 'Export: print(figPub, ''-dpdf'', ''-r300'', filename) or exportgraphics for vector output.';

voltImgTest_Analysis.snrTrialAverage = snrTrialAverage;
voltImgTest_Analysis.sdtSingleTrial = sdtSingleTrial;
voltImgTest_Analysis.publishSnrDprimeSideBySide = pubFig;

disp('--- SNR / SDT summary (trial-mean SNR uses single-trial prestim σ) ---');
disp(['  Prestim σ (filtered dF/F) = ', num2str(sigmaDffBaseline, '%.5g'), ...
    ' ; pooled baseline samples = ', num2str(nBaselineSamples)]);
for kk = 1:numel(snrCfg.snrThresholdList)
    thr = snrCfg.snrThresholdList(kk);
    fld = ['snr', strrep(num2str(thr, '%g'), '.', 'p')];
    if isfield(dVmAtSnr, fld)
        disp(['  Interpolated |ΔV| at SNR=', num2str(thr, '%g'), '  ≈ ', num2str(dVmAtSnr.(fld), '%.3f'), ' mV']);
    end
end
disp(['  Interpolated |ΔV| at d''=', num2str(snrCfg.dPrimeTarget, '%g'), ' (MF, mean pulses) ≈ ', num2str(dVmAtDprime, '%.3f'), ' mV']);
disp(['  Publishable side-by-side SNR vs d'': figure ', num2str(pubFig.figureNumber)]);

%% Save Analysis Results
% voltImgTest_Analysis includes: traces/ROI/MC paths, peak & CI summaries, SNR trial-average
% (snrTrialAverage), SDT single-trial (sdtSingleTrial), publish panel params (publishSnrDprimeSideBySide),
% and analysis config structs snrCfg / sdtCfg for reproducibility.
analysisResultsDirectory = fullfile('/Volumes/phoenixinthesky/Masato/Voltage Imaging Data_Phoenix', ...
    'voltImgTest', 'Analysis Results', 'ASAP7y');
if ~exist(analysisResultsDirectory, 'dir')
    mkdir(analysisResultsDirectory);
end

voltImgTest_Analysis.snrCfg = snrCfg;
voltImgTest_Analysis.sdtCfg = sdtCfg;
voltImgTest_Analysis.analysisResultsDirectory = analysisResultsDirectory;
% saveDirectory (set earlier) remains the motion-correction / per-mouse output root; .mat is under analysisResultsDirectory.

fileName = [voltImgTest_Analysis.mouseID, '.mat'];
matFullPath = fullfile(analysisResultsDirectory, fileName);
voltImgTest_Analysis.savedMatFullPath = matFullPath;

save(matFullPath, 'voltImgTest_Analysis', '-v7.3');

TimeNow = datetime('now', 'Format', 'dd-MMM-yyyy HH:mm:ss');
disp(['finished saving at: ', char(TimeNow)]);
disp(['  Saved: ', matFullPath]);

function [imageStack, d1] = slice_loadRawTiffStack(currImgPath)
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
    if numDirsTotal < 1
        error('No image directories in TIFF: %s', currImgPath);
    end
    t.setDirectory(1);
    firstFrame = t.read();
    [d1, d2] = size(firstFrame);
    imageStack = zeros(d1, d2, numDirsTotal, 'like', firstFrame);
    imageStack(:, :, 1) = firstFrame;
    for ki = 2:numDirsTotal
        t.setDirectory(ki);
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
