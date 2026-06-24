function trialMeanPlane = VoltImg_mapping_parallel_mcPass2OneTrial(tt, ctx)
%VOLTIMG_MAPPING_PARALLEL_MCPASS2ONETRIAL  Motion-correct one trial (NoRMCorre pass 2).
%   trialMeanPlane = HxW mean of the motion-corrected stack; all NaN if trial is in
%   ctx.excludeTrials (still writes _mc.tif like the serial script).

currImgPath = fullfile(ctx.ImgfolderContents(ctx.imagesIndex(tt)).folder, ...
    ctx.ImgfolderContents(ctx.imagesIndex(tt)).name);

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

if ctx.rawImgNChannels == 2
    dirList = 1:2:numDirsTotal;
else
    dirList = 1:numDirsTotal;
end
nKeep = numel(dirList);

t.setDirectory(dirList(1));
firstFrame = t.read();
[H, W] = size(firstFrame);

imageStack = zeros(H, W, nKeep, 'like', firstFrame);
imageStack(:, :, 1) = firstFrame;

for ki = 2:nKeep
    t.setDirectory(dirList(ki));
    imageStack(:, :, ki) = t.read();
end
t.close();

if ctx.useLaserRowArtifactFilter
    imageStack = VoltImg_applyLaserRowArtifactToStack(imageStack, ctx.laserArtifactGateColFirst, ...
        ctx.laserArtifactGateColLast, ctx.laserArtifactThreshMode, ctx.laserArtifactThreshParam, ctx.laserArtifactMcMode);
end

if ctx.mcUseGateColumnsOnly
    gateColFirstTrial = max(1, min(ctx.laserArtifactGateColFirst, ctx.laserArtifactGateColLast));
    gateColLastTrial = min(size(imageStack, 2), max(ctx.laserArtifactGateColFirst, ctx.laserArtifactGateColLast));
    Ygate = single(imageStack(:, gateColFirstTrial:gateColLastTrial, :));
    [d1Gate, d2Gate, ~] = size(Ygate);
    options_gate = NoRMCorreSetParms('d1', d1Gate, 'd2', d2Gate, 'bin_width', 15, 'max_shift', 4, 'us_fac', 50, 'init_batch', 1);
    [~, shifts] = normcorre(Ygate, options_gate, ctx.globalTemplate);

    Yfull = single(imageStack);
    [d1, d2, ~] = size(Yfull);
    options_full = NoRMCorreSetParms('d1', d1, 'd2', d2, 'bin_width', 15, 'max_shift', 4, 'us_fac', 50, 'init_batch', 1);
    imageStack_mc = apply_shifts(Yfull, shifts, options_full);
else
    Y = single(imageStack);
    [d1, d2, ~] = size(Y);
    options_rigid = NoRMCorreSetParms('d1', d1, 'd2', d2, 'bin_width', 15, 'max_shift', 4, 'us_fac', 50, 'init_batch', 1);
    [M_mc, ~] = normcorre(Y, options_rigid, ctx.globalTemplate);
    imageStack_mc = M_mc;
end

[~, rawName, ~] = fileparts(ctx.ImgfolderContents(ctx.imagesIndex(tt)).name);
mcName = [rawName, '_mc.tif'];
mcPath = fullfile(ctx.mcTiffFolder, mcName);

if ctx.laserArtifactMcSecondSweepForDff
    badRowMask = VoltImg_laserRowArtifact_badRowMaskStack(single(imageStack_mc), ...
        ctx.laserArtifactGateColFirst, ctx.laserArtifactGateColLast, ...
        ctx.laserArtifactThreshMode, ctx.laserArtifactThreshParam);
    badRowsMatPath = fullfile(ctx.mcTiffFolder, [rawName, '_mc_badRows.mat']);
    try
        save(badRowsMatPath, 'badRowMask', 'laserArtifactGateColFirst', 'laserArtifactGateColLast', ...
            'laserArtifactThreshMode', 'laserArtifactThreshParam', '-v7.3');
    catch
        save(badRowsMatPath, 'badRowMask', 'laserArtifactGateColFirst', 'laserArtifactGateColLast', ...
            'laserArtifactThreshMode', 'laserArtifactThreshParam');
    end
end

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

if ~isfield(ctx, 'maxDvFrameCap') || isempty(ctx.maxDvFrameCap)
    ctx.maxDvFrameCap = inf;
end
if ~isfield(ctx, 'maxDvTrialMask') || isempty(ctx.maxDvTrialMask)
    ctx.maxDvTrialMask = true(1, numel(ctx.imagesIndex));
end

nZ = size(imageStack_mc, 3);
nCap = min(nZ, ctx.maxDvFrameCap);
mcForMean = imageStack_mc(:, :, 1:nCap);
imageStackMean_mc = mean(mcForMean, 3);

if ismember(tt, ctx.excludeTrials)
    trialMeanPlane = nan(size(imageStackMean_mc), 'like', imageStackMean_mc);
elseif tt > numel(ctx.maxDvTrialMask) || ~ctx.maxDvTrialMask(tt)
    trialMeanPlane = nan(size(imageStackMean_mc), 'like', imageStackMean_mc);
else
    trialMeanPlane = imageStackMean_mc;
end

end
