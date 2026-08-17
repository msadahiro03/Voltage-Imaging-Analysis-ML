function [roughRoiXAllCells, roughRoiYAllCells, autoRoiReport] = ...
    voltimg_autoRoi_cellpose(meanFluorMaxDvStack, roiCfg, autoRoiDir)
%VOLTIMG_AUTOROI_CELLPOSE Automatic rough-ROI detection via Cellpose.
%   Drop-in replacement for the hand-drawn drawfreehand rough ROIs in the
%   MultiCell mapping pipeline. Writes meanFluorMaxDvStack to a TIFF, runs
%   cellpose_wrapper.py in a subprocess, applies quality filters and a
%   dilation that emulates the generous hand-drawn margins, and returns
%   per-cell pixel indices in the exact format the pipeline expects:
%
%   roughRoiXAllCells{nn} - row indices  (Nx1 double, 1-based)
%   roughRoiYAllCells{nn} - col indices  (Nx1 double, 1-based)
%   (identical to [x, y] = find(createMask(drawfreehand)))
%
%   autoRoiReport carries all params, per-label decisions, and timing; it is
%   also written to autoRoiDir as autoRoi_report.json/.mat together with a
%   QC overlay figure (autoRoi_overlay.png/.fig) and the Cellpose log.
%   QC outputs never block the run.

tTotal = tic;

if ~exist(autoRoiDir, 'dir'), mkdir(autoRoiDir); end

% ---- Config with defaults --------------------------------------------------
thisDir   = fileparts(mfilename('fullpath'));
pythonExe = voltimg_getcfg(roiCfg, 'pythonExe', fullfile(thisDir, '.venv_cellpose', 'bin', 'python'));
wrapper   = voltimg_getcfg(roiCfg, 'wrapperPath', fullfile(thisDir, 'cellpose_wrapper.py'));

model              = voltimg_getcfg(roiCfg, 'model', 'cpsam');
diameter           = voltimg_getcfg(roiCfg, 'diameter', 30);
flowThreshold      = voltimg_getcfg(roiCfg, 'flowThreshold', 0.4);
cellprobThreshold  = voltimg_getcfg(roiCfg, 'cellprobThreshold', 0.0);
useGpu             = voltimg_getcfg(roiCfg, 'useGpu', false);

minAreaPx        = voltimg_getcfg(roiCfg, 'minAreaPx', 150);
maxAreaPx        = voltimg_getcfg(roiCfg, 'maxAreaPx', 5000);
minSeparationPx  = voltimg_getcfg(roiCfg, 'minSeparationPx', 15);
maxCells         = voltimg_getcfg(roiCfg, 'maxCells', Inf);
excludeBorderPx  = voltimg_getcfg(roiCfg, 'excludeBorderPx', 0);
dilateRadiusPx   = voltimg_getcfg(roiCfg, 'dilateRadiusPx', 3);   % 3 px best matches hand-drawn margins (validation sweep 2026-08)

% ---- 1. Sanitize image and write exchange TIFF ------------------------------
img = double(meanFluorMaxDvStack);
finiteMask = isfinite(img);
if ~any(finiteMask(:))
    error('VoltImg:autoRoi', 'meanFluorMaxDvStack has no finite pixels.');
end
img(~finiteMask) = min(img(finiteMask));
[H, W] = size(img);

inputTif = fullfile(autoRoiDir, 'input_mean_image.tif');
imwrite(uint16(rescale(img) * 65535), inputTif);

% ---- 2. Write params.json and run the wrapper -------------------------------
params = struct( ...
    'model', model, ...
    'diameter', diameter, ...
    'flow_threshold', flowThreshold, ...
    'cellprob_threshold', cellprobThreshold, ...
    'use_gpu', logical(useGpu));
paramsFile = fullfile(autoRoiDir, 'params.json');
fid = fopen(paramsFile, 'w');
fwrite(fid, jsonencode(params));
fclose(fid);

cmd = sprintf('"%s" "%s" --image "%s" --params "%s" --outdir "%s"', ...
    pythonExe, wrapper, inputTif, paramsFile, autoRoiDir);
tCellpose = tic;
[status, cmdOut] = system(cmd);
cellposeElapsed = toc(tCellpose);

logFile = fullfile(autoRoiDir, 'cellpose_log.txt');
fid = fopen(logFile, 'w');
fwrite(fid, sprintf('%s\n\nexit status: %d\nelapsed: %.1f s\n\n%s', ...
    cmd, status, cellposeElapsed, cmdOut));
fclose(fid);

masksFile = fullfile(autoRoiDir, 'masks.tif');
if status ~= 0 || ~exist(masksFile, 'file')
    tail = cmdOut(max(1, numel(cmdOut) - 1500):end);
    error('VoltImg:autoRoi', ...
        ['Cellpose wrapper failed (exit %d). Log: %s\n' ...
         '--- output tail ---\n%s'], status, logFile, tail);
end

% ---- 3. Read label mask ------------------------------------------------------
labelMask = imread(masksFile);
if ~isequal(size(labelMask), [H, W])
    error('VoltImg:autoRoi', 'masks.tif size %s does not match image %s.', ...
        mat2str(size(labelMask)), mat2str([H, W]));
end
labelMask = double(labelMask);
nRaw = max(labelMask(:));

cellposeInfo = struct();
outJson = fullfile(autoRoiDir, 'cellpose_output.json');
if exist(outJson, 'file')
    try
        cellposeInfo = jsondecode(fileread(outJson));
    catch
    end
end

% ---- 4. Per-label stats + quality-filter chain (tight masks) -----------------
labels = struct('label', {}, 'areaPx', {}, 'centroidRow', {}, 'centroidCol', {}, ...
    'meanIntensity', {}, 'accepted', {}, 'rejectReason', {});
props = regionprops(labelMask, img, 'Area', 'Centroid', 'MeanIntensity', 'PixelIdxList');

for ll = 1:nRaw
    labels(ll).label = ll;
    if ll > numel(props) || isempty(props(ll).PixelIdxList)
        labels(ll).areaPx = 0;
        labels(ll).accepted = false;
        labels(ll).rejectReason = 'empty label';
        continue
    end
    labels(ll).areaPx = props(ll).Area;
    labels(ll).centroidRow = props(ll).Centroid(2);   % regionprops Centroid = [x y]
    labels(ll).centroidCol = props(ll).Centroid(1);
    labels(ll).meanIntensity = props(ll).MeanIntensity;
    labels(ll).accepted = true;
    labels(ll).rejectReason = '';

    % 4a. border exclusion
    if excludeBorderPx > 0
        [pr, pc] = ind2sub([H, W], props(ll).PixelIdxList);
        if any(pr <= excludeBorderPx | pr > H - excludeBorderPx | ...
               pc <= excludeBorderPx | pc > W - excludeBorderPx)
            labels(ll).accepted = false;
            labels(ll).rejectReason = sprintf('within %d px of border', excludeBorderPx);
            continue
        end
    end

    % 4b. area
    if props(ll).Area < minAreaPx
        labels(ll).accepted = false;
        labels(ll).rejectReason = sprintf('area %d < minAreaPx %d', props(ll).Area, minAreaPx);
    elseif props(ll).Area > maxAreaPx
        labels(ll).accepted = false;
        labels(ll).rejectReason = sprintf('area %d > maxAreaPx %d', props(ll).Area, maxAreaPx);
    end
end

% 4c. min separation: greedy keep by descending mean intensity
survivors = find([labels.accepted]);
[~, order] = sort([labels(survivors).meanIntensity], 'descend');
keptCentroids = zeros(0, 2);
for ii = order(:)'
    ll = survivors(ii);
    c = [labels(ll).centroidRow, labels(ll).centroidCol];
    if ~isempty(keptCentroids) && ...
            any(sqrt(sum((keptCentroids - c).^2, 2)) < minSeparationPx)
        labels(ll).accepted = false;
        labels(ll).rejectReason = sprintf('centroid < minSeparationPx (%g px) from a kept cell', minSeparationPx);
    else
        keptCentroids(end+1, :) = c; %#ok<AGROW>
    end
end

% 4d. maxCells cap: keep top-N by mean intensity
survivors = find([labels.accepted]);
if numel(survivors) > maxCells
    [~, order] = sort([labels(survivors).meanIntensity], 'descend');
    for ii = order(maxCells+1:end)
        ll = survivors(ii);
        labels(ll).accepted = false;
        labels(ll).rejectReason = sprintf('beyond maxCells cap (%g)', maxCells);
    end
end

% ---- 5. Deterministic ordering: centroid row, then column --------------------
accepted = find([labels.accepted]);
[~, order] = sortrows([[labels(accepted).centroidRow]' [labels(accepted).centroidCol]']);
accepted = accepted(order);
nCellsOut = numel(accepted);

% ---- 6-7. Dilate and emit the hand-drawn contract ----------------------------
roughRoiXAllCells = cell(nCellsOut, 1);
roughRoiYAllCells = cell(nCellsOut, 1);
dilatedMasks = cell(nCellsOut, 1);
se = strel('disk', dilateRadiusPx);
for nn = 1:nCellsOut
    ll = accepted(nn);
    cellMask = false(H, W);
    cellMask(props(ll).PixelIdxList) = true;
    if dilateRadiusPx > 0
        cellMask = imdilate(cellMask, se);
    end
    dilatedMasks{nn} = cellMask;
    [roughRoiX, roughRoiY] = find(cellMask);
    roughRoiXAllCells{nn} = roughRoiX;
    roughRoiYAllCells{nn} = roughRoiY;
end

% ---- 9. Report ----------------------------------------------------------------
autoRoiReport = struct();
autoRoiReport.paramsUsed = struct('model', model, 'diameter', diameter, ...
    'flowThreshold', flowThreshold, 'cellprobThreshold', cellprobThreshold, ...
    'useGpu', logical(useGpu), 'minAreaPx', minAreaPx, 'maxAreaPx', maxAreaPx, ...
    'minSeparationPx', minSeparationPx, 'maxCells', maxCells, ...
    'excludeBorderPx', excludeBorderPx, 'dilateRadiusPx', dilateRadiusPx, ...
    'pythonExe', pythonExe, 'wrapperPath', wrapper);
if isfield(cellposeInfo, 'cellposeVersion')
    autoRoiReport.cellposeVersion = cellposeInfo.cellposeVersion;
end
autoRoiReport.nRaw = nRaw;
autoRoiReport.nAccepted = nCellsOut;
autoRoiReport.acceptedLabelsInOrder = [labels(accepted).label];
autoRoiReport.labels = labels;
autoRoiReport.cellposeElapsedSec = cellposeElapsed;
autoRoiReport.totalElapsedSec = toc(tTotal);
autoRoiReport.autoRoiDir = autoRoiDir;

try
    fid = fopen(fullfile(autoRoiDir, 'autoRoi_report.json'), 'w');
    fwrite(fid, jsonencode(autoRoiReport));
    fclose(fid);
    save(fullfile(autoRoiDir, 'autoRoi_report.mat'), 'autoRoiReport');
catch reportErr
    warning('VoltImg:autoRoi', 'Could not save autoRoi report: %s', reportErr.message);
end

% ---- 8. QC overlay (never blocks) ---------------------------------------------
try
    fig = figure('Visible', 'off', 'Position', [100, 100, 1400, 900]);
    imagesc(img); colormap(gray); axis image; hold on;
    % raw Cellpose boundaries (yellow)
    for ll = 1:nRaw
        b = bwboundaries(labelMask == ll, 'noholes');
        for kk = 1:numel(b)
            plot(b{kk}(:, 2), b{kk}(:, 1), 'y-', 'LineWidth', 0.75);
        end
    end
    % rejected (red, with reason)
    for ll = find(~[labels.accepted])
        if isempty(labels(ll).rejectReason) || labels(ll).areaPx == 0, continue; end
        b = bwboundaries(labelMask == ll, 'noholes');
        for kk = 1:numel(b)
            plot(b{kk}(:, 2), b{kk}(:, 1), 'r-', 'LineWidth', 1.2);
        end
        text(labels(ll).centroidCol, labels(ll).centroidRow, labels(ll).rejectReason, ...
            'Color', 'r', 'FontSize', 7, 'Interpreter', 'none', ...
            'HorizontalAlignment', 'center');
    end
    % accepted dilated ROIs (green, numbered)
    for nn = 1:nCellsOut
        ll = accepted(nn);
        b = bwboundaries(dilatedMasks{nn}, 'noholes');
        for kk = 1:numel(b)
            plot(b{kk}(:, 2), b{kk}(:, 1), 'g-', 'LineWidth', 1.5);
        end
        text(labels(ll).centroidCol, labels(ll).centroidRow, sprintf('%d', nn), ...
            'Color', 'g', 'FontSize', 12, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center');
    end
    title(sprintf(['AutoROI: %d raw \\rightarrow %d accepted | model=%s diam=%g ' ...
        'dilate=%dpx area=[%g %g] sep=%gpx'], nRaw, nCellsOut, model, diameter, ...
        dilateRadiusPx, minAreaPx, maxAreaPx, minSeparationPx), 'Interpreter', 'tex');
    saveas(fig, fullfile(autoRoiDir, 'autoRoi_overlay.png'));
    savefig(fig, fullfile(autoRoiDir, 'autoRoi_overlay.fig'));
    close(fig);
catch qcErr
    warning('VoltImg:autoRoi', 'QC overlay failed (non-fatal): %s', qcErr.message);
end

% ---- 10. Zero-cell handling -----------------------------------------------------
if nCellsOut == 0
    reasons = '';
    for ll = 1:numel(labels)
        if ~isempty(labels(ll).rejectReason)
            reasons = sprintf('%s\n  label %d: %s', reasons, labels(ll).label, labels(ll).rejectReason);
        end
    end
    if nRaw == 0
        error('VoltImg:autoRoi', ...
            'Cellpose found 0 masks. See QC overlay + log in %s', autoRoiDir);
    else
        error('VoltImg:autoRoi', ...
            'All %d Cellpose masks rejected by quality filters:%s\nSee QC overlay in %s', ...
            nRaw, reasons, autoRoiDir);
    end
end

fprintf('AutoROI: %d/%d Cellpose masks accepted (%.1f s total). QC: %s\n', ...
    nCellsOut, nRaw, autoRoiReport.totalElapsedSec, autoRoiDir);
end
