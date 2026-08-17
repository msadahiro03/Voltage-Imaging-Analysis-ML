function results = voltimg_validate_autoRoi(matPath, outDir, dilateRadii)
%VOLTIMG_VALIDATE_AUTOROI Validate Cellpose auto-ROIs against hand-drawn ROIs.
%   results = VOLTIMG_VALIDATE_AUTOROI(matPath, outDir, dilateRadii)
%
%   Loads meanFluorMaxDvStack + hand-drawn roughRoiX/YAllCells (ground truth)
%   from a saved voltMapping analysis .mat (v7.3; partial load via matfile),
%   runs voltimg_autoRoi_cellpose at each dilation radius, greedy-matches
%   detections to hand-drawn cells by IoU, and reports per-cell IoU, centroid
%   distance, hit rate, and false-positive count. Saves a per-radius overlay
%   (green = auto dilated, magenta = hand-drawn) and a summary .mat/.txt.
%
%   Defaults: the repo sample data, <auto_roi>/validation_out, [0 3 6 9 12].

thisDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(thisDir);
if nargin < 1 || isempty(matPath)
    matPath = fullfile(repoRoot, 'MC Imaging Data Sample', 'Parity_Check_MLAnalysisData', ...
        ['voltMapping_Analysis_SCNNCRE_hSynyASAP7Kv_DIOChrom2s_IC_InVivo_MS26_21_A7P_Chrome2s_', ...
         '050726_FOV1_2PMapping_Day1_MultiCellAnalysis_MCfineROI_laserRowArtifact.mat']);
end
if nargin < 2 || isempty(outDir)
    outDir = fullfile(thisDir, 'validation_out');
end
if nargin < 3 || isempty(dilateRadii)
    dilateRadii = [0 3 6 9 12];
end
if ~exist(outDir, 'dir'), mkdir(outDir); end

% ---- Ground truth (partial load; the .mat can be >1 GB) ---------------------
m = matfile(matPath);
img = m.meanFluorMaxDvStack;
gtX = m.roughRoiXAllCells;   % hand-drawn: row indices per cell
gtY = m.roughRoiYAllCells;   % hand-drawn: col indices per cell
[H, W] = size(img);
nGt = numel(gtX);
gtMasks = cell(nGt, 1);
gtCentroids = zeros(nGt, 2);
for gg = 1:nGt
    mk = false(H, W);
    mk(sub2ind([H, W], gtX{gg}, gtY{gg})) = true;
    gtMasks{gg} = mk;
    gtCentroids(gg, :) = [mean(gtX{gg}), mean(gtY{gg})];
end
fprintf('Ground truth: %d hand-drawn cells on %dx%d image (%s)\n', nGt, H, W, matPath);

% ---- Run detector at each dilation radius ------------------------------------
hitIouThresh = 0.3;
hitCentroidPx = 15;
results = struct('dilateRadiusPx', {}, 'nDetected', {}, 'nHits', {}, ...
    'nMisses', {}, 'nFalsePos', {}, 'meanIouMatched', {}, 'perCell', {});

for rr = 1:numel(dilateRadii)
    rad = dilateRadii(rr);
    roiCfg = struct('dilateRadiusPx', rad);   % all other params = defaults
    runDir = fullfile(outDir, sprintf('AutoROI_dilate%02d', rad));
    [autoX, autoY, report] = voltimg_autoRoi_cellpose(img, roiCfg, runDir);
    nDet = numel(autoX);

    detMasks = cell(nDet, 1);
    for dd = 1:nDet
        mk = false(H, W);
        mk(sub2ind([H, W], autoX{dd}, autoY{dd})) = true;
        detMasks{dd} = mk;
    end

    % IoU matrix (GT x detections), then greedy matching (best IoU first)
    iouMat = zeros(nGt, nDet);
    for gg = 1:nGt
        for dd = 1:nDet
            iouMat(gg, dd) = nnz(gtMasks{gg} & detMasks{dd}) / nnz(gtMasks{gg} | detMasks{dd});
        end
    end
    matchDet = zeros(nGt, 1); matchIou = zeros(nGt, 1);
    iouWork = iouMat;
    while true
        [bestIou, idx] = max(iouWork(:));
        if isempty(bestIou) || bestIou <= 0, break; end
        [gg, dd] = ind2sub(size(iouWork), idx);
        matchDet(gg) = dd; matchIou(gg) = bestIou;
        iouWork(gg, :) = -1; iouWork(:, dd) = -1;
    end

    perCell = struct('gtCell', {}, 'matchedDet', {}, 'iou', {}, 'centroidDistPx', {}, 'hit', {});
    nHits = 0;
    for gg = 1:nGt
        perCell(gg).gtCell = gg;
        perCell(gg).matchedDet = matchDet(gg);
        perCell(gg).iou = matchIou(gg);
        if matchDet(gg) > 0
            dd = matchDet(gg);
            dCent = [mean(autoX{dd}), mean(autoY{dd})];
            perCell(gg).centroidDistPx = norm(gtCentroids(gg, :) - dCent);
        else
            perCell(gg).centroidDistPx = Inf;
        end
        perCell(gg).hit = matchIou(gg) > hitIouThresh || perCell(gg).centroidDistPx < hitCentroidPx;
        nHits = nHits + perCell(gg).hit;
    end
    nFalsePos = nDet - nnz(matchDet);

    results(rr).dilateRadiusPx = rad;
    results(rr).nDetected = nDet;
    results(rr).nHits = nHits;
    results(rr).nMisses = nGt - nHits;
    results(rr).nFalsePos = nFalsePos;
    matched = matchIou(matchDet > 0);
    results(rr).meanIouMatched = mean(matched);
    results(rr).perCell = perCell;

    fprintf(['dilate=%2d px: %2d detected | %d/%d hand cells hit | %d false pos | ' ...
        'mean IoU (matched) = %.3f\n'], rad, nDet, nHits, nGt, nFalsePos, results(rr).meanIouMatched);

    % Overlay: green = auto (dilated), magenta = hand-drawn
    try
        fig = figure('Visible', 'off', 'Position', [50, 50, 1600, 500]);
        imagesc(img); colormap(gray); axis image; hold on;
        for gg = 1:nGt
            b = bwboundaries(gtMasks{gg}, 'noholes');
            for kk = 1:numel(b), plot(b{kk}(:,2), b{kk}(:,1), 'm-', 'LineWidth', 1.5); end
            text(gtCentroids(gg,2), gtCentroids(gg,1), sprintf('G%d', gg), 'Color', 'm', ...
                'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        end
        for dd = 1:nDet
            b = bwboundaries(detMasks{dd}, 'noholes');
            for kk = 1:numel(b), plot(b{kk}(:,2), b{kk}(:,1), 'g-', 'LineWidth', 1.0); end
        end
        title(sprintf('dilate=%d px: %d/%d hit, %d FP (magenta=hand, green=auto)', ...
            rad, nHits, nGt, nFalsePos));
        saveas(fig, fullfile(outDir, sprintf('overlay_dilate%02d.png', rad)));
        close(fig);
    catch qcErr
        warning('validate: overlay failed: %s', qcErr.message);
    end
end

% ---- Summary ------------------------------------------------------------------
save(fullfile(outDir, 'validation_results.mat'), 'results', 'matPath', 'dilateRadii');
fid = fopen(fullfile(outDir, 'validation_summary.txt'), 'w');
fprintf(fid, 'Ground truth: %d cells | %s\nHit = IoU > %.2f or centroid < %d px\n\n', ...
    nGt, matPath, hitIouThresh, hitCentroidPx);
fprintf(fid, 'dilate  nDet  hits  misses  falsePos  meanIoU(matched)\n');
for rr = 1:numel(results)
    fprintf(fid, '%5d  %4d  %4d  %6d  %8d  %16.3f\n', results(rr).dilateRadiusPx, ...
        results(rr).nDetected, results(rr).nHits, results(rr).nMisses, ...
        results(rr).nFalsePos, results(rr).meanIouMatched);
end
fclose(fid);
fprintf('Saved: %s\n', fullfile(outDir, 'validation_summary.txt'));
end
