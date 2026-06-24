%% VoltImg_cellPool_snrDprime_thresholds
% Pool slice-test analyses at the *cell* level (one saved voltImgTest_Analysis
% per cell). Each file already contains within-cell SNR vs |ΔV|, d', and
% interpolated detection thresholds (same logic as
% VoltImg_slice_test_analysis_MCfineROI_TrialSpecROIandNeuropil_042026.m).
%
% This script:
%  1) Loads many .mat files (or uses structs already in workspace).
%  2) Builds per-cell (dVm, y) curves for SNR or d' (thin line per cell).
%  3) Overlays pooled linear fit ± SEM (standard error of the fitted mean) on a
%     common |ΔV| grid.
%  4) Strip / swarm plots of per-cell thresholds with median + IQR, optional
%     split by cohort group label.
%  5) Optional inclusion rule: keep cells whose reference-stimulus SNR ≥ X
%     (reports exclusion fraction). Raw thresholds are primary; inclusion only
%     filters which dots enter the summary.
%  6) Optional save: poolCfg.saveOutputDir nonempty writes cellPoolSummary.mat
%     and, if poolCfg.saveFigures, raster/vector exports for screen (and publish)
%     figure numbers.
%
% Y-axis modes:
%   'snr' — snrTrialAverage.snrCalib_meanAcrossPulses (trial-mean peak / prestim σ).
%   'peak_filt_dff' — mean(peakMeanFiltDff,2) without σ (expression confounded).
%   'peak_over_sigma' — same as 'snr' for the curve; documented alias.
%
% Thresholds are always the stored interpolated |ΔV| (mV) from each cell's own
% SNR–|ΔV| or d'–|ΔV| curve. When minReferenceSnr is set, curves/threshold
% summaries use included cells only; the script also prints the same
% medians for *all* loaded cells so exclusion is transparent.
%
% Figures: primary set uses figureOffset+1..4 (light theme for non-dark UI).
% Optional publish set: makePublishFigures + publishFigureOffset (default
% figureOffset+100), 14 pt axis labels, square plot regions (pbaspect).
% poolCfg.plotAbsDVmMax_mV: SNR/d' curves and sensitivity use |ΔV_m| < this only
% (NaN = all points). Threshold swarm plots still use saved thresholds.
% Sensitivity (figureOffset+4): each point is one pulse × condition when the loaded
% analysis includes snrTrialAverage.peakMeanFiltDff and dVmMeasMeanTrace (from pulse
% timing and imagingFreq / ephys Fs in the slice pipeline); otherwise one point per condition.

%% ---------------- user configuration --------------------------------
clearvars -except voltImgTest_Analysis  % keep one struct in base if debugging

% Either set matFilePaths, or leave empty and assign poolCfg.analysisStructs
% as a cell array of voltImgTest_Analysis structs (same order as groupLabels).
% Each matFilePaths entry may be either:
%   - a saved analysis .mat file path, or
%   - a directory path (all *.mat files inside are pooled).
matFilePaths = {
    '/Volumes/phoenixinthesky/Masato/Voltage Imaging Data_Phoenix/voltImgTest/Analysis Results/ASAP7y_selected_for_SNR_1exclusion';
    };

% Parallel optional group / condition label per matFilePaths entry.
% If an entry is a directory, its label is applied to all .mat files expanded
% from that directory.
% Example: {'ASAP7','Control'}
% Empty string '' => single cohort "all".
groupLabels = {
    % 'ASAP7y';
    };

poolCfg = struct();
poolCfg.analysisStructs = {};  % if nonempty, overrides matFilePaths
poolCfg.yMetricSnrPanel = 'snr';   % 'snr' | 'peak_filt_dff'
poolCfg.yMetricDprimePanel = 'dprime_mf'; % only 'dprime_mf' supported (saved calib)
poolCfg.nGridInterp = 64;          % common |ΔV| samples for mean ± SEM
poolCfg.figureOffset = 200;       % fig numbers: offset+1, offset+2, ...
% Second figure set for publication (light theme, larger fonts, square axes).
poolCfg.makePublishFigures = true;
poolCfg.publishFigureOffset = [];  % [] => figureOffset + 100
% Inclusion: NaN => no filter. Else require reference SNR >= this (finite refs only).
poolCfg.minReferenceSnr = NaN;
% Reference SNR default = SNR at the calibration amplitude with largest |ΔV|.
% Optional: set referenceCondIndex (1..nConds) to pin the reference to one step.
poolCfg.referenceCondIndex = [];   % e.g. 5 = 5th row in each cell's saved dVmCalib vector
poolCfg.dPrimeTargetOverride = []; % [] => read sdtSingleTrial.dPrimeTarget from each file
poolCfg.snrThresholdListOverride = []; % [] => read snrCfg.snrThresholdList from first valid file
% Curves / sensitivity: keep calibration points with |ΔV_m| strictly below this (mV).
% NaN => use all amplitudes. Threshold strip plots still use full saved thresholds.
poolCfg.plotAbsDVmMax_mV = 75;
% Optional disk output (empty saveOutputDir => no files written).
poolCfg.saveOutputDir = '/Volumes/phoenixinthesky/Masato/Voltage Imaging Data_Phoenix/voltImgTest/Analysis Results/Pooled Analysis Result';   % e.g. '/path/to/cellPool_run1'
poolCfg.saveFileStem = 'cellPool_snrDprime';  % <stem>.mat, <stem>_screen_*.png, ...
poolCfg.saveFigures = true;
poolCfg.saveFigureFormats = {'png'};  % e.g. {'png','pdf'}; pdf uses print -dpdf

%% ---------------- load structs --------------------------------------
if isempty(poolCfg.analysisStructs)
    [matFilePaths, groupLabels] = expandMatFileInputs(matFilePaths, groupLabels);
    nF = numel(matFilePaths);
    if nF == 0
        error(['Set matFilePaths to a cell array of saved analysis .mat files ', ...
            '(or directories containing .mat files), ', ...
            'or set poolCfg.analysisStructs to {V1, V2, ...}.']);
    end
    analysisStructsRaw = cell(nF, 1);
    keepValid = false(nF, 1);
    for fi = 1:nF
        S = load(matFilePaths{fi}, 'voltImgTest_Analysis');
        if ~isfield(S, 'voltImgTest_Analysis')
            warning('Skipping file missing voltImgTest_Analysis: %s', matFilePaths{fi});
            continue;
        end
        V = S.voltImgTest_Analysis;
        if ~isfield(V, 'snrTrialAverage') || ~isfield(V, 'sdtSingleTrial')
            [V, didReconstruct] = reconstructLegacySnrDprime(V);
            if didReconstruct
                warning(['Reconstructed snrTrialAverage/sdtSingleTrial from legacy ', ...
                    'fields: %s'], matFilePaths{fi});
            end
        end
        if ~isfield(V, 'snrTrialAverage') || ~isfield(V, 'sdtSingleTrial')
            warning(['Skipping incompatible file after fallback attempt ', ...
                '(missing snrTrialAverage or sdtSingleTrial): %s'], matFilePaths{fi});
            continue;
        end
        analysisStructsRaw{fi} = V;
        keepValid(fi) = true;
    end
    if ~any(keepValid)
        error('No valid analysis files remained after loading/validation.');
    end
    poolCfg.analysisStructs = analysisStructsRaw(keepValid);
    matFilePaths = matFilePaths(keepValid);
    groupLabels = groupLabels(keepValid);
    nF = numel(poolCfg.analysisStructs);
    fprintf('Loaded %d valid analysis files (skipped %d invalid).\n', nF, sum(~keepValid));
else
    nF = numel(poolCfg.analysisStructs);
    matFilePaths = repmat({''}, nF, 1);
end

if isempty(groupLabels)
    groupLabels = repmat({''}, nF, 1);
elseif numel(groupLabels) ~= nF
    error('groupLabels must be empty or same length as number of analyses (%d).', nF);
end

%% ---------------- extract per-cell tables ----------------------------
cellRows = repmat(voltImgStructTemplate(), nF, 1);
for fi = 1:nF
    V = poolCfg.analysisStructs{fi};
    cellRows(fi) = extractOneCell(V, matFilePaths{fi});
    cellRows(fi).groupLabel = groupLabels{fi};
end

snrThr0 = poolCfg.snrThresholdListOverride;
if isempty(snrThr0)
    snrThr0 = defaultSnrThresholdList(poolCfg.analysisStructs{1});
end
dpTgt0 = poolCfg.dPrimeTargetOverride;
if isempty(dpTgt0)
    dpTgt0 = defaultDprimeTarget(poolCfg.analysisStructs{1});
end

refSnr = computeReferenceSnrPerCell(cellRows, poolCfg);
for fi = 1:nF
    cellRows(fi).snrReference = refSnr(fi);
end

nTotal = nF;
if ~isnan(poolCfg.minReferenceSnr)
    keep = isfinite(refSnr) & refSnr >= poolCfg.minReferenceSnr;
else
    keep = true(nF, 1);
end
exclFrac = sum(~keep) / max(1, nTotal);
cellRowsUse = cellRows(keep);
nUse = numel(cellRowsUse);
if nUse < 1
    error('No cells passed inclusion criterion (minReferenceSnr=%g).', poolCfg.minReferenceSnr);
end

%% ---------------- y vectors for SNR-style panel ---------------------
[yCellSnr, yLabelSnr] = buildYCalibPerCell(cellRowsUse, poolCfg.yMetricSnrPanel);
[dVmCell, yCellDp] = packDvmYCells(cellRowsUse);
[dVmCell, yCellSnr, yCellDp] = filterCalibCellsForPlot(dVmCell, yCellSnr, yCellDp, poolCfg.plotAbsDVmMax_mV);

%% ---------------- threshold collection -------------------------------
threshPack = collectThresholds(cellRowsUse, snrThr0, dpTgt0);
snrThrForFig = 2;
fldSnrFig = snrFieldTag(snrThrForFig);
thrValsSnrFig = vectorSnrThreshVm(cellRowsUse, fldSnrFig);
thrValsDpFig = [cellRowsUse.dVmAtDprime];

%% ---------------- figures: population trend + threshold inset --------
f1 = poolCfg.figureOffset + 1;
styleScreen = defaultFigureStyleOpts();
snrPlotPack = plotMetricWithThresholdInset( ...
    f1, dVmCell, yCellSnr, yLabelSnr, snrThrForFig, thrValsSnrFig, ...
    sprintf('SNR vs |\\Delta V_m| (n=%d)', nUse, nTotal), ...
    sprintf('|\\Delta V_m| at SNR=%g', snrThrForFig), poolCfg.nGridInterp, styleScreen);
xg = snrPlotPack.gridX;
muY = snrPlotPack.gridMean;
seY = snrPlotPack.gridSem;

f2 = poolCfg.figureOffset + 2;
yLabDp = 'd'' (matched filter, mean pulses)';
if ~strcmp(poolCfg.yMetricDprimePanel, 'dprime_mf')
    warning('Only dprime_mf supported; using sdtSingleTrial.dPrimeCalib_matchedFilter_meanPulses.');
end
dpPlotPack = plotMetricWithThresholdInset( ...
    f2, dVmCell, yCellDp, yLabDp, dpTgt0, thrValsDpFig, ...
    sprintf('d'' vs |\\Delta V_m| (n=%d)', nUse), ...
    sprintf('|\\Delta V_m| at d''=%g', dpTgt0), poolCfg.nGridInterp, styleScreen);
xg2 = dpPlotPack.gridX;
muDp = dpPlotPack.gridMean;
seDp = dpPlotPack.gridSem;

%% ---------------- per-cell sensitivity: peak mean dF/F vs peak mean Vm ----
f4 = poolCfg.figureOffset + 4;
sensitivityPack = plotPerCellSensitivity(cellRowsUse, f4, styleScreen, poolCfg.plotAbsDVmMax_mV);

%% ---------------- strip / swarm plots -------------------------------
f3 = poolCfg.figureOffset + 3;
figure(f3); clf;
set(gcf, 'Color', 'w', 'Position', [160, 90, 640, 420]);
tiledlayout('flow', 'Padding', 'compact', 'TileSpacing', 'compact');
uniqG = unique(groupLabelsFromRows(cellRowsUse), 'stable');
if isscalar(uniqG) && strcmp(uniqG{1}, '')
    uniqG = {'all'};
    for fi = 1:nUse
        cellRowsUse(fi).groupLabel = 'all';
    end
end

gList = groupLabelsFromRows(cellRowsUse);

for ti = 1:numel(snrThr0)
    thr = snrThr0(ti);
    fld = snrFieldTag(thr);
    vals = vectorSnrThreshVm(cellRowsUse, fld);
    nexttile;
    swarmStripByGroup(vals, gList, uniqG, sprintf('|\\Delta V_m| at SNR=%g', thr), ...
        styleScreen.fontSize, styleScreen.showGrid);
end

nexttile;
valsDp = [cellRowsUse.dVmAtDprime];
swarmStripByGroup(valsDp, gList, uniqG, sprintf('|\\Delta V_m| at d''=%g (MF)', dpTgt0), ...
    styleScreen.fontSize, styleScreen.showGrid);
sgtitle(sprintf(['|\\Delta V_m| thresholds (median bar + IQR whisker). ', ...
    'Included n=%d / %d (excluded %.0f%%)'], nUse, nTotal, 100 * exclFrac), ...
    'FontWeight', 'normal', 'FontSize', styleScreen.fontSize + 1, 'Interpreter', 'tex');
applyLightThemeToFigure(figure(f3), styleScreen);

if poolCfg.makePublishFigures
    pubOff = poolCfg.publishFigureOffset;
    if isempty(pubOff)
        pubOff = poolCfg.figureOffset + 100;
    end
    stylePub = publishFigureStyleOpts();
    plotMetricWithThresholdInset( ...
        pubOff + 1, dVmCell, yCellSnr, yLabelSnr, snrThrForFig, thrValsSnrFig, ...
        sprintf('SNR vs |\\Delta V_m| (n=%d)', nUse, nTotal), ...
        sprintf('|\\Delta V_m| at SNR=%g', snrThrForFig), poolCfg.nGridInterp, stylePub);
    plotMetricWithThresholdInset( ...
        pubOff + 2, dVmCell, yCellDp, yLabDp, dpTgt0, thrValsDpFig, ...
        sprintf('d'' vs |\\Delta V_m| (n=%d)', nUse), ...
        sprintf('|\\Delta V_m| at d''=%g', dpTgt0), poolCfg.nGridInterp, stylePub);
    plotPerCellSensitivity(cellRowsUse, pubOff + 4, stylePub, poolCfg.plotAbsDVmMax_mV);
    figure(pubOff + 3); clf;
    set(gcf, 'Color', 'w', 'Position', [180, 70, 720, 480]);
    tiledlayout('flow', 'Padding', 'compact', 'TileSpacing', 'compact');
    for ti = 1:numel(snrThr0)
        thr = snrThr0(ti);
        fld = snrFieldTag(thr);
        vals = vectorSnrThreshVm(cellRowsUse, fld);
        nexttile;
        swarmStripByGroup(vals, gList, uniqG, sprintf('|\\Delta V_m| at SNR=%g', thr), ...
            stylePub.fontSize, stylePub.showGrid);
        if stylePub.squareAxes
            axis(gca, 'tight');
            pbaspect(gca, [1, 1, 1]);
        end
    end
    nexttile;
    swarmStripByGroup(valsDp, gList, uniqG, sprintf('|\\Delta V_m| at d''=%g (MF)', dpTgt0), ...
        stylePub.fontSize, stylePub.showGrid);
    if stylePub.squareAxes
        axis(gca, 'tight');
        pbaspect(gca, [1, 1, 1]);
    end
    sgtitle(sprintf(['|\\Delta V_m| thresholds (median bar + IQR whisker). ', ...
        'Included n=%d / %d (excluded %.0f%%)'], nUse, nTotal, 100 * exclFrac), ...
        'FontWeight', 'normal', 'FontSize', stylePub.fontSize + 1, 'Interpreter', 'tex');
    applyLightThemeToFigure(gcf, stylePub);
    poolCfg.publishFigureOffset = pubOff;
end

%% ---------------- console summary (paper-ready sentence fragments) -----
disp(' ');
if isnan(poolCfg.minReferenceSnr)
    disp('=== Cell-level pool: per-cell |ΔV| thresholds (median / IQR) ===');
    printThresholdSummary(cellRowsUse, snrThr0, dpTgt0);
else
    disp('=== Thresholds: ALL loaded cells (before reference-SNR gate) ===');
    printThresholdSummary(cellRows, snrThr0, dpTgt0);
    disp('=== Thresholds: INCLUDED cells (curves / strip plots use this set) ===');
    printThresholdSummary(cellRowsUse, snrThr0, dpTgt0);
    disp(sprintf(['  Gate: reference SNR >= %g  =>  excluded %.1f%% (%d/%d cells).'], ...
        poolCfg.minReferenceSnr, 100 * exclFrac, nTotal - nUse, nTotal));
end

%% ---------------- export struct --------------------------------------
cellPoolSummary = struct();
cellPoolSummary.poolCfg = poolCfg;
cellPoolSummary.snrThresholdList = snrThr0;
cellPoolSummary.dPrimeTarget = dpTgt0;
cellPoolSummary.nCellsTotal = nTotal;
cellPoolSummary.nCellsIncluded = nUse;
cellPoolSummary.exclusionFraction = exclFrac;
cellPoolSummary.cellRowsAll = cellRows;
cellPoolSummary.cellRowsIncluded = cellRowsUse;
cellPoolSummary.curveGrid.x = xg;
cellPoolSummary.curveGrid.snr_mean = muY;
cellPoolSummary.curveGrid.snr_sem = seY;
cellPoolSummary.curveGrid.dPrime_mean = muDp;
cellPoolSummary.curveGrid.dPrime_sem = seDp;
cellPoolSummary.thresholds = threshPack;
cellPoolSummary.sensitivityVmDff = sensitivityPack;
cellPoolSummary.yLabelSnrPanel = yLabelSnr;
cellPoolSummary.refSnrPerCell = refSnr;
cellPoolSummary.includedLogical = keep;
if poolCfg.makePublishFigures && ~isempty(poolCfg.publishFigureOffset)
    cellPoolSummary.publishFigureNumbers = poolCfg.publishFigureOffset + (1:4);
else
    cellPoolSummary.publishFigureNumbers = [];
end

%% ---------------- save results (optional) -------------------------------
saveCellPoolOutputs(poolCfg, cellPoolSummary);

%% ========================================================================
function saveCellPoolOutputs(poolCfg, cellPoolSummary)
    if ~isfield(poolCfg, 'saveOutputDir') || isempty(poolCfg.saveOutputDir)
        return;
    end
    outDir = char(strtrim(string(poolCfg.saveOutputDir)));
    if isempty(outDir)
        return;
    end
    if ~isfolder(outDir)
        mkdir(outDir);
    end
    stem = 'cellPool_snrDprime';
    if isfield(poolCfg, 'saveFileStem') && ~isempty(poolCfg.saveFileStem)
        stem = char(strtrim(string(poolCfg.saveFileStem)));
    end
    matPath = fullfile(outDir, [stem '.mat']);
    save(matPath, 'cellPoolSummary', '-v7.3');
    fprintf('Saved cell pool summary: %s\n', matPath);

    doFig = isfield(poolCfg, 'saveFigures') && poolCfg.saveFigures;
    if ~doFig
        return;
    end
    fmtList = {'png'};
    if isfield(poolCfg, 'saveFigureFormats') && ~isempty(poolCfg.saveFigureFormats)
        fmtList = poolCfg.saveFigureFormats;
        if ischar(fmtList)
            fmtList = cellstr(fmtList);
        elseif isstring(fmtList)
            fmtList = cellstr(fmtList);
        end
    end
    figLabels = {'snr_curve', 'dprime_curve', 'threshold_swarm', 'sensitivity_vm_dff'};
    baseScreen = poolCfg.figureOffset + (1:4);
    for ii = 1:4
        exportOneFigureByNumber(baseScreen(ii), fullfile(outDir, ...
            sprintf('%s_screen_%s', stem, figLabels{ii})), fmtList);
    end
    if isfield(poolCfg, 'makePublishFigures') && poolCfg.makePublishFigures ...
            && isfield(poolCfg, 'publishFigureOffset') && ~isempty(poolCfg.publishFigureOffset)
        pubOff = poolCfg.publishFigureOffset;
        for ii = 1:4
            exportOneFigureByNumber(pubOff + ii, fullfile(outDir, ...
                sprintf('%s_publish_%s', stem, figLabels{ii})), fmtList);
        end
    end
end

function exportOneFigureByNumber(figNum, pathWithoutExt, formats)
    fh = findall(0, 'Type', 'figure', 'Number', figNum);
    if isempty(fh)
        fprintf('Save figures: no figure with Number=%g, skipped %s\n', figNum, pathWithoutExt);
        return;
    end
    fh = fh(1);
    for k = 1:numel(formats)
        fmt = lower(strtrim(char(string(formats{k}))));
        outFile = [pathWithoutExt '.' fmt];
        if strcmp(fmt, 'pdf')
            print(fh, outFile, '-dpdf', '-painters');
        elseif strcmp(fmt, 'png')
            print(fh, outFile, '-dpng', '-r300');
        elseif strcmp(fmt, 'eps')
            print(fh, outFile, '-depsc', '-painters');
        elseif strcmp(fmt, 'svg')
            print(fh, outFile, '-dsvg');
        else
            warning('VoltImg_cellPool_snrDprime_thresholds:UnknownFigFormat', ...
                'Unknown saveFigureFormats entry "%s", skipping.', fmt);
        end
    end
    fprintf('Wrote figure %g -> %s.*\n', figNum, pathWithoutExt);
end

%% ========================================================================
function T = voltImgStructTemplate()
    T = struct('dVmCalib', [], 'snrCalib', [], 'peakCalib', [], 'dPrimeCalibMF', [], ...
        'dVmPerPulseMeas', [], 'peakFiltDffPerPulse', [], ...
        'sigmaDffBaseline', NaN, 'dVmAtSnr', struct(), 'dVmAtDprime', NaN, ...
        'snrReference', NaN, 'mouseID', '', 'sourcePath', '', 'groupLabel', '');
end

function [expandedPaths, expandedLabels] = expandMatFileInputs(matFilePaths, groupLabels)
    if isempty(matFilePaths)
        expandedPaths = {};
        expandedLabels = {};
        return;
    end

    nInput = numel(matFilePaths);
    if isempty(groupLabels)
        groupLabels = repmat({''}, nInput, 1);
    elseif numel(groupLabels) ~= nInput
        error(['groupLabels must be empty or same length as matFilePaths ', ...
            '(%d).'], nInput);
    end

    expandedPaths = {};
    expandedLabels = {};
    for ii = 1:nInput
        p = matFilePaths{ii};
        g = groupLabels{ii};

        if isfolder(p)
            D = dir(fullfile(p, '*.mat'));
            if isempty(D)
                warning('No .mat files found in directory: %s', p);
                continue;
            end
            for kk = 1:numel(D)
                expandedPaths{end+1,1} = fullfile(D(kk).folder, D(kk).name); %#ok<AGROW>
                expandedLabels{end+1,1} = g; %#ok<AGROW>
            end
        elseif isfile(p)
            expandedPaths{end+1,1} = p; %#ok<AGROW>
            expandedLabels{end+1,1} = g; %#ok<AGROW>
        else
            error('matFilePaths entry is not a valid file or directory: %s', p);
        end
    end
end

function [V, didReconstruct] = reconstructLegacySnrDprime(V)
    didReconstruct = false;
    if isfield(V, 'snrTrialAverage') && isfield(V, 'sdtSingleTrial')
        return;
    end
    req = {'dvToTest', 'filtdffAllConds'};
    if ~all(cellfun(@(f) isfield(V, f), req))
        return;
    end

    dVmCalib = V.dvToTest(:);
    nCond = numel(dVmCalib);
    if ~iscell(V.filtdffAllConds) || numel(V.filtdffAllConds) ~= nCond
        return;
    end

    peakByCond = cell(nCond, 1);
    baselineSigmaByCond = nan(nCond, 1);
    baselinePool = [];
    for ii = 1:nCond
        A = V.filtdffAllConds{ii};
        if isempty(A) || ndims(A) ~= 2
            continue;
        end
        nFrames = size(A, 1);
        nBase = max(5, round(0.2 * nFrames));
        idxBase = 1:min(nBase, nFrames);
        baseVals = A(idxBase, :);
        baselinePool = [baselinePool; baseVals(:)]; %#ok<AGROW>
        baselineSigmaByCond(ii) = std(baseVals(:), 'omitnan');
        peakByCond{ii} = max(A, [], 1, 'omitnan');
    end
    if isempty(baselinePool)
        return;
    end

    peakMean = nan(nCond, 1);
    for ii = 1:nCond
        if ~isempty(peakByCond{ii})
            peakMean(ii) = mean(peakByCond{ii}, 'omitnan');
        end
    end
    sigmaDffBaseline = median(baselineSigmaByCond, 'omitnan');
    if ~isfinite(sigmaDffBaseline) || sigmaDffBaseline <= 0
        sigmaDffBaseline = std(baselinePool, 'omitnan');
    end
    if ~isfinite(sigmaDffBaseline) || sigmaDffBaseline <= 0
        return;
    end
    snrCalib = peakMean ./ sigmaDffBaseline;

    snrThrList = [2, 3];
    if isfield(V, 'snrCfg') && isfield(V.snrCfg, 'snrThresholdList')
        snrThrList = V.snrCfg.snrThresholdList(:).';
    end
    dVmAtSnrThreshold = struct();
    for kk = 1:numel(snrThrList)
        fld = snrFieldTag(snrThrList(kk));
        dVmAtSnrThreshold.(fld) = thresholdCrossingVm(dVmCalib, snrCalib, snrThrList(kk));
    end

    dPrimeTarget = 1;
    if isfield(V, 'sdtCfg') && isfield(V.sdtCfg, 'dPrimeTarget')
        dPrimeTarget = V.sdtCfg.dPrimeTarget;
    end
    muBase = mean(baselinePool, 'omitnan');
    varBase = var(baselinePool, 'omitnan');
    dPrimeCalib = nan(nCond, 1);
    for ii = 1:nCond
        y = peakByCond{ii};
        if isempty(y)
            continue;
        end
        muY = mean(y, 'omitnan');
        varY = var(y, 'omitnan');
        den = sqrt(0.5 * (varY + varBase));
        if den > 0 && isfinite(den)
            dPrimeCalib(ii) = (muY - muBase) / den;
        end
    end
    dVmAtDprimeTarget = thresholdCrossingVm(dVmCalib, dPrimeCalib, dPrimeTarget);

    V.snrTrialAverage = struct();
    V.snrTrialAverage.sigmaDffBaseline = sigmaDffBaseline;
    V.snrTrialAverage.peakMeanFiltDff = peakMean;
    V.snrTrialAverage.dVmNominal = dVmCalib;
    V.snrTrialAverage.dVmCalib = dVmCalib;
    V.snrTrialAverage.snrCalib_meanAcrossPulses = snrCalib;
    V.snrTrialAverage.dVmAtSnrThreshold = dVmAtSnrThreshold;
    V.snrTrialAverage.snrThresholdList = snrThrList;
    V.snrTrialAverage.note = ['Reconstructed from legacy fields (dvToTest + ', ...
        'filtdffAllConds).'];

    V.sdtSingleTrial = struct();
    V.sdtSingleTrial.dPrimeCalib_matchedFilter_meanPulses = dPrimeCalib;
    V.sdtSingleTrial.dVmCalib = dVmCalib;
    V.sdtSingleTrial.dVmNominal = dVmCalib;
    V.sdtSingleTrial.dVmAtDprimeTarget = dVmAtDprimeTarget;
    V.sdtSingleTrial.dPrimeTarget = dPrimeTarget;
    V.sdtSingleTrial.note = ['Approximate d'' reconstructed from peak response vs ', ...
        'baseline distributions.'];
    didReconstruct = true;
end

function xCross = thresholdCrossingVm(x, y, yTarget)
    xCross = NaN;
    x = x(:);
    y = y(:);
    good = isfinite(x) & isfinite(y);
    if nnz(good) < 2
        return;
    end
    x = x(good);
    y = y(good);
    [x, ord] = sort(x, 'ascend');
    y = y(ord);
    for ii = 1:(numel(x) - 1)
        y1 = y(ii);
        y2 = y(ii + 1);
        if y1 == y2
            if y1 >= yTarget
                xCross = x(ii);
                return;
            end
            continue;
        end
        if (y1 - yTarget) * (y2 - yTarget) <= 0
            t = (yTarget - y1) / (y2 - y1);
            xCross = x(ii) + t * (x(ii + 1) - x(ii));
            return;
        end
    end
end

function R = extractOneCell(V, srcPath)
    R = voltImgStructTemplate();
    if isfield(V, 'mouseID')
        R.mouseID = V.mouseID;
    end
    R.sourcePath = srcPath;
    if ~isfield(V, 'snrTrialAverage') || ~isfield(V, 'sdtSingleTrial')
        error('voltImgTest_Analysis missing snrTrialAverage or sdtSingleTrial (%s).', R.mouseID);
    end
    sn = V.snrTrialAverage;
    sd = V.sdtSingleTrial;
    R.dVmCalib = sn.dVmCalib(:);
    R.snrCalib = sn.snrCalib_meanAcrossPulses(:);
    R.sigmaDffBaseline = sn.sigmaDffBaseline;
    if isfield(sn, 'dVmAtSnrThreshold')
        R.dVmAtSnr = sn.dVmAtSnrThreshold;
    else
        R.dVmAtSnr = struct();
    end
    R.dPrimeCalibMF = sd.dPrimeCalib_matchedFilter_meanPulses(:);
    if isfield(sd, 'dVmAtDprimeTarget')
        R.dVmAtDprime = sd.dVmAtDprimeTarget;
    else
        R.dVmAtDprime = NaN;
    end
    n = numel(R.dVmCalib);
    if isfield(sn, 'peakMeanFiltDff')
        R.peakCalib = mean(sn.peakMeanFiltDff, 2, 'omitnan');
    else
        R.peakCalib = nan(n, 1);
    end
    % Per-pulse pairs for sensitivity scatter (nConds×nPulses); same windows as slice analysis
    % (pulse timing × imagingFreq / ephys Fs → peakMeanFiltDff, dVmMeasMeanTrace).
    R.dVmPerPulseMeas = [];
    R.peakFiltDffPerPulse = [];
    if isfield(sn, 'peakMeanFiltDff') && isfield(sn, 'dVmMeasMeanTrace')
        Pk = sn.peakMeanFiltDff;
        dV = sn.dVmMeasMeanTrace;
        if isnumeric(Pk) && isnumeric(dV) && isequal(size(Pk), size(dV)) && ~isempty(Pk)
            R.peakFiltDffPerPulse = Pk;
            R.dVmPerPulseMeas = dV;
        end
    end
    if numel(R.snrCalib) ~= n || numel(R.dPrimeCalibMF) ~= n || numel(R.peakCalib) ~= n
        error('Calibration vector length mismatch for %s.', R.mouseID);
    end
end

function snrThr = defaultSnrThresholdList(V)
    snrThr = [2, 3];
    if isfield(V, 'snrCfg') && isfield(V.snrCfg, 'snrThresholdList')
        snrThr = V.snrCfg.snrThresholdList(:).';
    end
end

function dp = defaultDprimeTarget(V)
    dp = 1;
    if isfield(V, 'sdtSingleTrial') && isfield(V.sdtSingleTrial, 'dPrimeTarget')
        dp = V.sdtSingleTrial.dPrimeTarget;
    elseif isfield(V, 'snrCfg') && isfield(V.snrCfg, 'dPrimeTarget')
        dp = V.snrCfg.dPrimeTarget;
    end
end

function ref = computeReferenceSnrPerCell(cellRows, poolCfg)
    n = numel(cellRows);
    ref = nan(n, 1);
    for fi = 1:n
        dv = cellRows(fi).dVmCalib(:);
        sn = cellRows(fi).snrCalib(:);
        ok = ~isnan(dv) & ~isnan(sn);
        dv = dv(ok);
        sn = sn(ok);
        if isempty(sn)
            continue
        end
        if ~isempty(poolCfg.referenceCondIndex)
            ix = poolCfg.referenceCondIndex;
            if ix >= 1 && ix <= numel(sn)
                ref(fi) = sn(ix);
            end
        else
            [~, j] = max(abs(dv));
            ref(fi) = sn(j);
        end
    end
end

function [yCell, yLabel] = buildYCalibPerCell(rows, metric)
    n = numel(rows);
    yCell = cell(n, 1);
    for fi = 1:n
        switch lower(strrep(metric, ' ', '_'))
            case {'snr', 'peak_over_sigma'}
                yCell{fi} = rows(fi).snrCalib(:);
                yLabel = 'SNR (trial-mean peak \Delta F/F / \sigma_{pre})';
            case 'peak_filt_dff'
                yCell{fi} = rows(fi).peakCalib(:);
                yLabel = 'trial-mean peak filtered \Delta F/F (a.u.; not divided by prestim \sigma)';
            otherwise
                error('Unknown yMetric: %s', metric);
        end
    end
end

function plotThinCellLines(dVmCell, yCell, col, lw)
    for k = 1:numel(dVmCell)
        x = dVmCell{k}(:);
        y = yCell{k}(:);
        m = ~isnan(x) & ~isnan(y);
        x = x(m);
        y = y(m);
        if numel(x) < 2
            continue
        end
        [x, ord] = sort(x);
        y = y(ord);
        plot(x, y, '-', 'Color', col, 'LineWidth', lw);
    end
end

function [xg, mu, se] = meanSemOnCommonGrid(dVmCell, yCell, nGrid)
    % Pooled ordinary least squares on all finite (|ΔV_m|, y) points, evaluated on
    % xg; se is the standard error of the fitted mean response (not across-cell
    % interp SEM). Per-cell colored fits in plotMetricWithThresholdInset stay linear.
    allX = [];
    allY = [];
    for k = 1:numel(dVmCell)
        x = dVmCell{k}(:);
        y = yCell{k}(:);
        m = isfinite(x) & isfinite(y);
        x = x(m);
        y = y(m);
        allX = [allX; x]; %#ok<AGROW>
        allY = [allY; y]; %#ok<AGROW>
    end
    if isempty(allX)
        xg = linspace(0, 1, nGrid);
        mu = nan(size(xg));
        se = nan(size(xg));
        return
    end
    lo = min(allX);
    hi = max(allX);
    if hi <= lo
        hi = lo + eps;
    end
    xg = linspace(lo, hi, nGrid);
    nPts = numel(allX);
    if nPts < 2
        mu = nan(size(xg));
        se = nan(size(xg));
        return
    end
    p = polyfit(allX, allY, 1);
    mu = polyval(p, xg);
    yhat = polyval(p, allX);
    res = allY - yhat;
    df = max(nPts - 2, 1);
    s = sqrt(sum(res.^2) / df);
    xbar = mean(allX);
    Sxx = sum((allX - xbar).^2);
    if Sxx < eps
        se = (s ./ sqrt(max(nPts, 1))) * ones(size(xg));
    else
        se = s .* sqrt(1 ./ nPts + (xg - xbar).^2 ./ Sxx);
    end
end

function s = defaultFigureStyleOpts()
    s = struct('fontSize', 11, 'squareAxes', false, 'hideAxesToolbar', true, 'showGrid', true, ...
        'showThresholdInset', true, 'xAxisZeroMin', false);
end

function s = publishFigureStyleOpts()
    s = struct('fontSize', 14, 'squareAxes', true, 'hideAxesToolbar', true, 'showGrid', false, ...
        'showThresholdInset', false, 'xAxisZeroMin', true);
end

function applyLightThemeToFigure(fig, styleOpts)
    if nargin < 2 || isempty(styleOpts)
        styleOpts = defaultFigureStyleOpts();
    end
    if ~isgraphics(fig)
        return
    end
    set(fig, 'Color', 'w');
    axList = findall(fig, 'Type', 'axes');
    for k = 1:numel(axList)
        applyLightThemeToAxes(axList(k), styleOpts);
    end
end

function applyLightThemeToAxes(ax, styleOpts)
    if ~isgraphics(ax) || ~strcmp(ax.Type, 'axes')
        return
    end
    fs = styleOpts.fontSize;
    set(ax, 'Color', 'w', ...
        'XColor', [0.12 0.12 0.12], ...
        'YColor', [0.12 0.12 0.12], ...
        'ZColor', [0.12 0.12 0.12], ...
        'GridColor', [0.78 0.78 0.82], ...
        'MinorGridColor', [0.88 0.88 0.92], ...
        'GridAlpha', 0.9, ...
        'Box', 'on', ...
        'LineWidth', 1.5, ...
        'TickDir', 'in', ...
        'FontSize', fs);
    if isfield(styleOpts, 'hideAxesToolbar') && styleOpts.hideAxesToolbar
        try
            if isprop(ax, 'Toolbar') && ~isempty(ax.Toolbar) && isprop(ax.Toolbar, 'Visible')
                ax.Toolbar.Visible = 'off';
            end
        catch %#ok<CTCH>
        end
    end
    if isfield(styleOpts, 'showGrid') && ~styleOpts.showGrid
        grid(ax, 'off');
    else
        grid(ax, 'on');
    end
end

function plotPack = plotMetricWithThresholdInset(figNum, dVmCell, yCell, yLabel, yThresh, thrVals, mainTitle, insetTitle, nGrid, styleOpts)
    if nargin < 10 || isempty(styleOpts)
        styleOpts = defaultFigureStyleOpts();
    end
    fs = styleOpts.fontSize;
    fsInset = max(8, round(0.85 * fs));
    showInset = true;
    if isfield(styleOpts, 'showThresholdInset') && ~styleOpts.showThresholdInset
        showInset = false;
    end

    figure(figNum); clf;
    set(gcf, 'Color', 'w', 'Position', [100, 100, 680, 460]);

    if showInset
        axMain = axes('Position', [0.10, 0.13, 0.62, 0.80]); %#ok<LAXES>
    else
        axMain = axes('Position', [0.12, 0.13, 0.86, 0.80]); %#ok<LAXES>
    end
    hold(axMain, 'on');

    [xg, mu, se] = meanSemOnCommonGrid(dVmCell, yCell, nGrid);
    okp = isfinite(xg) & isfinite(mu) & isfinite(se);
    xgv = xg(okp);
    muv = mu(okp);
    sev = se(okp);
    colMean = [0.5 0.5 0.53];
    colSem = [0.74 0.74 0.78];
    if numel(xgv) >= 1
        plot(axMain, xgv, muv + sev, '--', 'Color', colSem, 'LineWidth', 1.05);
        plot(axMain, xgv, muv - sev, '--', 'Color', colSem, 'LineWidth', 1.05);
        plot(axMain, xgv, muv, '-', 'Color', colMean, 'LineWidth', 1.85);
    end

    n = numel(dVmCell);
    cols = lines(max(n, 7));
    for ii = 1:n
        x = dVmCell{ii}(:);
        y = yCell{ii}(:);
        ok = isfinite(x) & isfinite(y);
        x = x(ok);
        y = y(ok);
        if isempty(x)
            continue
        end
        cLine = max(0, min(1, cols(ii, :)));
        scatter(axMain, x, y, 44, cLine, 'filled', ...
            'MarkerFaceAlpha', 1, 'MarkerEdgeColor', [0.06 0.06 0.06], 'LineWidth', 0.75);
        if numel(x) >= 2
            p = polyfit(x, y, 1);
            xx = linspace(min(x), max(x), max(120, 40 * numel(x)));
            yy = polyval(p, xx);
            plot(axMain, xx, yy, '-', 'Color', cLine, 'LineWidth', 2.8);
        end
    end
    yline(axMain, yThresh, ':', 'LineWidth', 1.2, 'Color', [0.2 0.2 0.2]);
    xlabel(axMain, '|\Delta V_m| calibration (mV)', 'FontSize', fs, 'Color', [0.05 0.05 0.05], 'Interpreter', 'tex');
    ylabel(axMain, yLabel, 'FontSize', fs, 'Color', [0.05 0.05 0.05], 'Interpreter', 'tex');
    title(axMain, mainTitle, 'FontWeight', 'normal', 'FontSize', fs, 'Color', [0.05 0.05 0.05], 'Interpreter', 'tex');
    hold(axMain, 'off');

    if showInset
        axInset = axes('Position', [0.76, 0.57, 0.20, 0.31]); %#ok<LAXES>
        hold(axInset, 'on');
        v = thrVals(:);
        v = v(isfinite(v));
        if isempty(v)
            text(axInset, 0.5, 0.5, 'No valid thresholds', 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontSize', fsInset, 'Color', [0.05 0.05 0.05]);
        else
            jitter = 0.10 * (rand(numel(v), 1) - 0.5);
            scatter(axInset, 1 + jitter, v, 28, [0.08 0.08 0.08], 'filled', 'MarkerFaceAlpha', 1, ...
                'MarkerEdgeColor', [0.02 0.02 0.02], 'LineWidth', 0.45);
            muV = mean(v, 'omitnan');
            seV = std(v, 0, 'omitnan') ./ sqrt(numel(v));
            errorbar(axInset, 1.17, muV, seV, 'k', 'LineWidth', 1.3, 'CapSize', 8, ...
                'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', [0.1 0.1 0.1], 'Color', [0.05 0.05 0.05]);
        end
        xlim(axInset, [0.75, 1.35]);
        set(axInset, 'XTick', []);
        ylabel(axInset, '|\Delta V_m| threshold (mV)', 'FontSize', fsInset, 'Color', [0.05 0.05 0.05], 'Interpreter', 'tex');
        title(axInset, sprintf('%s (mean \\pm SEM)', insetTitle), 'FontSize', fsInset, ...
            'FontWeight', 'normal', 'Color', [0.05 0.05 0.05], 'Interpreter', 'tex');
        hold(axInset, 'off');
    end

    if styleOpts.squareAxes
        axis(axMain, 'tight');
        pbaspect(axMain, [1, 1, 1]);
        if showInset
            axis(axInset, 'tight');
            pbaspect(axInset, [1, 1, 1]);
        end
    end
    if isfield(styleOpts, 'xAxisZeroMin') && styleOpts.xAxisZeroMin
        axis(axMain, 'tight');
        xl = xlim(axMain);
        xlim(axMain, [0, max(xl(2), 1e-6)]);
    end

    applyLightThemeToFigure(figure(figNum), styleOpts);

    plotPack = struct();
    plotPack.gridX = xg;
    plotPack.gridMean = mu;
    plotPack.gridSem = se;
end

function swarmStripByGroup(values, groupCell, uniqGroups, ttl, fontSize, showGrid)
    if nargin < 5 || isempty(fontSize)
        fontSize = 11;
    end
    if nargin < 6 || isempty(showGrid)
        showGrid = true;
    end
    axesH = gca;
    hold(axesH, 'on');
    ng = numel(uniqGroups);
    xv = 1:ng;
    cols = lines(max(ng, 3));
    for gi = 1:ng
        gname = uniqGroups{gi};
        idx = strcmp(groupCell, gname);
        v = values(idx);
        v = v(:);
        v = v(~isnan(v));
        if isempty(v)
            continue
        end
        nP = numel(v);
        jitter = zeros(nP, 1);
        cG = max(0, min(1, cols(gi, :)));
        scatter(xv(gi) + jitter, v, 40, cG, 'filled', ...
            'MarkerFaceAlpha', 1, 'MarkerEdgeColor', [0.06 0.06 0.06], 'LineWidth', 0.55);
    end
    set(axesH, 'XTick', xv, 'XTickLabel', uniqGroups);
    ylabel(axesH, '|\Delta V_m| (mV)', 'FontSize', fontSize, 'Color', [0.05 0.05 0.05], 'Interpreter', 'tex');
    title(axesH, ttl, 'FontSize', fontSize, 'Color', [0.05 0.05 0.05], 'Interpreter', 'tex');
    if showGrid
        grid(axesH, 'on');
    else
        grid(axesH, 'off');
    end
    hold(axesH, 'off');
    set(axesH, 'FontSize', fontSize, 'XColor', [0.12 0.12 0.12], 'YColor', [0.12 0.12 0.12], ...
        'GridColor', [0.78 0.78 0.82], 'Color', 'w', 'TickDir', 'in', 'Box', 'on', 'LineWidth', 1.5);
end

function fld = snrFieldTag(thr)
    fld = ['snr', strrep(num2str(thr, '%g'), '.', 'p')];
end

function pack = collectThresholds(rows, snrThr, dpTgt)
    pack = struct();
    pack.dPrimeTarget = dpTgt;
    pack.dVmAtDprime = [rows.dVmAtDprime];
    for ti = 1:numel(snrThr)
        fld = snrFieldTag(snrThr(ti));
        pack.(fld) = vectorSnrThreshVm(rows, fld);
    end
end

function printThresholdSummary(rows, snrThr, dpTgt)
    n = numel(rows);
    fprintf('  n cells = %d\n', n);
    for ti = 1:numel(snrThr)
        fld = snrFieldTag(snrThr(ti));
        if n < 1 || ~isfield(rows(1).dVmAtSnr, fld)
            continue
        end
        v = vectorSnrThreshVm(rows, fld);
        v = v(~isnan(v));
        printOneMetric(sprintf('SNR=%g', snrThr(ti)), v);
    end
    vdp = [rows.dVmAtDprime];
    vdp = vdp(~isnan(vdp));
    printOneMetric(sprintf('d''=%g (MF)', dpTgt), vdp);
end

function printOneMetric(nameStr, v)
    if isempty(v)
        fprintf('  %s: no finite values\n', nameStr);
        return
    end
    md = median(v);
    q1 = prctile(v, 25);
    q3 = prctile(v, 75);
    fprintf(['  %s: median |ΔV| = %.2f mV (IQR %.2f–%.2f) over n=%d cells\n'], ...
        nameStr, md, q1, q3, numel(v));
end

function v = vectorSnrThreshVm(rows, fld)
    n = numel(rows);
    v = nan(n, 1);
    for ii = 1:n
        if isfield(rows(ii).dVmAtSnr, fld)
            v(ii) = rows(ii).dVmAtSnr.(fld);
        end
    end
end

function [dVmCell, yCellDp] = packDvmYCells(rows)
    n = numel(rows);
    dVmCell = cell(n, 1);
    yCellDp = cell(n, 1);
    for ii = 1:n
        dVmCell{ii} = rows(ii).dVmCalib(:);
        yCellDp{ii} = rows(ii).dPrimeCalibMF(:);
    end
end

function [dVmCell, ySnr, yDp] = filterCalibCellsForPlot(dVmCell, ySnr, yDp, maxAbs_mV)
    if nargin < 4 || isempty(maxAbs_mV) || ~isfinite(maxAbs_mV)
        return
    end
    n = numel(dVmCell);
    for k = 1:n
        dv = dVmCell{k}(:);
        ys = ySnr{k}(:);
        yd = yDp{k}(:);
        nn = min([numel(dv), numel(ys), numel(yd)]);
        if nn < 1
            dVmCell{k} = [];
            ySnr{k} = [];
            yDp{k} = [];
            continue
        end
        dv = dv(1:nn);
        ys = ys(1:nn);
        yd = yd(1:nn);
        m = isfinite(dv) & isfinite(ys) & isfinite(yd) & (abs(dv) < maxAbs_mV);
        dVmCell{k} = dv(m);
        ySnr{k} = ys(m);
        yDp{k} = yd(m);
    end
end

function gList = groupLabelsFromRows(rows)
    n = numel(rows);
    gList = cell(n, 1);
    for ii = 1:n
        gList{ii} = rows(ii).groupLabel;
    end
end

function sensitivityPack = plotPerCellSensitivity(rows, figNum, styleOpts, plotAbsDVmMax_mV)
    if nargin < 3 || isempty(styleOpts)
        styleOpts = defaultFigureStyleOpts();
    end
    if nargin < 4
        plotAbsDVmMax_mV = NaN;
    end
    fs = styleOpts.fontSize;
    n = numel(rows);
    figure(figNum); clf;
    set(gcf, 'Color', 'w', 'Position', [190, 120, 650, 460]);
    ax = gca;
    hold(ax, 'on');

    cols = lines(max(n, 7));
    allX = [];
    allY = [];
    slope = nan(n, 1);
    intercept = nan(n, 1);
    r2 = nan(n, 1);
    nPts = zeros(n, 1);

    pooledFit = [NaN, NaN];
    pooledR2 = NaN;
    if numel(rows) >= 1
        for ii = 1:n
            [x, y] = sensitivityScatterXYForRow(rows(ii));
            ok = isfinite(x) & isfinite(y);
            x = x(ok);
            y = y(ok);
            if isfinite(plotAbsDVmMax_mV)
                m = abs(x) < plotAbsDVmMax_mV;
                x = x(m);
                y = y(m);
            end
            if numel(x) >= 2
                allX = [allX; x]; %#ok<AGROW>
                allY = [allY; y]; %#ok<AGROW>
            end
        end
    end
    colPool = [0.52 0.52 0.55];
    if numel(allX) >= 2
        pooledFit = polyfit(allX, allY, 1);
        xxAll = linspace(min(allX), max(allX), 200);
        yyAll = polyval(pooledFit, xxAll);
        plot(ax, xxAll, yyAll, '-', 'Color', colPool, 'LineWidth', 1.9);

        yHatAll = polyval(pooledFit, allX);
        ssResAll = sum((allY - yHatAll).^2);
        ssTotAll = sum((allY - mean(allY)).^2);
        if ssTotAll > 0
            pooledR2 = 1 - (ssResAll / ssTotAll);
        end
    end

    for ii = 1:n
        [x, y] = sensitivityScatterXYForRow(rows(ii));
        ok = isfinite(x) & isfinite(y);
        x = x(ok);
        y = y(ok);
        if isfinite(plotAbsDVmMax_mV)
            m = abs(x) < plotAbsDVmMax_mV;
            x = x(m);
            y = y(m);
        end
        nPts(ii) = numel(x);
        if isempty(x)
            continue
        end

        cLine = max(0, min(1, cols(ii, :)));
        ms = 42;
        if nPts(ii) > 25
            ms = 28;
        end
        scatter(ax, x, y, ms, cLine, 'filled', ...
            'MarkerFaceAlpha', 0.88, 'MarkerEdgeColor', [0.06 0.06 0.06], 'LineWidth', 0.65);

        if numel(x) >= 2
            pLin = polyfit(x, y, 1);
            slope(ii) = pLin(1);
            intercept(ii) = pLin(2);
            yHat = polyval(pLin, x);
            ssRes = sum((y - yHat).^2);
            ssTot = sum((y - mean(y)).^2);
            if ssTot > 0
                r2(ii) = 1 - (ssRes / ssTot);
            end
            xx = linspace(min(x), max(x), max(120, 40 * numel(x)));
            yy = polyval(pLin, xx);
            plot(ax, xx, yy, '-', 'Color', cLine, 'LineWidth', 2.85);
        end
    end

    showG = true;
    if isfield(styleOpts, 'showGrid')
        showG = styleOpts.showGrid;
    end
    if showG
        grid(ax, 'on');
    else
        grid(ax, 'off');
    end
    xlabel(ax, 'peak |\Delta V_m| per pulse (mV)', 'FontSize', fs, 'Color', [0.05 0.05 0.05], 'Interpreter', 'tex');
    ylabel(ax, 'peak filtered \Delta F/F per pulse', 'FontSize', fs, 'Color', [0.05 0.05 0.05], 'Interpreter', 'tex');
    title(ax, sprintf(['\\Delta F/F vs |\\Delta V_m| (n=%d cells).'], n), ...
        'FontWeight', 'normal', 'FontSize', fs, 'Color', [0.05 0.05 0.05], 'Interpreter', 'tex');
    hold(ax, 'off');

    if styleOpts.squareAxes
        axis(ax, 'tight');
        pbaspect(ax, [1, 1, 1]);
    end
    if isfield(styleOpts, 'xAxisZeroMin') && styleOpts.xAxisZeroMin
        axis(ax, 'tight');
        xl = xlim(ax);
        xlim(ax, [0, max(xl(2), 1e-6)]);
    end

    applyLightThemeToFigure(figure(figNum), styleOpts);

    sensitivityPack = struct();
    sensitivityPack.cellSlope = slope;
    sensitivityPack.cellIntercept = intercept;
    sensitivityPack.cellR2 = r2;
    sensitivityPack.cellNumPoints = nPts;
    sensitivityPack.pooledFit = pooledFit;
    sensitivityPack.pooledR2 = pooledR2;
end

function [x, y] = sensitivityScatterXYForRow(row)
    % One scatter point per (condition × pulse): measured peak ΔVm vs peak filtered dF/F
    % when saved by slice analysis; else one point per condition (mean-over-pulses calibration).
    hasPulse = isfield(row, 'dVmPerPulseMeas') && isfield(row, 'peakFiltDffPerPulse') ...
        && ~isempty(row.dVmPerPulseMeas) && ~isempty(row.peakFiltDffPerPulse);
    if hasPulse
        x = row.dVmPerPulseMeas(:);
        y = row.peakFiltDffPerPulse(:);
    else
        x = row.dVmCalib(:);
        y = row.peakCalib(:);
    end
end
