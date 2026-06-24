F%% VoltMapping ConnMatrix Analysis
clearvars -except CCwsAllSessions

%% Load Analysis Results
% Run this section after loading the specific cell analysis file, then
% re-run the above sections to regenerate figures

names = fieldnames(voltMapping);
for i = 1:numel(names)
    assignin('caller', names{i}, voltMapping.(names{i}));
end
nConds = length(outParams.power);

names = fieldnames(ephys);
for i = 1:numel(names)
    assignin('caller', names{i}, ephys.(names{i}));
end

Fs = voltMapping.daqParams.Fs;
imagingFreq = voltMapping.imagingFreq;
trialTime = voltMapping.daqParams.maxSweepLengthSec;
nTrials = length(voltMapping.trialCond); % ALTERNATIVELY "length(find(cellfun(@isempty, ExpStruct.inputs)==0))". Instead of "length(ExpStruct.inputs)" this puts out true number of trials successfully recorded
powers = voltMapping.outParams.power; % ALTERNATIVELY "unique(ExpStruct.trialCond)" what powers were used
nConds = length(voltMapping.outParams.power); % ALTERNATIVELY "length(unique(ExpStruct.trialCond))" total number of powers used
nHolos = voltMapping.holoStimParams.nHolos; % number of holograms in grid
pulseDurs = unique(voltMapping.outParams.pulseDur);
nPulses = unique(voltMapping.outParams.nPulses);
ipi = voltMapping.outParams.ipi;F
totalPulses = nHolos*nPulses;
% SpotCoordinates = SortedData.holoRequest.targets;
startTime = (voltMapping.holoStimParams.startTime)/1000;
imagesIndex = voltMapping.imagesIndex;
UpOrDown = voltMapping.UpOrDown;
ephysFilePath = voltMapping.ephysFilePath;
ImgsFilePath = voltMapping.ImgsFilePath;
cellID = voltMapping.cellID;

cutOffFreq = 50;   % Cutoff frequency
[b, a] = butter(4, cutOffFreq/(imagingFreq/2));

%% Publication (light mode): filt mean + 95% CI + ephys, shared y-scale per condition, scale bars
% One figure per (condition, hologram); advance with pause like the four-panel section. For each
% condition, all holos share identical xlim, left ylim (dF/F percent), and right ylim (rel. Vm, mV).
% Pre-stim mean subtracted from filt + ephys so baselines meet at y=0; symmetric ylims align zeros.
% Scale bars: dF/F + time form an L-corner; Vm vertical offset in x (0.5 mV).
cellIdxPub = double(input('which cell number (publication filt + scale bars)? '));

filtPubMean     = voltMapping.(cellID{cellIdxPub}).filtHoloSortedImagingMean_commonF0;
filtPubCI       = voltMapping.(cellID{cellIdxPub}).filtCIDffAllConds_commonF0;
FsEphysPub      = voltMapping.daqParams.Fs;
plotEphysPub    = exist('ePhysAvail', 'var') && ePhysAvail == 1 ...
    && isfield(voltMapping, 'ephys') && isfield(voltMapping.ephys, 'holoSortedDataMean');
scaleBarPct     = 0.5;  % dF/F scale bar height (percent)
scaleBar_mV     = 0.5;  % Vm scale bar height (mV)
timeBarSec      = 0.01; % horizontal time scale bar: 10 ms
baselinePrePulseMarginSec = 1e-3; % pre-stim window ends this far before 1st pulse (s)
pubFluorRgb     = [0, 0.45, 0.2];
pubFillRgb      = [0.88, 0.96, 0.88];
pubCiLineRgb    = [0.62, 0.82, 0.62];
pubEphysRgb     = [0.15, 0.15, 0.15];
pubStimRgb      = [1, 0.45, 0.45];
pubFont         = 'Arial';
pubFontSize     = 12;
pubLineW        = 2.5;    % filt mean + ephys mean (pt)

for ccPub = 3 %1:nConds
    nHP = nHolos(ccPub);
    if nHP < 1
        continue
    end

    % --- Max duration and global y-ranges (identical for every holo in this condition) ---
    maxTEnd = 0;
    yFluorMin = inf;
    yFluorMax = -inf;
    yVmMin = inf;
    yVmMax = -inf;
    pdMsP = voltMapping.outParams.pulseDur(:);
    pulseDurMsP = pdMsP(min(ccPub, numel(pdMsP)));
    if pulseDurMsP <= 0 && exist('pulseDurs', 'var') && ~isempty(pulseDurs)
        pulseDurMsP = pulseDurs(1);
    end
    pulseDurSecP = pulseDurMsP / 1000;

    for hhP = 1:nHP
        nT = numel(filtPubMean{ccPub}{hhP});
        tAxG = linspace(0, nT / imagingFreq, nT);
        maxTEnd = max(maxTEnd, nT / imagingFreq);
        mPct = filtPubMean{ccPub}{hhP}(:) * 100;
        ciP = filtPubCI{ccPub}{hhP, 1} * 100;
        if isempty(nPulseCoords)
            idxB = false(size(tAxG));
            nb = min(max(3, round(0.08 * nT)), nT);
            idxB(1:nb) = true;
        else
            t1 = nPulseCoords(1) / FsEphysPub;
            tHi = max(tAxG(1), t1 - baselinePrePulseMarginSec);
            idxB = tAxG <= tHi;
            if sum(idxB) < 3
                nb = min(max(3, round(0.08 * nT)), nT);
                idxB = false(size(tAxG));
                idxB(1:nb) = true;
            end
        end
        blFl = mean(mPct(idxB), 'omitnan');
        mAdj = mPct - blFl;
        ciLo = ciP(:, 1) - blFl;
        ciHi = ciP(:, 2) - blFl;
        yFluorMin = min([yFluorMin; mAdj(:); ciLo(:); ciHi(:)], [], 'omitnan');
        yFluorMax = max([yFluorMax; mAdj(:); ciLo(:); ciHi(:)], [], 'omitnan');

        if plotEphysPub && hhP <= size(voltMapping.ephys.holoSortedDataMean{ccPub}, 2)
            eMC = voltMapping.ephys.holoSortedDataMean{ccPub}(:, hhP);
            nEp = numel(eMC);
            if nEp > 1
                tEp = linspace(0, nEp / FsEphysPub, nEp);
                eL = interp1(tEp(:), double(eMC(:)), tAxG(:), 'linear', 'extrap');
            elseif nEp == 1
                eL = repmat(double(eMC), numel(tAxG), 1);
            else
                eL = [];
            end
            if ~isempty(eL)
                blE = mean(eL(idxB), 'omitnan');
                eAdj = eL(:) - blE;
                yVmMin = min(yVmMin, min(eAdj, [], 'omitnan'));
                yVmMax = max(yVmMax, max(eAdj, [], 'omitnan'));
            end
        end
    end

    aL = max(abs([yFluorMin; yFluorMax]));
    padF = 0.06 * max(aL, eps);
    yLimL = [-aL - padF, aL + padF];
    if plotEphysPub && isfinite(yVmMin) && isfinite(yVmMax)
        aV = max(abs([yVmMin; yVmMax]));
        padV = 0.06 * max(aV, eps);
        yLimR = [-aV - padV, aV + padV];
    else
        yLimR = [];
    end
    xLimAll = [0, maxTEnd];

    for hhP = 1:nHP
        figPub = figure(43000 + ccPub * 1000 + hhP);
        clf(figPub)
        set(figPub, 'Color', 'w', 'InvertHardcopy', 'off', 'Position', [120, 120, 560, 560]);
        axh = axes('Parent', figPub, 'Position', [0.14, 0.14, 0.72, 0.72]);
        set(axh, 'Color', 'w', 'Box', 'off', ...
            'FontName', pubFont, 'FontSize', pubFontSize, ...
            'XColor', 'none', 'YColor', 'none', 'TickDir', 'out', ...
            'XTick', [], 'YTick', [], 'XTickLabel', [], 'YTickLabel', []);

        nT = numel(filtPubMean{ccPub}{hhP});
        tAx = linspace(0, nT / imagingFreq, nT);
        mPct = filtPubMean{ccPub}{hhP}(:) * 100;
        ciP = filtPubCI{ccPub}{hhP, 1} * 100;
        if isempty(nPulseCoords)
            idxB = false(size(tAx));
            nb = min(max(3, round(0.08 * nT)), nT);
            idxB(1:nb) = true;
        else
            t1 = nPulseCoords(1) / FsEphysPub;
            tHi = max(tAx(1), t1 - baselinePrePulseMarginSec);
            idxB = tAx <= tHi;
            if sum(idxB) < 3
                nb = min(max(3, round(0.08 * nT)), nT);
                idxB = false(size(tAx));
                idxB(1:nb) = true;
            end
        end
        blFl = mean(mPct(idxB), 'omitnan');
        mPlot = mPct - blFl;
        ciLo = ciP(:, 1) - blFl;
        ciHi = ciP(:, 2) - blFl;

        ephysLine = [];
        if plotEphysPub && hhP <= size(voltMapping.ephys.holoSortedDataMean{ccPub}, 2)
            eMC = voltMapping.ephys.holoSortedDataMean{ccPub}(:, hhP);
            nEp = numel(eMC);
            if nEp > 1
                tEp = linspace(0, nEp / FsEphysPub, nEp);
                ephysLine = interp1(tEp(:), double(eMC(:)), tAx(:), 'linear', 'extrap');
            elseif nEp == 1
                ephysLine = repmat(double(eMC), numel(tAx), 1);
            end
        end
        if ~isempty(ephysLine)
            blE = mean(ephysLine(idxB), 'omitnan');
            ephysLine = ephysLine(:) - blE;
        end

        yyaxis(axh, 'left')
        hold(axh, 'on')
        xlim(axh, xLimAll)
        ylim(axh, yLimL)
        fill(axh, [tAx, fliplr(tAx)], [ciLo', fliplr(ciHi')], pubFillRgb, ...
            'EdgeColor', pubFillRgb, 'FaceAlpha', 0.9);
        for pulseIdxP = 1:length(nPulseCoords)
            tOnP = nPulseCoords(pulseIdxP) / FsEphysPub;
            patch(axh, [tOnP, tOnP + pulseDurSecP, tOnP + pulseDurSecP, tOnP], ...
                [yLimL(1), yLimL(1), yLimL(2), yLimL(2)], pubStimRgb, ...
                'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HitTest', 'off');
        end
        plot(axh, tAx, ciLo, '--', 'Color', pubCiLineRgb, 'LineWidth', 0.8);
        plot(axh, tAx, ciHi, '--', 'Color', pubCiLineRgb, 'LineWidth', 0.8);
        plot(axh, tAx, mPlot, '-', 'Color', pubFluorRgb, 'LineWidth', pubLineW);

        yyaxis(axh, 'left')
        ylblL = ylabel(axh, '');
        set(ylblL, 'Visible', 'off');

        if ~isempty(ephysLine) && ~isempty(yLimR)
            yyaxis(axh, 'right')
            plot(axh, tAx, ephysLine(:).', '-', 'Color', pubEphysRgb, 'LineWidth', pubLineW);
            ylim(axh, yLimR)
            xlim(axh, xLimAll)
            ylblR = ylabel(axh, '');
            set(ylblR, 'Visible', 'off');
            axh.YAxis(2).Color = 'k';
            axh.YAxis(2).TickValues = [];
        end

        yyaxis(axh, 'left')
        xlabel(axh, '')
        try
            pbaspect(axh, [1, 1, 1])
        catch %#ok<*CTCH>
        end

        % --- Scale bars: dF/F (vertical) + time (horizontal) share a corner; Vm vertical offset in x ---
        yyaxis(axh, 'left')
        ylim(axh, yLimL)
        xlim(axh, xLimAll)
        xlR = xLimAll(2) - xLimAll(1);
        ylRng = yLimL(2) - yLimL(1);
        xCorner = xLimAll(1) + 0.045 * xlR;
        yCorner = yLimL(1) + 0.16 * ylRng;
        xBarVmL = xCorner + 0.06 * xlR;
        plot(axh, [xCorner, xCorner], yCorner + [0, scaleBarPct], '-', 'Color', 'k', ...
            'LineWidth', 1.2, 'Clipping', 'off', 'HandleVisibility', 'off');
        text(axh, xCorner - 0.026 * xlR, yCorner + 0.5 * scaleBarPct, sprintf('%.1f%s dF/F', scaleBarPct, '%'), ...
            'FontName', pubFont, 'FontSize', pubFontSize, 'Color', 'k', ...
            'Rotation', 90, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', ...
            'Interpreter', 'none', 'Clipping', 'off');

        tBarW = min(timeBarSec, 0.85 * xlR);
        plot(axh, xCorner + [0, tBarW], [yCorner, yCorner], '-', 'Color', 'k', ...
            'LineWidth', 1.2, 'Clipping', 'off', 'HandleVisibility', 'off');
        text(axh, xCorner + 0.5 * tBarW, yCorner - 0.035 * ylRng, sprintf('%.0f ms', tBarW * 1000), ...
            'FontName', pubFont, 'FontSize', pubFontSize - 0.5, 'Color', 'k', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'Interpreter', 'none', 'Clipping', 'off');

        if ~isempty(ephysLine) && ~isempty(yLimR)
            ylRngV = yLimR(2) - yLimR(1);
            vFracCorner = (yCorner - yLimL(1)) / max(ylRng, eps);
            yBarBotV = yLimR(1) + vFracCorner * ylRngV;
            yyaxis(axh, 'right')
            ylim(axh, yLimR)
            xlim(axh, xLimAll)
            plot(axh, [xBarVmL, xBarVmL], yBarBotV + [0, scaleBar_mV], '-', 'Color', 'k', ...
                'LineWidth', 1.2, 'Clipping', 'off', 'HandleVisibility', 'off');
            ymidVm = yBarBotV + 0.5 * scaleBar_mV;
            vFrac = (ymidVm - yLimR(1)) / max(yLimR(2) - yLimR(1), eps);
            yLblVm = yLimL(1) + vFrac * ylRng;
            yyaxis(axh, 'left')
            text(axh, xBarVmL + 0.020 * xlR, yLblVm, sprintf('%.1f mV', scaleBar_mV), ...
                'FontName', pubFont, 'FontSize', pubFontSize, 'Color', 'k', ...
                'Rotation', 90, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', ...
                'Interpreter', 'none', 'Clipping', 'off');
        end

        yyaxis(axh, 'left')
        hold(axh, 'off')

        axh.XAxis.Visible = 'off';
        axh.YAxis(1).Visible = 'off';
        if numel(axh.YAxis) > 1
            axh.YAxis(2).Visible = 'off';
        end

        % sgtitle(figPub, sprintf(['Cell %d — cond %d, holo %d (filt+CI, ephys; pre-stim mean removed, ', ...
        %     'symmetric axes; scales shared within cond)'], cellIdxPub, ccPub, hhP), ...
        %     'FontName', pubFont, 'FontSize', pubFontSize + 1, 'FontWeight', 'normal', ...
        %     'Interpreter', 'none');
        pause
    end
end

%% Connection Matrix
connMatrix = cell(nConds, 1);
connPvalMatrix = cell(nConds, 1);
respMatrix = cell(nConds, 1);
respMatrixConnectedCells = cell(nConds, 1);
percentConn = cell(nConds, 1);
connMatrixEphys = cell(nConds, 1);
connPvalEphys = cell(nConds, 1);
respMatrixEphys = cell(nConds, 1);
respMatrixEphysConnected = cell(nConds, 1);
percentConnEphys = cell(nConds, 1);
for cc = 1:nConds
    connMatrix{cc} = nan(nCells, nHolos(cc));
    connPvalMatrix{cc} = nan(nCells, nHolos(cc));
    respMatrix{cc} = nan(nCells, nHolos(cc));
    respMatrixConnectedCells{cc} = nan(nCells, nHolos(cc));
    connMatrixEphys{cc} = nan(1, nHolos(cc));
    connPvalEphys{cc} = nan(1, nHolos(cc));
    respMatrixEphys{cc} = nan(1, nHolos(cc));
    respMatrixEphysConnected{cc} = nan(1, nHolos(cc));
end

% --- Imaging-only connection stats (simplified) ---
% Paired one-tailed test: mean(stim window) > mean(prestim window) per trial.
% connPvalMatrix: raw one-sided p-values; connMatrix: 1 if p < alphaConn (increase), else 0; NaN = not testable.

preStimIdx  = 1:floor(preStimWindow/1000*imagingFreq);
stimStart = floor(preStimWindow/1000*imagingFreq) + 1;
stimIdx = stimStart:(floor(preStimWindow/1000*imagingFreq+ipi*(nPulses-3)/1000*imagingFreq));

alphaConn = 0.05;

for nn = 1:nCells
    clear exclFiltHoloSortedImagingAllTrials
    exclFiltHoloSortedImagingAllTrials = voltMapping.(cellID{nn}).exclFiltHoloSortedImagingAllTrials_commonF0;
    % exclFiltHoloSortedImagingAllTrials = voltMapping.(cellID{nn}).exclFiltHoloSortedImagingAllTrials;

    for cc = 1:nConds
        for hh = 1:nHolos(cc)
            testData = exclFiltHoloSortedImagingAllTrials{cc}{hh};
            if isempty(testData) || max([preStimIdx stimIdx]) > size(testData, 1)
                connPvalMatrix{cc}(nn, hh) = NaN;
                connMatrix{cc}(nn, hh) = NaN;
                respMatrix{cc}(nn, hh) = NaN;
                respMatrixConnectedCells{cc}(nn, hh) = NaN;
                continue;
            end

            analysisRows = unique([preStimIdx stimIdx]);
            validCols = all(~isnan(testData(analysisRows, :)), 1);
            validData = testData(:, validCols);
            if size(validData, 2) < 2
                connPvalMatrix{cc}(nn, hh) = NaN;
                connMatrix{cc}(nn, hh) = NaN;
                respMatrix{cc}(nn, hh) = NaN;
                respMatrixConnectedCells{cc}(nn, hh) = NaN;
                continue;
            end

            meanBaseline = mean(validData(preStimIdx, :), 1);
            meanStim = mean(validData(stimIdx, :), 1);
            deltaTrial = meanStim - meanBaseline;
            respMatrix{cc}(nn, hh) = mean(deltaTrial, "omitnan");

            validPair = ~isnan(meanBaseline) & ~isnan(meanStim);
            if sum(validPair) < 2
                connPvalMatrix{cc}(nn, hh) = NaN;
                connMatrix{cc}(nn, hh) = NaN;
                respMatrix{cc}(nn, hh) = NaN;
                respMatrixConnectedCells{cc}(nn, hh) = NaN;
                continue;
            end

            [~, connPvalMatrix{cc}(nn, hh)] = ttest(meanStim(validPair) - meanBaseline(validPair), 0, ...
                'Alpha', alphaConn, 'Tail', 'right');

            connMatrix{cc}(nn, hh) = double(connPvalMatrix{cc}(nn, hh) < alphaConn);
            if connMatrix{cc}(nn, hh) == 1
                respMatrixConnectedCells{cc}(nn, hh) = respMatrix{cc}(nn, hh);
            else
                respMatrixConnectedCells{cc}(nn, hh) = NaN;
            end
        end

        pRow = connPvalMatrix{cc}(nn, :);
        analyzableMask = ~isnan(pRow);
        if any(analyzableMask)
            percentConn{cc}(nn, 1) = sum(connMatrix{cc}(nn, analyzableMask) == 1) / sum(analyzableMask);
        else
            percentConn{cc}(nn, 1) = NaN;
        end
    end
end

%% Plot imaging connectivity matrix — gridded (same layout as rows=cells, cols=presynaptic holos)
for cc = 1:nConds
    figConn = figure(199 + cc);
    clf(figConn);
    set(figConn, 'Position', [100, 100, 1600, 300], 'Color', 'w', 'InvertHardcopy', 'off');

    M = connMatrix{cc};
    Mplot = M;
    Mplot(isnan(Mplot)) = 0.5; % mid-gray: not testable (NaN); 0 = not sig, 1 = sig

    imagesc(Mplot);
    colormap(gray);
    clim([0 1]);
    axis equal tight;

    hold on
    [nRows, nCols] = size(M);
    for r = 0.5:1:nRows + 0.5
        plot([0.5, nCols + 0.5], [r, r], 'k-', 'LineWidth', 0.75);
    end
    for c = 0.5:1:nCols + 0.5
        plot([c, c], [0.5, nRows + 0.5], 'k-', 'LineWidth', 0.75);
    end
    hold off

    cbConnImg = colorbar;
    set(cbConnImg, 'FontSize', 9, 'Color', 'k');
    if isprop(cbConnImg, 'FontColor')
        cbConnImg.FontColor = 'k';
    end
    set(gca, 'TickLength', [0 0], 'FontSize', 16, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    xticks(1:nCols);
    yticks(1:nRows);
    xlabel('Presynaptic target (hologram)');
    ylabel('Cell (ROI)');
    title(sprintf(['Connectivity — paired  t, ', ...
        'p < 0.05 '], cc, alphaConn), 'Color', 'k');
end

%% Trial raster viewer (choose condition/hologram)
% Visualize trial-wise imaging traces for a selected holo as a raster.
% Imaging uses one selected ROI (set rasterCellToPlot).
showTrialRaster = true;
if showTrialRaster
    rasterCondToPlot = 1; % set [] to choose interactively in command window
    rasterHoloToPlot = 1; % set [] to choose interactively in command window
    rasterCellToPlot = 1; % imaging ROI index

    rasterCellToPlot = min(max(1, rasterCellToPlot), nCells);

    if isempty(rasterCondToPlot)
        promptCond = sprintf('Condition index to view (1-%d): ', nConds);
        rasterCondToPlot = input(promptCond);
    end
    rasterCondToPlot = round(rasterCondToPlot);
    if rasterCondToPlot < 1 || rasterCondToPlot > nConds
        warning('Trial raster skipped: condition index %d is out of range.', rasterCondToPlot);
    else
        if isempty(rasterHoloToPlot)
            promptHolo = sprintf('Hologram index to view for cond %d (1-%d): ', rasterCondToPlot, nHolos(rasterCondToPlot));
            rasterHoloToPlot = input(promptHolo);
        end
        rasterHoloToPlot = round(rasterHoloToPlot);

        if rasterHoloToPlot < 1 || rasterHoloToPlot > nHolos(rasterCondToPlot)
            warning('Trial raster skipped: hologram index %d is out of range for cond %d.', ...
                rasterHoloToPlot, rasterCondToPlot);
        else
            imgTrials = voltMapping.(cellID{rasterCellToPlot}).exclFiltHoloSortedImagingAllTrials_commonF0{rasterCondToPlot}{rasterHoloToPlot};

            if isempty(imgTrials)
                warning('Trial raster skipped: no imaging trial data for cond %d holo %d.', ...
                    rasterCondToPlot, rasterHoloToPlot);
            else
                validImgCols = all(~isnan(imgTrials), 1);
                imgTrials = imgTrials(:, validImgCols);

                figure(105);
                clf;
                set(gcf, 'Position', [100, 100, 400, 520]);
                set(gcf, 'Color', 'w', 'InvertHardcopy', 'off');

                nSImg = size(imgTrials, 1);
                nTrialImg = size(imgTrials, 2);
                tMsImg = (0:nSImg-1) / imagingFreq * 1000;
                firstPulseMsImg = (stimStart - 1) / imagingFreq * 1000;
                imagesc(tMsImg, 1:nTrialImg, imgTrials');
                clim([0 0.1]);
                hold on;
                xline(firstPulseMsImg, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
                hold off;
                xlabel('Time (ms)');
                ylabel('Trial');
                title(sprintf('Imaging raster (ROI %d: %s; cond %d, holo %d)', ...
                    rasterCellToPlot, cellID{rasterCellToPlot}, rasterCondToPlot, rasterHoloToPlot));
                axImg = gca;
                axImg.Title.Position(2) = axImg.Title.Position(2) - 1.5;
                set(axImg, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontSize', 12);
                axImg.Title.Color = 'k';
                colormap(gca, 'parula');
                cbImg = colorbar;
                set(cbImg, 'FontSize', 9, 'Color', 'k');
                if isprop(cbImg, 'FontColor')
                    cbImg.FontColor = 'k';
                end
            end
        end
    end
end