%% VoltMapping ConnMatrix Analysis
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
ipi = voltMapping.outParams.ipi;
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
% Scale bars: 1%% dF/F, 0.5 mV Vm (both vertical bars on the left), 10 ms horizontal.
cellIdxPub = double(input('which cell number (publication filt + scale bars)? '));

filtPubMean     = voltMapping.(cellID{cellIdxPub}).filtHoloSortedImagingMean;
filtPubCI       = voltMapping.(cellID{cellIdxPub}).filtCIDffAllConds;
FsEphysPub      = voltMapping.daqParams.Fs;
plotEphysPub    = exist('ePhysAvail', 'var') && ePhysAvail == 1 ...
    && isfield(voltMapping, 'ephys') && isfield(voltMapping.ephys, 'holoSortedDataMean');
scaleBarPct     = 1;    % dF/F scale bar height (percent)
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

for ccPub = 2 %1:nConds
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

        % --- Scale bars: both vertical segments on left; time bar horizontal ---
        yyaxis(axh, 'left')
        ylim(axh, yLimL)
        xlim(axh, xLimAll)
        xlR = xLimAll(2) - xLimAll(1);
        ylRng = yLimL(2) - yLimL(1);
        xBarPct = xLimAll(1) + 0.04 * xlR;
        xBarVmL = xLimAll(1) + 0.10 * xlR;
        yBarBotL = yLimL(1) + 0.08 * ylRng;
        plot(axh, [xBarPct, xBarPct], yBarBotL + [0, scaleBarPct], '-', 'Color', 'k', ...
            'LineWidth', 1.2, 'Clipping', 'off', 'HandleVisibility', 'off');
        text(axh, xBarPct + 0.020 * xlR, yBarBotL + 0.5 * scaleBarPct, sprintf('1%s dF/F', '%'), ...
            'FontName', pubFont, 'FontSize', pubFontSize, 'Color', 'k', ...
            'Rotation', 90, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', ...
            'Interpreter', 'none', 'Clipping', 'off');

        tBarW = min(timeBarSec, 0.85 * xlR);
        yBarTime = yLimL(1) + 0.03 * ylRng;
        plot(axh, xBarPct + [0, tBarW], [yBarTime, yBarTime], '-', 'Color', 'k', ...
            'LineWidth', 1.2, 'Clipping', 'off', 'HandleVisibility', 'off');
        text(axh, xBarPct + 0.5 * tBarW, yBarTime - 0.04 * ylRng, sprintf('%.0f ms', tBarW * 1000), ...
            'FontName', pubFont, 'FontSize', pubFontSize - 0.5, 'Color', 'k', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'Interpreter', 'none', 'Clipping', 'off');

        if ~isempty(ephysLine) && ~isempty(yLimR)
            ylRngV = yLimR(2) - yLimR(1);
            yBarBotV = yLimR(1) + 0.08 * ylRngV;
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

        sgtitle(figPub, sprintf(['Cell %d — cond %d, holo %d (filt+CI, ephys; pre-stim mean removed, ', ...
            'symmetric axes; scales shared within cond)'], cellIdxPub, ccPub, hhP), ...
            'FontName', pubFont, 'FontSize', pubFontSize + 1, 'FontWeight', 'normal', ...
            'Interpreter', 'none');
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

preStimIdx  = 1:floor(preStimWindow/1000*imagingFreq);
stimStart = floor(preStimWindow/1000*imagingFreq) + 1;
stimIdx = stimStart:(floor(preStimWindow/1000*imagingFreq+ipi*(nPulses)/1000*imagingFreq));
% Raw p-values from paired t-tests; decisions use BH-FDR across holos (per cell & condition).
alphaConn = 0.05;
fdrQConn = 0.01;
% Optional: require mean trial delta (stim - baseline) >= this to call "connected" (NaN = no minimum).
connMinDeltaResp = 0.003;

% Robust ephys availability check from struct + flag.
ePhysAvailInStruct = isfield(voltMapping, 'ePhysAvail') && isequal(voltMapping.ePhysAvail, 1);
ePhysAvailInWorkspace = exist('ePhysAvail', 'var') && isequal(ePhysAvail, 1);
hasEphysField = isfield(voltMapping, 'ephys') && ~isempty(voltMapping.ephys);
useEphysForConn = hasEphysField && (ePhysAvailInStruct || ePhysAvailInWorkspace);

ephysTrialsByCond = [];
if useEphysForConn
    if isfield(voltMapping.ephys, 'exclFiltHoloSortedDataAllTrials') && ~isempty(voltMapping.ephys.exclFiltHoloSortedDataAllTrials)
        ephysTrialsByCond = voltMapping.ephys.exclFiltHoloSortedDataAllTrials;
    elseif isfield(voltMapping.ephys, 'exclHoloSortedDataAllTrials') && ~isempty(voltMapping.ephys.exclHoloSortedDataAllTrials)
        ephysTrialsByCond = voltMapping.ephys.exclHoloSortedDataAllTrials;
    elseif isfield(voltMapping.ephys, 'holoSortedDataAllTrials') && ~isempty(voltMapping.ephys.holoSortedDataAllTrials)
        ephysTrialsByCond = voltMapping.ephys.holoSortedDataAllTrials;
    else
        useEphysForConn = false;
        warning('Ephys flagged as available but no trial-wise ephys arrays were found in voltMapping.ephys.');
    end
end

if useEphysForConn && isfield(voltMapping.ephys, 'preStimWindow') && ~isempty(voltMapping.ephys.preStimWindow)
    preStimWindowEphys = voltMapping.ephys.preStimWindow;
else
    preStimWindowEphys = preStimWindow;
end
preStimIdxEphys = 1:floor(preStimWindowEphys/1000*Fs);
stimStartEphys = floor(preStimWindowEphys/1000*Fs) + 1;
stimIdxEphys = stimStartEphys:(ceil(preStimWindowEphys/1000*Fs+ipi*(nPulses-2)/1000*Fs));

for nn = 1:nCells
    clear exclFiltHoloSortedImagingAllTrials
    exclFiltHoloSortedImagingAllTrials = voltMapping.(cellID{nn}).exclFiltHoloSortedImagingAllTrials_commonF0;
    for cc = 1:nConds
        for hh = 1:nHolos(cc)
            testData = exclFiltHoloSortedImagingAllTrials{cc}{hh};
            if isempty(testData) || max([preStimIdx stimIdx]) > size(testData, 1)
                connPvalMatrix{cc}(nn, hh) = NaN;
                respMatrix{cc}(nn, hh) = NaN;
                respMatrixConnectedCells{cc}(nn, hh) = NaN;
                continue;
            end

            % Keep only trials with finite values in the analysis windows.
            analysisRows = unique([preStimIdx stimIdx]);
            validCols = all(~isnan(testData(analysisRows, :)), 1);
            validData = testData(:, validCols);
            if size(validData, 2) < 2
                connPvalMatrix{cc}(nn, hh) = NaN;
                respMatrix{cc}(nn, hh) = NaN;
                respMatrixConnectedCells{cc}(nn, hh) = NaN;
                continue;
            end

            meanBaseline = mean(validData(preStimIdx, :), 1);
            meanStim = mean(validData(stimIdx, :), 1);
            deltaResp = meanStim - meanBaseline;
            respMatrix{cc}(nn, hh) = mean(deltaResp, "omitnan");

            validPair = ~isnan(meanBaseline) & ~isnan(meanStim);
            if sum(validPair) < 2
                connPvalMatrix{cc}(nn, hh) = NaN;
                respMatrix{cc}(nn, hh) = NaN;
                respMatrixConnectedCells{cc}(nn, hh) = NaN;
                continue;
            end

            [~, connPvalMatrix{cc}(nn, hh)] = ttest(meanBaseline(validPair), meanStim(validPair));
        end

        % Benjamini-Hochberg FDR across holograms for this cell & condition.
        pRow = connPvalMatrix{cc}(nn, :);
        ampRow = respMatrix{cc}(nn, :);
        validTest = ~isnan(pRow);
        rejRow = connMatrix_bhFdrReject(pRow, fdrQConn);
        sigRow = rejRow;
        if ~isnan(connMinDeltaResp)
            sigRow = sigRow & (ampRow >= connMinDeltaResp);
        end
        connMatrix{cc}(nn, :) = NaN;
        connMatrix{cc}(nn, validTest) = double(sigRow(validTest));
        respMatrixConnectedCells{cc}(nn, :) = NaN;
        respMatrixConnectedCells{cc}(nn, validTest) = 0;
        maskConn = validTest & sigRow;
        respMatrixConnectedCells{cc}(nn, maskConn) = ampRow(maskConn);

        analyzableMask = ~isnan(pRow);
        if any(analyzableMask)
            percentConn{cc}(nn, 1) = sum(connMatrix{cc}(nn, analyzableMask) == 1) / sum(analyzableMask);
        else
            percentConn{cc}(nn, 1) = NaN;
        end
    end
end

% Ephys-only connectivity matrix (one patched cell trace per hologram).
if useEphysForConn
    for cc = 1:nConds
        for hh = 1:nHolos(cc)
            if cc > numel(ephysTrialsByCond) || hh > numel(ephysTrialsByCond{cc})
                continue;
            end
            testDataE = ephysTrialsByCond{cc}{hh};
            if isempty(testDataE) || max([preStimIdxEphys stimIdxEphys]) > size(testDataE, 1)
                continue;
            end

            analysisRowsE = unique([preStimIdxEphys stimIdxEphys]);
            validColsE = all(~isnan(testDataE(analysisRowsE, :)), 1);
            validDataE = testDataE(:, validColsE);
            if size(validDataE, 2) < 2
                continue;
            end

            meanBaselineE = mean(validDataE(preStimIdxEphys, :), 1);
            meanStimE = mean(validDataE(stimIdxEphys, :), 1);
            deltaRespE = meanStimE - meanBaselineE;
            respMatrixEphys{cc}(1, hh) = mean(deltaRespE, "omitnan");

            validPairE = ~isnan(meanBaselineE) & ~isnan(meanStimE);
            if sum(validPairE) < 2
                connPvalEphys{cc}(1, hh) = NaN;
                respMatrixEphys{cc}(1, hh) = NaN;
                continue;
            end

            [~, connPvalEphys{cc}(1, hh)] = ttest(meanBaselineE(validPairE), meanStimE(validPairE));
        end

        pRowE = connPvalEphys{cc}(1, :);
        ampRowE = respMatrixEphys{cc}(1, :);
        validTestE = ~isnan(pRowE);
        rejE = connMatrix_bhFdrReject(pRowE, fdrQConn);
        sigE = rejE;
        if ~isnan(connMinDeltaResp)
            sigE = sigE & (ampRowE >= connMinDeltaResp);
        end
        connMatrixEphys{cc}(1, :) = nan(1, nHolos(cc));
        connMatrixEphys{cc}(1, validTestE) = double(sigE(validTestE));
        respMatrixEphysConnected{cc}(1, :) = nan(1, nHolos(cc));
        respMatrixEphysConnected{cc}(1, validTestE) = 0;
        maskConnE = validTestE & sigE;
        respMatrixEphysConnected{cc}(1, maskConnE) = ampRowE(maskConnE);

        analyzableMaskE = ~isnan(pRowE);
        if any(analyzableMaskE)
            percentConnEphys{cc}(1, 1) = sum(connMatrixEphys{cc}(1, analyzableMaskE) == 1) / sum(analyzableMaskE);
        else
            percentConnEphys{cc}(1, 1) = NaN;
        end
    end
end

% Compare imaging vs ephys connectivity on paired ROI when available.
connAgreement = struct();
connAgreement.alpha = alphaConn;
connAgreement.fdrQ = fdrQConn;
connAgreement.connMinDeltaResp = connMinDeltaResp;
connAgreement.multipleTesting = 'BH-FDR across holograms per cell (imaging) or per row (ephys)';
connAgreement.useEphysForConn = useEphysForConn;
connAgreement.pairedRoiIdx = NaN;
connAgreement.perCondition = cell(nConds, 1);
connAgreement.global = struct();

if useEphysForConn
    if isfield(voltMapping.ephys, 'pairedRoiIndex') && isscalar(voltMapping.ephys.pairedRoiIndex) ...
            && ~isnan(voltMapping.ephys.pairedRoiIndex)
        pairedRoiIdx = round(voltMapping.ephys.pairedRoiIndex);
    elseif nCells == 1
        pairedRoiIdx = 1;
    else
        pairedRoiIdx = NaN;
    end

    connAgreement.pairedRoiIdx = pairedRoiIdx;
    allImgBinary = [];
    allEphysBinary = [];
    allImgAmp = [];
    allEphysAmp = [];

    if ~isnan(pairedRoiIdx) && pairedRoiIdx >= 1 && pairedRoiIdx <= nCells
        for cc = 1:nConds
            imgConn = connMatrix{cc}(pairedRoiIdx, :);
            ephConn = connMatrixEphys{cc}(1, :);
            validMask = ~isnan(imgConn) & ~isnan(ephConn);
            imgBin = imgConn(validMask) == 1;
            ephBin = ephConn(validMask) == 1;

            TP = sum(imgBin & ephBin);
            TN = sum(~imgBin & ~ephBin);
            FP = sum(imgBin & ~ephBin);
            FN = sum(~imgBin & ephBin);
            nValid = numel(imgBin);

            condStats = struct();
            condStats.nCompared = nValid;
            condStats.TP = TP;
            condStats.TN = TN;
            condStats.FP = FP;
            condStats.FN = FN;
            if nValid > 0
                condStats.accuracy = (TP + TN) / nValid;
                condStats.sensitivity = TP / max(TP + FN, eps);
                condStats.specificity = TN / max(TN + FP, eps);
                condStats.precision = TP / max(TP + FP, eps);
                condStats.f1 = 2 * TP / max(2 * TP + FP + FN, eps);
            else
                condStats.accuracy = NaN;
                condStats.sensitivity = NaN;
                condStats.specificity = NaN;
                condStats.precision = NaN;
                condStats.f1 = NaN;
            end

            discordant = FP + FN;
            condStats.mcnemarB = FP;
            condStats.mcnemarC = FN;
            condStats.mcnemarChi2 = NaN;
            condStats.mcnemarP = NaN;
            if discordant > 0
                condStats.mcnemarChi2 = (abs(FP - FN) - 1)^2 / discordant;
                condStats.mcnemarP = min(1, 2 * binocdf(min(FP, FN), discordant, 0.5));
            end

            imgAmp = respMatrix{cc}(pairedRoiIdx, :);
            ephAmp = respMatrixEphys{cc}(1, :);
            ampMask = ~isnan(imgAmp) & ~isnan(ephAmp);
            condStats.imgAmp = imgAmp(ampMask);
            condStats.ephysAmp = ephAmp(ampMask);
            condStats.nAmpCompared = sum(ampMask);
            if condStats.nAmpCompared >= 3
                zImg = zscore(condStats.imgAmp);
                zEphys = zscore(condStats.ephysAmp);
                [~, condStats.pairedTtestP] = ttest(zImg, zEphys);
                condStats.signrankP = signrank(zImg, zEphys);
                condStats.ampCorrR = corr(condStats.imgAmp(:), condStats.ephysAmp(:), "Rows", "complete");
            else
                condStats.pairedTtestP = NaN;
                condStats.signrankP = NaN;
                condStats.ampCorrR = NaN;
            end

            connAgreement.perCondition{cc} = condStats;
            allImgBinary = [allImgBinary, imgBin]; %#ok<AGROW>
            allEphysBinary = [allEphysBinary, ephBin]; %#ok<AGROW>
            allImgAmp = [allImgAmp, condStats.imgAmp]; %#ok<AGROW>
            allEphysAmp = [allEphysAmp, condStats.ephysAmp]; %#ok<AGROW>
        end

        TP = sum(allImgBinary & allEphysBinary);
        TN = sum(~allImgBinary & ~allEphysBinary);
        FP = sum(allImgBinary & ~allEphysBinary);
        FN = sum(~allImgBinary & allEphysBinary);
        nValidAll = numel(allImgBinary);
        connAgreement.global.nCompared = nValidAll;
        connAgreement.global.TP = TP;
        connAgreement.global.TN = TN;
        connAgreement.global.FP = FP;
        connAgreement.global.FN = FN;
        if nValidAll > 0
            connAgreement.global.accuracy = (TP + TN) / nValidAll;
            connAgreement.global.sensitivity = TP / max(TP + FN, eps);
            connAgreement.global.specificity = TN / max(TN + FP, eps);
            connAgreement.global.precision = TP / max(TP + FP, eps);
            connAgreement.global.f1 = 2 * TP / max(2 * TP + FP + FN, eps);
        else
            connAgreement.global.accuracy = NaN;
            connAgreement.global.sensitivity = NaN;
            connAgreement.global.specificity = NaN;
            connAgreement.global.precision = NaN;
            connAgreement.global.f1 = NaN;
        end

        discordantAll = FP + FN;
        connAgreement.global.mcnemarB = FP;
        connAgreement.global.mcnemarC = FN;
        connAgreement.global.mcnemarChi2 = NaN;
        connAgreement.global.mcnemarP = NaN;
        if discordantAll > 0
            connAgreement.global.mcnemarChi2 = (abs(FP - FN) - 1)^2 / discordantAll;
            connAgreement.global.mcnemarP = min(1, 2 * binocdf(min(FP, FN), discordantAll, 0.5));
        end

        connAgreement.global.nAmpCompared = numel(allImgAmp);
        if connAgreement.global.nAmpCompared >= 3
            zImgAll = zscore(allImgAmp);
            zEphysAll = zscore(allEphysAmp);
            [~, connAgreement.global.pairedTtestP] = ttest(zImgAll, zEphysAll);
            connAgreement.global.signrankP = signrank(zImgAll, zEphysAll);
            connAgreement.global.ampCorrR = corr(allImgAmp(:), allEphysAmp(:), "Rows", "complete");
        else
            connAgreement.global.pairedTtestP = NaN;
            connAgreement.global.signrankP = NaN;
            connAgreement.global.ampCorrR = NaN;
        end

        fprintf('Imaging vs ephys (paired ROI %d): accuracy=%.3f, sens=%.3f, spec=%.3f, McNemar p=%.4g\n', ...
            pairedRoiIdx, connAgreement.global.accuracy, connAgreement.global.sensitivity, ...
            connAgreement.global.specificity, connAgreement.global.mcnemarP);
        fprintf('Amplitude comparison (z-scored holo means): paired t-test p=%.4g, signrank p=%.4g, corr r=%.3f\n', ...
            connAgreement.global.pairedTtestP, connAgreement.global.signrankP, connAgreement.global.ampCorrR);
    else
        warning('Ephys available, but no valid paired ROI index found. Cross-modality stats were skipped.');
    end
end

%%
for nn = 1:nCells
    figure(100+nn);
    clf;
    set(gcf, 'Position',  [100, 100, 1600, 300]);
    set(gcf, 'Color', 'w', 'InvertHardcopy', 'off');
    imagesc(connMatrix{2, 1});
    colormap((gray));
    clim([0 1])
    axis equal tight;         % keep cells square
    
    % Draw grid lines
    hold on
    [nRows, nCols] = size(connMatrix{1, 1});
    for r = 0.5:1:nRows+0.5
        plot([0.5, nCols+0.5], [r r], 'k-', 'LineWidth', 0.75);
    end
    for c = 0.5:1:nCols+0.5
        plot([c c], [0.5, nRows+0.5], 'k-', 'LineWidth', 0.75);
    end
    hold off
    
    cbConnImg = colorbar;
    set(cbConnImg, 'FontSize', 9, 'Color', 'k');
    if isprop(cbConnImg, 'FontColor')
        cbConnImg.FontColor = 'k';
    end
    set(gca, 'TickLength', [0 0]); 
    set(gca, 'fontsize', 16);
    set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    xticks(1:nCols);
    yticks(1:nRows);
    xlabel('Presynaptic Target');
    % ylabel('Cell');
    title(sprintf('Connectivity according to Imaging (paired t-test p, BH-FDR q=%.3g across holos)', fdrQConn));
    set(get(gca, 'Title'), 'Color', 'k');
    
    % Response amplitude-based connectivity matrix
    figure(200+nn);
    clf;
    set(gcf, 'Position',  [100, 100, 1600, 300]);
    set(gcf, 'Color', 'w', 'InvertHardcopy', 'off');
    imagesc(respMatrix{1, 1});
    colormap default;
    axis equal tight;         % keep cells square
    
    % Draw grid lines
    hold on
    [nRows, nCols] = size(respMatrix{1, 1});
    for r = 0.5:1:nRows+0.5
        plot([0.5, nCols+0.5], [r r], 'k-', 'LineWidth', 0.75);
    end
    for c = 0.5:1:nCols+0.5
        plot([c c], [0.5, nRows+0.5], 'k-', 'LineWidth', 0.75);
    end
    hold off
    
    cbAmpImg = colorbar;
    set(cbAmpImg, 'FontSize', 9, 'Color', 'k');
    if isprop(cbAmpImg, 'FontColor')
        cbAmpImg.FontColor = 'k';
    end
    caxis([0 0.015])
    set(gca, 'TickLength', [0 0]); 
    set(gca, 'fontsize', 16);
    set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    xticks(1:nCols);
    yticks(1:nRows);
    % xlabel('Presynaptic Target');
    % ylabel('Cell');
    title('Postsynaptic Response Amplitude (dF/F)');
    set(get(gca, 'Title'), 'Color', 'k');
    
    if useEphysForConn
        figure(300+nn);
        clf;
        set(gcf, 'Position',  [100, 450, 1600, 180]);
        set(gcf, 'Color', 'w', 'InvertHardcopy', 'off');
        imagesc(connMatrixEphys{1, 1});
        colormap((gray));
        clim([0 1]);
        axis equal tight;
        hold on
        [nRowsE, nColsE] = size(connMatrixEphys{1, 1});
        for r = 0.5:1:nRowsE+0.5
            plot([0.5, nColsE+0.5], [r r], 'k-', 'LineWidth', 0.75);
        end
        for c = 0.5:1:nColsE+0.5
            plot([c c], [0.5, nRowsE+0.5], 'k-', 'LineWidth', 0.75);
        end
        hold off
        cbConnEph = colorbar;
        set(cbConnEph, 'FontSize', 9, 'Color', 'k');
        if isprop(cbConnEph, 'FontColor')
            cbConnEph.FontColor = 'k';
        end
        set(gca, 'TickLength', [0 0]);
        set(gca, 'fontsize', 16);
        set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
        xticks(1:nColsE);
        yticks(1:nRowsE);
        xlabel('Presynaptic Target');
        % ylabel('Patched Cell1');
        title(sprintf('Connectivity according to Electrophysiology (paired t-test p, BH-FDR q=%.3g across holos)', fdrQConn));
        set(get(gca, 'Title'), 'Color', 'k');
    
        % Ephys response amplitudes (same layout as fig. 102; one row = patched cell).
        figure(400+nn);
        clf;
        set(gcf, 'Position',  [100, 640, 1600, 180]);
        set(gcf, 'Color', 'w', 'InvertHardcopy', 'off');
        imagesc(respMatrixEphys{1, 1});
        colormap default;
        axis equal tight;
    
        hold on
        [nRowsE2, nColsE2] = size(respMatrixEphys{1, 1});
        for r = 0.5:1:nRowsE2+0.5
            plot([0.5, nColsE2+0.5], [r r], 'k-', 'LineWidth', 0.75);
        end
        for c = 0.5:1:nColsE2+0.5
            plot([c c], [0.5, nRowsE2+0.5], 'k-', 'LineWidth', 0.75);
        end
        hold off
    
        cbAmpEph = colorbar;
        set(cbAmpEph, 'FontSize', 9, 'Color', 'k');
        if isprop(cbAmpEph, 'FontColor')
            cbAmpEph.FontColor = 'k';
        end
        veFin = respMatrixEphys{1, 1}(isfinite(respMatrixEphys{1, 1}));
        if ~isempty(veFin)
            vLo = min(veFin);
            vHi = max(veFin);
            if vHi > vLo
                caxis([vLo vHi]);
            end
        end
        set(gca, 'TickLength', [0 0]);
        set(gca, 'fontsize', 16);
        set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
        xticks(1:nColsE2);
        yticks(1:nRowsE2);
        xlabel('Presynaptic Target');
        % ylabel('Ephys (patched cell)');
        title('Postsynaptic response amplitude (mV)');
        set(get(gca, 'Title'), 'Color', 'k');
    end
end

%% Trial raster viewer (choose condition/hologram)
% Visualize trial-wise traces for a selected holo as side-by-side rasters.
% Imaging uses one selected ROI (default: paired ROI if available).
showTrialRaster = true;
if showTrialRaster
    rasterCondToPlot = 1; % set [] to choose interactively in command window
    rasterHoloToPlot = []; % set [] to choose interactively in command window
    rasterCellToPlot = 1; % imaging ROI index; overridden by paired ROI when available

    if useEphysForConn && isfield(voltMapping.ephys, 'pairedRoiIndex') ...
            && ~isempty(voltMapping.ephys.pairedRoiIndex) && ~isnan(voltMapping.ephys.pairedRoiIndex)
        rasterCellToPlot = round(voltMapping.ephys.pairedRoiIndex);
    end
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
            imgTrials = voltMapping.(cellID{rasterCellToPlot}).exclFiltHoloSortedImagingAllTrials{rasterCondToPlot}{rasterHoloToPlot};
            ephTrials = [];
            if useEphysForConn && rasterCondToPlot <= numel(ephysTrialsByCond) ...
                    && rasterHoloToPlot <= numel(ephysTrialsByCond{rasterCondToPlot})
                ephTrials = ephysTrialsByCond{rasterCondToPlot}{rasterHoloToPlot};
            end

            if isempty(imgTrials) && isempty(ephTrials)
                warning('Trial raster skipped: no imaging/ephys trial data for cond %d holo %d.', ...
                    rasterCondToPlot, rasterHoloToPlot);
            else
                nTrialsImg = size(imgTrials, 2);
                nTrialsEph = size(ephTrials, 2);
                nTrialCommon = min(nTrialsImg, nTrialsEph);
                if nTrialCommon > 0
                    imgTrials = imgTrials(:, 1:nTrialCommon);
                    ephTrials = ephTrials(:, 1:nTrialCommon);

                    % Keep only trials that are finite in both modalities.
                    validImgCols = all(~isnan(imgTrials), 1);
                    validEphCols = all(~isnan(ephTrials), 1);
                    validBothCols = validImgCols & validEphCols;
                    imgTrials = imgTrials(:, validBothCols);
                    ephTrials = ephTrials(:, validBothCols);
                end

                figure(105);
                clf;
                set(gcf, 'Position', [100, 100, 520, 520]);
                set(gcf, 'Color', 'w', 'InvertHardcopy', 'off');

                subplot(1, 2, 1);
                if ~isempty(ephTrials)
                    nSEph = size(ephTrials, 1);
                    nTrialEph = size(ephTrials, 2);
                    tMsEph = (0:nSEph-1) / Fs * 1000;
                    firstPulseMsEph = (stimStartEphys - 1) / Fs * 1000;
                    imagesc(tMsEph, 1:nTrialEph, ephTrials');
                    clim([0 20]);
                    hold on;
                    xline(firstPulseMsEph, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
                    hold off;
                    xlabel('Time (ms)');
                    ylabel('Trial');
                    title(sprintf('Ephys raster (cond %d, holo %d)', rasterCondToPlot, rasterHoloToPlot));
                    axEph = gca;
                    axEph.Title.Position(2) = axEph.Title.Position(2) - 1.5;
                    set(axEph, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontSize', 12);
                    axEph.Title.Color = 'k';
                    colormap(gca, 'parula');
                    cbEph = colorbar;
                    set(cbEph, 'FontSize', 9, 'Color', 'k');
                    if isprop(cbEph, 'FontColor')
                        cbEph.FontColor = 'k';
                    end
                else
                    axis off;
                    text(0.5, 0.5, 'No ephys trial data', 'HorizontalAlignment', 'center');
                end

                subplot(1, 2, 2);
                if ~isempty(imgTrials)
                    nSImg = size(imgTrials, 1);
                    nTrialImg = size(imgTrials, 2);
                    tMsImg = (0:nSImg-1) / imagingFreq * 1000;
                    firstPulseMsImg = (stimStart - 1) / imagingFreq * 1000;
                    imagesc(tMsImg, 1:nTrialImg, imgTrials');
                    clim([0 0.15]);
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
                else
                    axis off;
                    text(0.5, 0.5, 'No imaging trial data', 'HorizontalAlignment', 'center');
                end
                
            end
        end
    end
end

%% Pearson correlation coefficient (within session)
nIter = 1000; % number of times random halves will be compared. Choose from 100~500
CCwsAllCells = cell(nCells, 1);
for nn = 1:nCells
    CCwsAllConds = cell(nConds, 1);
    clear exclFiltHoloSortedImagingAllTrials
    exclFiltHoloSortedImagingAllTrials = voltMapping.(cellID{nn}).exclFiltHoloSortedImagingAllTrials;
    for cc = 1:nConds
        CCwsAllHolos = nan(nHolos(cc), 1);
        for hh = 1:nHolos(cc)
            testData = exclFiltHoloSortedImagingAllTrials{cc}{hh};
            if isempty(testData)
                continue;
            end

            validCols = all(~isnan(testData), 1);
            validData = testData(:, validCols);
            nTrials = size(validData, 2);
            halfTrials = floor(nTrials/2);
            if halfTrials < 1
                continue;
            end

            CCAllIter = nan(nIter, 1);
            for ii = 1:nIter
                randIdx = randperm(nTrials);
                firstHalf = mean(validData(:, randIdx(1:halfTrials)), 2, "omitnan");
                secondHalf = mean(validData(:, randIdx(halfTrials+1:2*halfTrials)), 2, "omitnan");
                CCThisIter = corr(firstHalf, secondHalf, "Rows", "complete");
                if isnan(CCThisIter)
                    CCThisIter = 0;
            end
                CCAllIter(ii) = max(CCThisIter, 0);
            end

            CCwsAllHolos(hh) = mean(CCAllIter, 'omitnan');
        end
        CCwsAllConds{cc} = CCwsAllHolos;
    end
    CCwsAllCells{nn} = CCwsAllConds;
end

%%
% CCwsAllSessions = struct()
CCwsAllSessions.D6 = CCwsAllCells;

%% Save Analysis Results
expID = 'MS25_9';
voltMapping.mouseID = ['voltMapping_Analysis_', expID, '_MultiCellAnalysis', '_CC'];
% directory = '/Volumes/ExData2/Voltage Imaging/VoltMapping/Analysis Results';
directory = 'E:\Voltage Imaging\VoltMapping\Analysis Results';
fileName = [num2str(voltMapping.mouseID), '.mat'];
save(fullfile(directory, fileName), ...
    'CCwsAllSessions', ...
    'connMatrix', 'respMatrix', 'respMatrixConnectedCells', 'percentConn', ...
    'connMatrixEphys', 'respMatrixEphys', 'respMatrixEphysConnected', 'percentConnEphys', ...
    'connAgreement', 'alphaConn', ...
    '-v7.3');

TimeNow = datetime('now', 'Format', 'dd-MMM-yyyy HH:mm:ss');
disp(['finished saving at: ' char(TimeNow)])

%% Plot CCws distributions across multiple sessions
% Concatenate sessions of one cell into one array
CCwsThisCell = [CCwsAllSessions.D1{4, 1}{1, 1};...
    CCwsAllSessions.D2{4, 1}{2, 1};...
    CCwsAllSessions.D3{4, 1}{1, 1};...
    CCwsAllSessions.D4{4, 1}{1, 1};...
    CCwsAllSessions.D5{4, 1}{1, 1};...
    CCwsAllSessions.D6{4, 1}{1, 1}];
labels = [ones(size(CCwsAllSessions.D1{4, 1}{1, 1}));
          2*ones(size(CCwsAllSessions.D2{4, 1}{1, 1}));
          3*ones(size(CCwsAllSessions.D3{4, 1}{1, 1}));
          4*ones(size(CCwsAllSessions.D4{4, 1}{1, 1}));
          5*ones(size(CCwsAllSessions.D5{4, 1}{1, 1}));
          6*ones(size(CCwsAllSessions.D6{4, 1}{1, 1}))];

% Create box plot
figure(103); 
set(gcf, 'Position',  [100, 100, 560, 420])
clf
hold on;
boxplot(CCwsThisCell, labels, 'Colors',[0 0 0], 'Symbol','');
swarmchart(labels, CCwsThisCell, 12, 'k','filled','MarkerFaceAlpha',0.3);
set(gca,'XTick',1:6,'XTickLabel',{'Day 1','Day 4','Day 7','Day 10','Day 14', 'Day 17'});
set(gca, 'fontsize', 12);
ylabel('CC within-session');
title('Cell4: Within-session reliability across days');
ylim([0 0.6]); 
box on; 
% grid on;

% Optionally overlay mean � SEM for each session
CCwsThisCell2 = [CCwsAllSessions.D1{4, 1}{1, 1},...
    CCwsAllSessions.D2{4, 1}{2, 1},...
    CCwsAllSessions.D3{4, 1}{1, 1},...
    CCwsAllSessions.D4{4, 1}{1, 1},...
    CCwsAllSessions.D5{4, 1}{1, 1},...
    CCwsAllSessions.D6{4, 1}{1, 1}];

for ii = 1:numel(CCwsThisCell2)
    m = mean(CCwsThisCell2(:, ii),'omitnan');
    s = std(CCwsThisCell2(:, ii),'omitnan')/sqrt(numel(CCwsThisCell2(:, ii)));
    plot(ii, m, 'ko','MarkerFaceColor','k');
    line([ii ii], [m-s, m+s], 'Color','r','LineWidth',1.5);
end
hold off;

%% Stats (Kurskal Wallis)
CCwsThisCell = [CCwsAllSessions.D1{1, 1}{1, 1};...
    CCwsAllSessions.D2{1, 1}{2, 1};...
    CCwsAllSessions.D3{1, 1}{1, 1};...
    CCwsAllSessions.D4{1, 1}{1, 1};...
    CCwsAllSessions.D5{1, 1}{1, 1};...
    CCwsAllSessions.D6{1, 1}{1, 1}];
labels = [ones(size(CCwsAllSessions.D1{1, 1}{1, 1}));
          2*ones(size(CCwsAllSessions.D2{1, 1}{2, 1}));
          3*ones(size(CCwsAllSessions.D3{1, 1}{1, 1}));
          4*ones(size(CCwsAllSessions.D4{1, 1}{1, 1}));
          5*ones(size(CCwsAllSessions.D5{1, 1}{1, 1}));
          6*ones(size(CCwsAllSessions.D6{1, 1}{1, 1}))];

[p,tbl,stats] = kruskalwallis(CCwsThisCell, labels, 'off');

fprintf('Kruskal�Wallis p = %.4f\n', p);
if p > 0.05
    disp(' No significant change in CCws across sessions: within-session reliability is stable.');
else
    disp(' Significant difference: investigate experimental drift.');
end