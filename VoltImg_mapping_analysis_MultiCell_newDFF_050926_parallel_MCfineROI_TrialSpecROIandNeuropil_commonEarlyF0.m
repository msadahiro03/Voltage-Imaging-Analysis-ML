%% Voltage Imaging Mapping Analysis 050926 newDFFcalc and MCfineROI (common early-trial F0) — parallel MC pass 2 + parallel dF/F trials
% Parallel variant: parfor over NoRMCorre pass 2 (after global template) and parfor over imaging
% trials in the dF/F section. Helpers: VoltImg_mapping_parallel_mcPass2OneTrial,
% VoltImg_mapping_parallel_dffOneTrial_commonEarlyF0.
%
% F0 for all holograms in a trial/cell is shared: within the last commonEarlyF0BaselineMs
% before the first pulse of the first hologram, a rolling window (commonEarlyF0RollingWinMs)
% with minimum variance defines F0 (mean of that sub-window), separately for raw and filtered traces.
%
% Copy of TrialSpecROIandNeuropil pipeline with optional pre-NoRMCorre laser row cleaning
% (VoltImg_applyLaserRowArtifactToStack, VoltImg_mapping_removeArtifact_v2).
%
% 1) Per-trial NoRMCorre motion correction of raw image stacks, saved to a
%    separate folder.
% 2) maxDvStack: per-trial mean planes from MC stacks — NaN for ephys-excluded trials and for trials
%    not drawn into the random 50%% eligible-trial subsample; each contributing plane uses only the
%    first max(1,floor(imagingFreq*4)) frames. Grand mean uses 'omitnan'.
% 3) Hand-drawn rough ROIs on mean of motion-corrected trials.
% 4) Per-trial "fine" ROIs computed inside each rough ROI using the same
%    gaussian + fibermetric procedure as the original script, but applied
%    to that trial's mean image.
% 5) F, F0, dF, dF/F0 from motion-corrected stacks and per-trial fine ROIs;
%    F0 is trial-common (early prestim before first holo; min-var rolling window).

%%
clear all
%close all

%% Load files and setup
% Step 1: Read the ephys file
% ephysFilePath = char(uigetdir('E:\Voltage Imaging\VoltMapping\Ephys Data')); % Select and set root folder where all experiments with cells you want to analyze are located
% ephysFilePath = char(uigetdir('/Volumes/ExData2/Voltage Imaging/VoltMapping/Ephys Data')); % Select and set root folder where all experiments with cells you want to analyze are located
ephysFilePath = char(uigetdir('/Volumes/phoenixinthesky/Masato/Voltage Imaging Data_Phoenix/SliceMapping/ASAP7y Original Slice Experiment/DAQ Ephys Data/'));

ephysFileDir = dir(ephysFilePath);
load([ephysFileDir(end).folder, '/', ephysFileDir(end).name]);
disp(ephysFileDir(end).name);
 
% Step 2: Identify the folder containing the imaging files correspdonding to the ephys file
% ImgsFilePath = char(uigetdir('E:\Voltage Imaging\VoltMapping\Imaging Data')); % Select and set root folder where all experiments with cells you want to analyze are located
% ImgsFilePath = char(uigetdir('/Volumes/ExData2/Voltage Imaging/VoltMapping/Imaging Data')); % Select and set root folder where all experiments with cells you want to analyze are located
ImgsFilePath = char(uigetdir('/Users/masatosadahiro/Documents/Data/Voltage Imaging/Voltage Imaging/Slice Mapping/'));

ImgfolderContents = dir(ImgsFilePath);
disp(ImgfolderContents(end).name);

% Step 2a: Avoid hidden files and non image files
fileNames = [];
fileType = '.tif';
for ii = 1:length(ImgfolderContents)
    if ~ImgfolderContents(ii).isdir && ~startsWith(ImgfolderContents(ii).name, '.') && endsWith(ImgfolderContents(ii).name, fileType)
        fileNames{ii, 1} = ImgfolderContents(ii).name;
    end
end
imagesIndex = find(~cellfun(@isempty, fileNames));

% Step 3: NoRMCorre path
normcorrePath = 'C:\Users\lamia\OneDrive\Documents\MATLAB\NoRMCorre-master';
addpath(normcorrePath);

% --- Laser row artifact (applied to raw stacks before global template + NoRMCorre) ---
useLaserRowArtifactFilter = false;
laserArtifactGateColFirst = 130;
laserArtifactGateColLast  = 382;
laserArtifactThreshMode   = 'mad';       % 'fixed' | 'mad' | 'percentile'
laserArtifactThreshParam  = 5;           % fixed: variance cutoff; mad: k (lower = cuts more aggressively); percentile: P in [0 100]
laserArtifactMcMode       = 'fill_for_mc'; % 'fill_for_mc' (finite, MC-safe) | 'nan' (experimental)
mcUseGateColumnsOnly      = true;       % If true, estimate shifts from gate columns only, then apply shifts to full frames.
laserArtifactMcSecondSweepForDff = false; % If true: after MC, save HxT bad-row mask; dF/F omits those ROI/neuropil pixels per frame.

% --- Parallel execution (Parallel Computing Toolbox) ---
useParallelMcPass2 = false;   % parfor NoRMCorre pass 2 after global template is built
useParallelDffTrials = false; % parfor over trials in dF/F section

% Step 4: Setup struct
voltMapping = ExpStruct;

if exist('ExpStruct2', 'var') % If the experiment has a second patch electrode
    voltMapping.ExpStruct2 = ExpStruct2;
end

% Stimulation properties
imagingFreq = voltMapping.sampleFreq;
Fs = voltMapping.daqParams.Fs;
trialTime = voltMapping.daqParams.maxSweepLengthSec;
nTrials = length(voltMapping.trialCond);
powers = voltMapping.outParams.power;
nConds = length(voltMapping.outParams.sequence);
nHolos = voltMapping.holoStimParams.nHolos;
nHolos(1) = max(nHolos); % Hack for 0 holos conditions because of 0mW trials involved

% Stimulation properties - assume all of these properties are the same across all conditions/powers
pulseDurs = unique(voltMapping.outParams.pulseDur); % Assume same pulse duration or stim rate for all conditions
    pulseDurs = nonzeros(pulseDurs); % Hack for 0 pulse duration conditions because of 0mW trials involved
nPulses = unique(voltMapping.outParams.nPulses); % Assume every target in every condition receives same # of pulses (remove the no '0 pulses' for the 0mW conditions)
    nPulses = nonzeros(nPulses); % Hack for 0mW condition involved
ipi = unique(voltMapping.outParams.ipi); % Assume same ipi or stim rate for all conditions
    ipi = nonzeros(ipi); % Hack for 0mW condition involved
nextHoloDelay = unique(voltMapping.holoStimParams.nextHoloDelay); % Assume same delay between holos for all conditions
    nextHoloDelay = nonzeros(nextHoloDelay); % Hack for 0mW condition involved
startTime = (voltMapping.holoStimParams.startTime)/1000;

% Step 5: Identify GEVI type
UpOrDown = input('1 for upward GEVI, 2 for downward GEVI: ', 's');
ePhysAvail = input('1 if ephys readout avail, 2 if none: ');

% Auto-detect whether raw TIFF stacks are single-channel or 2-color interleaved.
% Assumes all stacks in this folder were acquired with identical settings,
% so probing one stack is sufficient.
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
    % Too few pages to robustly detect interleaving; default to single-channel.
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

    % Metric 1: odd/even mean-image correlation (lower suggests different channels)
    rImg = corrcoef(double(oddMeanImg(:)), double(evenMeanImg(:)));
    if numel(rImg) >= 4
        oddEvenImgCorr = rImg(1, 2);
    else
        oddEvenImgCorr = 1;
    end

    % Metric 2: alternating-frame intensity pattern
    if numel(frameMeans) >= 3
        lag1 = corrcoef(frameMeans(1:end-1), frameMeans(2:end));
        lag2 = corrcoef(frameMeans(1:end-2), frameMeans(3:end));
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
        sameChanDiff = mean(abs(frameMeans(3:end) - frameMeans(1:end-2)));
    else
        lag1Corr = 0;
        lag2Corr = 0;
        altStepDiff = 0;
        sameChanDiff = 0;
    end

    % Decide 2-channel interleaving only when multiple signs agree.
    isInterleaved = (oddEvenImgCorr < 0.90) && ...
                    ((lag2Corr > lag1Corr + 0.10) || (altStepDiff > 1.15 * max(sameChanDiff, eps)));
    if isInterleaved
        rawImgNChannels = 2;
    else
        rawImgNChannels = 1;
    end

    disp(['Auto-detected raw TIFF channel mode from ', ImgfolderContents(imagesIndex(1)).name, ...
          ': rawImgNChannels = ', num2str(rawImgNChannels), ...
          ' (odd-even mean image corr = ', num2str(oddEvenImgCorr, '%.3f'), ...
          ', lag1 = ', num2str(lag1Corr, '%.3f'), ', lag2 = ', num2str(lag2Corr, '%.3f'), ').']);
end

% Step 6: Save directory (later overwritten near save step)
mouseID = ExpStruct.mouseID;
directory = 'D:\Data\Voltage Imaging\voltMapping\Analysis Results\';
fileName = ['voltMapping ', num2str(ExpStruct.mouseID), '.mat'];

voltMapping.imagesIndex   = imagesIndex;
voltMapping.imagingFreq   = imagingFreq;
voltMapping.UpOrDown      = UpOrDown;
voltMapping.rawImgNChannels = rawImgNChannels;
voltMapping.ephysFilePath = ephysFilePath;
voltMapping.ImgsFilePath  = ImgsFilePath;
voltMapping.nConds        = nConds;
voltMapping.nHolos        = nHolos;
voltMapping.pulseDurs     = pulseDurs;
voltMapping.nPulses       = nPulses;
voltMapping.ipi           = ipi;
voltMapping.nextHoloDelay = nextHoloDelay;

%% Trial exclusion and baselining
vThreshold = -60; %set threshold Vm here
mappingInputs = ExpStruct.inputs;

baseVoltAllTrials = [];
excludeTrials = [];
mappingInputsBaselined =[];
for tt = 1:nTrials
    baseVolt = mean(mappingInputs{tt}(1:startTime*Fs));
    baseVoltAllTrials = [baseVoltAllTrials; baseVolt];
    mappingInputsBaselined{tt, 1} = mappingInputs{tt} - baseVolt;
    
    if baseVoltAllTrials(tt) > vThreshold
        excludeTrials = [excludeTrials, tt];
    end
end

if ePhysAvail == 2
    excludeTrials = [];
end

voltMapping.excludeTrials                = excludeTrials;
voltMapping.ephys.baseVoltAllTrials      = baseVoltAllTrials;
voltMapping.ephys.mappingInputsBaselined = mappingInputsBaselined;

%% Break apart each continuous ephys trace into stim windows and rearrange according to hologram sequence and conditions
cutOffFreq = 480;   % Cutoff frequency for Butterworth filter
[blp, alp] = butter(4, cutOffFreq/(Fs/2), 'low');

holoSeqIndex = cell(nConds, 1);
holoSortedDataAllTrials = cell(nConds, 1);
for cc = 1:nConds
    holoSortedDataAllTrials{cc} = cell(nHolos(cc), 1);
end

postStimWindow = 50; % time(ms) to add to stim window after final pulse + ipi
preStimWindow = nextHoloDelay - postStimWindow; % time(ms) window before first pulse 

nPulseCoords = [];
for pp = 1:nPulses
    nPulseCoords = [nPulseCoords, ((pp-1)*ipi/1000*Fs) + preStimWindow/1000*Fs];
end

nPulseCoordsImaging = [];
for pp = 1:nPulses
    nPulseCoordsImaging = [nPulseCoordsImaging, ((pp-1)*ipi/1000*imagingFreq) + preStimWindow/1000*imagingFreq];
end

% Prestim baseline for holo-sorted sweeps: samples strictly before first pulse
% (same idea as imaging: 1:(preStimWindow*freq)-1), with ms edge trim, linear
% detrend on prestim only, then median offset (robust; aligns trials with drift).
nFirstPulseSample = round(nPulseCoords(1));
iLastPreStim = max(1, nFirstPulseSample - 1);
preStimBaselineTrimMs = 2;
edgeSampPreBaseline = round(preStimBaselineTrimMs/1000*Fs);
preStimBaselineIdx1 = min(iLastPreStim, 1 + edgeSampPreBaseline);
preStimBaselineIdx2 = max(preStimBaselineIdx1, iLastPreStim - edgeSampPreBaseline);
if preStimBaselineIdx2 < preStimBaselineIdx1
    preStimBaselineIdx1 = 1;
    preStimBaselineIdx2 = iLastPreStim;
end
preStimBaselineIdx = preStimBaselineIdx1:preStimBaselineIdx2;

% Unified holo ephys length: use the shortest nominal window across holos/conds
% (floor/ceil on fst makes lengths differ by ~1; truncating longer holos avoids
% NaN tail padding). Then cap by the smallest overlap with the recorded trace
% across non-excluded trials so the common length is always fully sampled.
ephysHoloSweepLenNominal = inf;
for cc = 1:nConds
    fstVec = voltMapping.outParams.firstStimTimes{cc};
    if isempty(fstVec)
        fstVec = voltMapping.outParams.firstStimTimes{1, 2};
    end
    for hh = 1:nHolos(cc)
        fst = fstVec(min(hh, numel(fstVec)));
        iLo = floor((fst - preStimWindow/1000) * Fs);
        iHi = ceil((fst - preStimWindow/1000) * Fs) + ...
            ceil((ipi*nPulses + preStimWindow + postStimWindow) / 1000 * Fs);
        ephysHoloSweepLenNominal = min(ephysHoloSweepLenNominal, iHi - iLo + 1);
    end
end
if ~(ephysHoloSweepLenNominal < inf) || ephysHoloSweepLenNominal < 1
    ephysHoloSweepLenNominal = ceil((ipi*nPulses + preStimWindow + postStimWindow) / 1000 * Fs) + 2;
end

minAvailAcrossTrials = inf;
for tt = 1:nTrials
    if ismember(tt, excludeTrials)
        continue
    end
    nSweep = numel(mappingInputsBaselined{tt});
    cc = voltMapping.trialCond(tt, 1);
    fstVec = voltMapping.outParams.firstStimTimes{cc};
    if isempty(fstVec)
        fstVec = voltMapping.outParams.firstStimTimes{1, 2};
    end
    for hh = 1:nHolos(cc)
        fst = fstVec(min(hh, numel(fstVec)));
        iEphysLo = floor((fst - preStimWindow/1000) * Fs);
        iEphysHi = ceil((fst - preStimWindow/1000) * Fs) + ...
            ceil((ipi*nPulses + preStimWindow + postStimWindow) / 1000 * Fs);
        srcLo = max(1, iEphysLo);
        srcHi = min(nSweep, iEphysHi);
        if srcLo <= srcHi
            minAvailAcrossTrials = min(minAvailAcrossTrials, srcHi - srcLo + 1);
        end
    end
end

ephysHoloSweepLen = ephysHoloSweepLenNominal;
if minAvailAcrossTrials < inf && minAvailAcrossTrials >= 1
    ephysHoloSweepLen = min(ephysHoloSweepLenNominal, minAvailAcrossTrials);
end
if ephysHoloSweepLen < 1
    ephysHoloSweepLen = 1;
end

% Create dummy sequence for 0mV trials (if any)
if any(voltMapping.outParams.power == 0)
    zeroDummySequence = voltMapping.outParams.sequence{1, 2};
end

condSortedInputs = cell(nConds, 1);
for tt = 1:nTrials
    
    if isempty(voltMapping.outParams.sequenceThisTrial{tt})
        voltMapping.outParams.sequenceThisTrial{tt} = zeroDummySequence;
    end

    holoSeqThisTrial = unique(voltMapping.outParams.sequenceThisTrial{tt}, 'stable') - min(unique(voltMapping.outParams.sequenceThisTrial{tt}, 'stable')) + 1;
    holoSeqIndex{voltMapping.trialCond(tt, 1)} = [holoSeqIndex{voltMapping.trialCond(tt, 1)}, holoSeqThisTrial'];

    sortingIndex = holoSeqIndex{voltMapping.trialCond(tt, 1)}(:, end) - min(holoSeqIndex{voltMapping.trialCond(tt, 1)})+1;
    
    condSortedInputs{voltMapping.trialCond(tt, 1)} = [condSortedInputs{voltMapping.trialCond(tt, 1)}, mappingInputsBaselined{tt}];
    
    if ~ismember(tt, excludeTrials)
        sweepThisTrial = filtfilt(blp, alp, mappingInputsBaselined{tt});
    end
    
    for hh = 1:nHolos(voltMapping.trialCond(tt, 1))
        if isempty(voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)})
            voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)} = voltMapping.outParams.firstStimTimes{1, 2};
        end

        % Integer sample indices (same window construction as imaging iHoloLo:iHoloHi)
        fst = voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh);
        iEphysLo = floor((fst - preStimWindow/1000) * Fs);
        if ismember(tt, excludeTrials)
            thisHoloSweep = nan(ephysHoloSweepLen, 1);
        else
            nSweep = numel(sweepThisTrial);
            thisHoloSweep = nan(ephysHoloSweepLen, 1);
            for ii = 1:ephysHoloSweepLen
                srcIdx = iEphysLo + ii - 1;
                if srcIdx >= 1 && srcIdx <= nSweep
                    thisHoloSweep(ii) = sweepThisTrial(srcIdx);
                end
            end
            v = thisHoloSweep(:);
            nS = numel(v);
            % Prestim within this holo chunk (trial-relative indices are wrong here; mirror imaging 1:preStim-1)
            nPreSamp = max(1, round(preStimWindow/1000 * Fs) - 1);
            preStimBaselineIdxLocal = 1:min(nPreSamp, nS);
            tRel = ((1:nS)' - mean(preStimBaselineIdxLocal)) / Fs;
            if numel(preStimBaselineIdxLocal) >= 3
                pLin = polyfit(tRel(preStimBaselineIdxLocal), v(preStimBaselineIdxLocal), 1);
                v = v - polyval(pLin, tRel);
            end
            v = v - median(v(preStimBaselineIdxLocal), 'omitnan');
            thisHoloSweep = reshape(v, size(thisHoloSweep));
        end

        holoSortedDataAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, end), 1} = [holoSortedDataAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, end), 1}, thisHoloSweep];
    end
end

% Show average traces for each hologram
holoSortedDataMean = cell(nConds, 1);
for cc = 1:nConds
    for hh = 1:nHolos(cc)
        holoSortedDataMean{cc}(:, hh) = nanmean(holoSortedDataAllTrials{cc}{hh}, 2);
    end
end

CIephysAllConds = cell(nConds, 1);
for cc = 1:nConds 
    for hh = 1:nHolos(cc)
        confidence_level = 0.95;
        means = nanmean(holoSortedDataAllTrials{cc}{hh, 1}, 2);
        std_errors = std(holoSortedDataAllTrials{cc}{hh, 1}, 0, 2, "omitnan") / sqrt(size(holoSortedDataAllTrials{cc}{hh, 1}, 2));
    
        t_score = tinv((1 + confidence_level) / 2, size(holoSortedDataAllTrials{cc}{hh, 1}, 2) - 1);
        margin_of_error = t_score * std_errors;
        lower_bounds = means - margin_of_error;
        upper_bounds = means + margin_of_error;
        if UpOrDown == '2'
            CIephysAllConds{cc}{hh, 1} = [lower_bounds, upper_bounds];
        elseif UpOrDown =='1'
            CIephysAllConds{cc}{hh, 1} = [-lower_bounds, -upper_bounds];
        end
    end
end

voltMapping.ephys.holoSeqIndex            = holoSeqIndex; 
voltMapping.ephys.holoSortedDataAllTrials = holoSortedDataAllTrials;
voltMapping.ephys.nPulseCoords            = nPulseCoords;
voltMapping.ephys.CIephysAllConds         = CIephysAllConds;
voltMapping.ephys.holoSortedDataMean      = holoSortedDataMean;
voltMapping.ephys.preStimWindow           = preStimWindow;
voltMapping.ephys.postStimWindow          = postStimWindow;
voltMapping.ephys.condSortedInputs        = condSortedInputs;
voltMapping.ephys.preStimBaselineIdx      = preStimBaselineIdx;
voltMapping.ephys.ephysHoloSweepLen       = ephysHoloSweepLen;
voltMapping.ephys.ephysHoloSweepLenNominal = ephysHoloSweepLenNominal;
voltMapping.ephys.minAvailAcrossTrials    = minAvailAcrossTrials;
voltMapping.ephys.preStimBaselineTrimMs   = preStimBaselineTrimMs;
voltMapping.ephys.nFirstPulseSample       = nFirstPulseSample;
voltMapping.ephys.iLastPreStim            = iLastPreStim;

%% Motion-correct all imaging trials with NoRMCorre and build maxDvStack
input('Running NoRMCorre on all imaging trials, saving motion-corrected stacks, and building maxDvStack (continue or ctrl+c to stop!)');

% Setup save directory for motion-corrected images
expID = num2str(mouseID);
voltMapping.mouseID = ['voltMapping_Analysis_', expID, '_MultiCellAnalysis_MCfineROI_laserRowArtifact_parallel'];
savePath = '/Volumes/X10 Pro/MC Imaging Data Parallel';
saveDirectory = fullfile(savePath, num2str(voltMapping.mouseID));
if ~exist(saveDirectory, 'dir')
    mkdir(saveDirectory);
end

% Create subfolder for motion-corrected TIFFs
mcTiffFolder = fullfile(saveDirectory, 'Motion_Corrected_Tiffs');
if ~exist(mcTiffFolder, 'dir')
    mkdir(mcTiffFolder);
end

% Random 50%% of eligible (non-excluded) imaging trials + first 4 s of frames for maxDv reference planes
eligibleTrialTT = setdiff(1:length(imagesIndex), excludeTrials(:).');
[maxDvTrialMask, maxDvFrameCap] = VoltImg_mapping_maxDvStackSamplingPlan(length(imagesIndex), imagingFreq, eligibleTrialTT);

% Preallocate maxDvStack using first TIFF dimensions
firstImgPath = fullfile(ImgfolderContents(imagesIndex(1)).folder, ImgfolderContents(imagesIndex(1)).name);
infoFirst = imfinfo(firstImgPath);
% Third dim is trial index; excluded / non-subsampled trials are NaN (mean(...,'omitnan')).
maxDvStack = nan(infoFirst(1).Height, infoFirst(1).Width, length(imagesIndex));

% Pass 1: build one global template from trial mean images so all trials
% are registered into the same reference frame.
if mcUseGateColumnsOnly
    gateColFirstMc = max(1, min(laserArtifactGateColFirst, laserArtifactGateColLast));
    gateColLastMc = min(infoFirst(1).Width, max(laserArtifactGateColFirst, laserArtifactGateColLast));
    if gateColLastMc < gateColFirstMc
        error('Invalid gate column range for mcUseGateColumnsOnly.');
    end
    globalTemplateAccum = zeros(infoFirst(1).Height, gateColLastMc-gateColFirstMc+1, 'single');
else
    globalTemplateAccum = zeros(infoFirst(1).Height, infoFirst(1).Width, 'single');
end
nTemplateTrials = 0;
for tt = 1:length(imagesIndex)
    disp(['Global template registering trial ', num2str(tt), ' of ', num2str(length(imagesIndex))]);
    % Trial index tt must match excludeTrials (ephys trial 1:nTrials), same as dF/F section.
    if ismember(tt, excludeTrials)
        continue
    end

    currImgPath = fullfile(ImgfolderContents(imagesIndex(tt)).folder, ...
                           ImgfolderContents(imagesIndex(tt)).name);
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
        % Interleaved two-color TIFF: green pages are 1, 3, 5, ...
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

    if useLaserRowArtifactFilter
        imageStack = VoltImg_applyLaserRowArtifactToStack(imageStack, laserArtifactGateColFirst, ...
            laserArtifactGateColLast, laserArtifactThreshMode, laserArtifactThreshParam, laserArtifactMcMode);
    end

    if mcUseGateColumnsOnly
        gateColFirstTrial = max(1, min(laserArtifactGateColFirst, laserArtifactGateColLast));
        gateColLastTrial = min(size(imageStack,2), max(laserArtifactGateColFirst, laserArtifactGateColLast));
        globalTemplateAccum = globalTemplateAccum + mean(single(imageStack(:, gateColFirstTrial:gateColLastTrial, :)), 3);
    else
        globalTemplateAccum = globalTemplateAccum + mean(single(imageStack), 3);
    end
    nTemplateTrials = nTemplateTrials + 1;
end

if nTemplateTrials == 0
    error('No non-excluded trials available to build a global motion-correction template.');
end
globalTemplate = globalTemplateAccum ./ nTemplateTrials;

% ---- Crash-safe checkpoint: globalTemplate built ----
% Overwritten every run; allows resuming motion-correction if MATLAB
% crashes during the second pass.
globalTemplateCheckpointFile = fullfile(tempdir, 'VoltImg_mapping_MCfineROI_globalTemplate_progress.mat');
checkpointVars_globalTemplate = {'globalTemplate', 'nTemplateTrials', 'imagesIndex', 'excludeTrials'};
try
    save(globalTemplateCheckpointFile, checkpointVars_globalTemplate{:}, '-v7.3');
catch
    save(globalTemplateCheckpointFile, checkpointVars_globalTemplate{:});
end
disp(['Checkpoint saved (global template ready): ', globalTemplateCheckpointFile]);

disp(['Checkpoint saved (global template ready): ', globalTemplateCheckpointFile]);

if (useParallelMcPass2 || useParallelDffTrials) && isempty(gcp('nocreate'))
    try
        parpool;
    catch ME
        warning('VoltImg:parpool', 'Could not start parallel pool (%s). parfor may run on the client only.', ME.message);
    end
end

mcPass2Ctx = struct();
mcPass2Ctx.ImgfolderContents = ImgfolderContents;
mcPass2Ctx.imagesIndex = imagesIndex;
mcPass2Ctx.rawImgNChannels = rawImgNChannels;
mcPass2Ctx.useLaserRowArtifactFilter = useLaserRowArtifactFilter;
mcPass2Ctx.laserArtifactGateColFirst = laserArtifactGateColFirst;
mcPass2Ctx.laserArtifactGateColLast = laserArtifactGateColLast;
mcPass2Ctx.laserArtifactThreshMode = laserArtifactThreshMode;
mcPass2Ctx.laserArtifactThreshParam = laserArtifactThreshParam;
mcPass2Ctx.laserArtifactMcMode = laserArtifactMcMode;
mcPass2Ctx.mcUseGateColumnsOnly = mcUseGateColumnsOnly;
mcPass2Ctx.globalTemplate = globalTemplate;
mcPass2Ctx.excludeTrials = excludeTrials;
mcPass2Ctx.mcTiffFolder = mcTiffFolder;
mcPass2Ctx.laserArtifactMcSecondSweepForDff = laserArtifactMcSecondSweepForDff;
mcPass2Ctx.maxDvFrameCap = maxDvFrameCap;
mcPass2Ctx.maxDvTrialMask = maxDvTrialMask;

if useParallelMcPass2
    mcPass2ProgressQ = parallel.pool.DataQueue;
    nMcPass2Trials = length(imagesIndex);
    afterEach(mcPass2ProgressQ, @(ttDone) fprintf(1, 'MC pass 2 completed trial %d of %d\n', ttDone, nMcPass2Trials));
    parfor tt = 1:nMcPass2Trials
        maxDvStack(:, :, tt) = VoltImg_mapping_parallel_mcPass2OneTrial(tt, mcPass2Ctx);
        send(mcPass2ProgressQ, tt);
    end
else
    for tt = 1:length(imagesIndex)
        disp(['Motion correcting trial ', num2str(tt), ' of ', num2str(length(imagesIndex))]);
        maxDvStack(:, :, tt) = VoltImg_mapping_parallel_mcPass2OneTrial(tt, mcPass2Ctx);
    end
end

meanMaxDvStack = mean(maxDvStack, 3, 'omitnan'); % Grand mean of motion-corrected trial means (non-excluded only)
meanFluorMaxDvStack = meanMaxDvStack;

voltMapping.mcTiffFolder           = mcTiffFolder;
voltMapping.maxDvStack         = maxDvStack;
voltMapping.meanFluorMaxDvStack = meanFluorMaxDvStack;
voltMapping.globalMcTemplate = globalTemplate;
voltMapping.maxDvTrialMask = maxDvTrialMask;
voltMapping.maxDvFrameCap = maxDvFrameCap;
voltMapping.laserRowArtifact.useFilter = useLaserRowArtifactFilter;
voltMapping.laserRowArtifact.gateColFirst = laserArtifactGateColFirst;
voltMapping.laserRowArtifact.gateColLast = laserArtifactGateColLast;
voltMapping.laserRowArtifact.threshMode = laserArtifactThreshMode;
voltMapping.laserRowArtifact.threshParam = laserArtifactThreshParam;
voltMapping.laserRowArtifact.mcMode = laserArtifactMcMode;
voltMapping.laserRowArtifact.mcUseGateColumnsOnly = mcUseGateColumnsOnly;
voltMapping.laserRowArtifact.mcSecondSweepForDff = laserArtifactMcSecondSweepForDff;

%% Calculate rough ROI mask on motion-corrected maxDvStack
input('Viewing motion-corrected mean image for selecting rough ROIs (continue or ctrl+c to stop!)');

figure(9); set(gcf, 'Position',  [100, 100, 1800, 900]); clf; 
colormap('winter'); imagesc(meanFluorMaxDvStack); axis equal; axis image; colorbar; set(gca, 'fontsize', 12);
nCells = input('How many neurons to analyze?: ', 's');
nCells = str2double(nCells);

% Hand select cell or area of interest, by freehand drawing (rough ROIs)
roughRoiXAllCells = cell(nCells, 1);
roughRoiYAllCells = cell(nCells, 1);
roiXAllCells_global = cell(nCells, 1);
roiYAllCells_global = cell(nCells, 1);
bkgrndRoiXAllCells = cell(nCells, 1);
bkgrndRoiYAllCells = cell(nCells, 1);

for nn = 1:nCells
    f1 = figure(10);
    set(gcf, 'Position',  [100, 100, 1800, 900]);
    clf
    colormap('winter'); imagesc(meanFluorMaxDvStack); axis equal; axis image; colorbar; set(gca, 'fontsize', 12);
    roughRoiX = []; roughRoiY = [];
    roiHandSelect = drawfreehand;
    roiHandSelectMask = createMask(roiHandSelect);
    [roughRoiX, roughRoiY] = find(roiHandSelectMask);

    roughRoiXAllCells{nn} = roughRoiX;
    roughRoiYAllCells{nn} = roughRoiY;
    close(f1);
end

% Derive one "global" fine ROI per cell from meanFluorMaxDvStack using
% the same gaussian + fibermetric approach, primarily for visualization
for nn = 1:nCells
    figure(10+(nn));
    set(gcf, 'Position',  [100, 100, 1600, 800]);
    clf
    subplot(4,1,1)
    colormap('winter'); imagesc(meanFluorMaxDvStack); axis equal; axis image; set(gca, 'fontsize', 12);

    roiMeanMaxDvStack = zeros(size(maxDvStack, 1), size(maxDvStack, 2));
    for rr = 1:length(roughRoiXAllCells{nn})
        roiMeanMaxDvStack(roughRoiXAllCells{nn}(rr), roughRoiYAllCells{nn}(rr)) = mean(maxDvStack(roughRoiXAllCells{nn}(rr), roughRoiYAllCells{nn}(rr),:), 3, 'omitnan');
    end
    
    subplot(4,1,2);
    imagesc(roiMeanMaxDvStack); axis equal; axis image; set(gca, 'fontsize', 12);

    % Same new ROI selection method as original script
    roiMaxStackDouble   = im2double(roiMeanMaxDvStack);
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
        thr  = prctile(valsR, 50);
    else
        thr = 0;
    end
    roiMaxStackRidgeReduced = roiMaxStackRidge;
    roiMaxStackRidgeReduced(roiMaxStackRidgeReduced < thr) = 0;
    roiMaxStackRidgeReduced(roiMaxStackRidgeReduced > 0) = 1;

    subplot(4,1,3)
    imagesc(roiMaxStackRidge); axis equal; axis image; set(gca, 'fontsize', 12);
    hold on
    [x, y] = find(roiMaxStackRidgeReduced);
    plot(y, x,'w.','MarkerSize',6);

    [roiX, roiY] = find(roiMaxStackRidgeReduced);
    if isempty(roiX)
        roiX = roughRoiXAllCells{nn};
        roiY = roughRoiYAllCells{nn};
    end
    roiXAllCells_global{nn} = roiX;
    roiYAllCells_global{nn} = roiY;

    subplot(4,1,4);
    imagesc(roiMaxStackRidge); axis equal; axis image; set(gca, 'fontsize', 12);
    hold on

    % Background ROI (unchanged logic from original script)
    innerBuffer = 2;
    ringWidth   = 3;
    minArea     = 50;
    
    innerSelect = imdilate(roiMaxStackRidge, strel('disk', innerBuffer));
    outerSelect = imdilate(roiMaxStackRidge, strel('disk', innerBuffer + ringWidth));
    backgroundRing  = outerSelect & ~innerSelect;
    
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

    [bkgrndRoiX, bkgrndRoiY] = find(ringClean);    
    bkgrndRoiXAllCells{nn} = bkgrndRoiX;
    bkgrndRoiYAllCells{nn} = bkgrndRoiY;

    [yBk,xBk] = find(ringClean);
    plot(xBk, yBk,'r.','MarkerSize',6);
    [xFine, yFine] = find(roiMaxStackRidgeReduced);
    plot(yFine, xFine,'w.','MarkerSize',6);
    hold off
end

% Summary ROI maps
allRois = zeros(size(meanFluorMaxDvStack, 1), size(meanFluorMaxDvStack, 2));
centerXY = [];
for nn = 1:nCells
    thisRoi = zeros(size(meanFluorMaxDvStack, 1), size(meanFluorMaxDvStack, 2)); 
    for rr = 1:length(roiXAllCells_global{nn})
        thisRoi(roiXAllCells_global{nn}(rr), roiYAllCells_global{nn}(rr)) = mean(meanFluorMaxDvStack(roiXAllCells_global{nn}(rr), roiYAllCells_global{nn}(rr),:), 3);
    end
    allRois = allRois + thisRoi;
    centerXY(nn, 1) = (min(roiXAllCells_global{nn}) + max(roiXAllCells_global{nn}))/2;
    centerXY(nn, 2) = (min(roiYAllCells_global{nn}) + max(roiYAllCells_global{nn}))/2;
end
figure(20); set(gcf, 'Position',  [100, 100, 1800, 800]); clf;
subplot(2, 1, 1);
colormap('winter'); imagesc(meanFluorMaxDvStack); axis equal; axis image; colorbar; set(gca, 'fontsize', 12);
subplot(2, 1, 2);
colormap('winter'); axis equal; axis image; colorbar; set(gca, 'fontsize', 12);
imagesc(allRois); axis equal; axis image; colorbar;
hold on;
for nn = 1:nCells
    plot(centerXY(nn, 2), centerXY(nn, 1), 'r+', 'LineWidth', 2, 'MarkerSize', 3);
    dx = 5; dy = 5;
    text(centerXY(nn, 2)+dx, centerXY(nn, 1)+dy, num2str(nn), ...
        'Color','w', 'FontSize', 18, 'FontWeight','bold', ...
        'HorizontalAlignment','left', 'VerticalAlignment','bottom');
end
hold off

allBkgrndRois = zeros(size(meanFluorMaxDvStack, 1), size(meanFluorMaxDvStack, 2));
centerXY = [];
for nn = 1:nCells
    thisRoi = zeros(size(meanFluorMaxDvStack, 1), size(meanFluorMaxDvStack, 2)); 
    for rr = 1:length(bkgrndRoiXAllCells{nn})
        thisRoi(bkgrndRoiXAllCells{nn}(rr), bkgrndRoiYAllCells{nn}(rr)) = mean(meanFluorMaxDvStack(bkgrndRoiXAllCells{nn}(rr), bkgrndRoiYAllCells{nn}(rr),:), 3);
    end
    allBkgrndRois = allBkgrndRois + thisRoi;
    centerXY(nn, 1) = (min(bkgrndRoiXAllCells{nn}) + max(bkgrndRoiXAllCells{nn}))/2;
    centerXY(nn, 2) = (min(bkgrndRoiYAllCells{nn}) + max(bkgrndRoiYAllCells{nn}))/2;
end
figure(21); set(gcf, 'Position',  [100, 100, 1800, 300]); clf; 
colormap('winter'); axis equal; axis image; colorbar; set(gca, 'fontsize', 12);
imagesc(allRois); axis equal; axis image; colorbar;
hold on;
[yBk2,xBk2] = find(allBkgrndRois);
plot(xBk2, yBk2,'r.','MarkerSize',6);
hold off;

voltMapping.nCells                 = nCells;
voltMapping.roughRoiXAllCells      = roughRoiXAllCells;
voltMapping.roughRoiYAllCells      = roughRoiYAllCells;
voltMapping.roiXAllCells_global    = roiXAllCells_global;
voltMapping.roiYAllCells_global    = roiYAllCells_global;
voltMapping.bkgrndRoiXAllCells     = bkgrndRoiXAllCells;
voltMapping.bkgrndRoiYAllCells     = bkgrndRoiYAllCells;
voltMapping.allRois                = allRois;

%% F, F0, dF, dF/F0 Calculation using motion-corrected stacks and per-trial fine ROIs
input('This step uses motion-corrected movies, computes trial-specific fine ROIs, and calculates dF/F (ctrl+c to stop!)');

% Preallocate containers for per-trial fine ROIs
fineRoiXAllCells = cell(nCells, 1);
fineRoiYAllCells = cell(nCells, 1);
bkgrndRoiXAllCells_trial = cell(nCells, 1);
bkgrndRoiYAllCells_trial = cell(nCells, 1);
for nn = 1:nCells
    fineRoiXAllCells{nn} = cell(length(imagesIndex), 1);
    fineRoiYAllCells{nn} = cell(length(imagesIndex), 1);
    bkgrndRoiXAllCells_trial{nn} = cell(length(imagesIndex), 1);
    bkgrndRoiYAllCells_trial{nn} = cell(length(imagesIndex), 1);
end

% VoltMapping "unpack to workspace" can define analysisStruct as a non-struct if that name
% exists as a voltMapping field; dF/F requires a struct here.
if exist('analysisStruct', 'var') && ~isstruct(analysisStruct)
    warning('VoltImg:analysisStructNotStruct', ...
        'analysisStruct was a %s, not a struct; clearing so dF/F can build a fresh struct.', class(analysisStruct));
    clear analysisStruct
end

% Preallocate data for each cell
holoSortedImagingCellNames = cell(nCells, 1);
filtHoloSortedImagingCellNames = cell(nCells, 1);
for nn = 1:nCells
    F0CellNames{nn}                    = ['F0AllTrials_', 'cell', num2str(nn)];
    roiMeanFCellNames{nn}              = ['roiMeanF_', 'cell', num2str(nn)];
    bkgrndMeanFCellNames{nn}           = ['bkgrndMeanF_', 'cell', num2str(nn)];
    subScalarCellNames{nn}             = ['subScalar_', 'cell', num2str(nn)];
    roiMeanFCorrectedCellNames{nn}     = ['roiMeanFCorrected_', 'cell', num2str(nn)];
    globalF0CellNames{nn}              = ['globalF0_', 'cell', num2str(nn)];
    dFCellNames{nn}                    = ['dF_', 'cell', num2str(nn)];
    % Note: field is literally e.g. dFFcell1 (no underscore before "cell"); use dFCellNames{nn} pattern for dF_/globalF0_.
    dFFCellNames{nn}                   = ['dFF', 'cell', num2str(nn)];
    holoSortedImagingCellNames{nn}     = ['holoSortedImagingAllTrials_', 'cell', num2str(nn)];
    filtHoloSortedImagingCellNames{nn} = ['filtHoloSortedImagingAllTrials_', 'cell', num2str(nn)];

    analysisStruct.(F0CellNames{nn})                    = cell(nConds, 1);
    analysisStruct.(holoSortedImagingCellNames{nn})     = cell(nConds, 1);
    analysisStruct.(filtHoloSortedImagingCellNames{nn}) = cell(nConds, 1);
    
    for cc = 1:nConds
        analysisStruct.(F0CellNames{nn}){cc}                    = cell(nHolos(cc), 1);
        analysisStruct.(holoSortedImagingCellNames{nn}){cc}     = cell(nHolos(cc), 1);
        analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc} = cell(nHolos(cc), 1);
    end
end

% Parallel-safe dF/F: pack each trial on workers, then merge in trial order (columns = trial index).
if exist('zeroDummySequence', 'var')
    zds = zeroDummySequence;
else
    zds = voltMapping.outParams.sequence{1, 2};
end

% Common early-trial F0 (all holograms per trial/cell); see VoltImg_mapping_parallel_dffOneTrial_commonEarlyF0.
% optional: set commonEarlyF0BaselineMs = startTime*1000 to match neuropil baseline length.
commonEarlyF0BaselineMs   = voltMapping.holoStimParams.startTime; % in ms
commonEarlyF0RollingWinMs = 50; % in ms

dffCtx = struct();
dffCtx.mcTiffFolder = mcTiffFolder;
dffCtx.imagesIndex = imagesIndex;
dffCtx.ImgfolderContents = ImgfolderContents;
dffCtx.nCells = nCells;
dffCtx.roughRoiXAllCells = roughRoiXAllCells;
dffCtx.roughRoiYAllCells = roughRoiYAllCells;
dffCtx.bkgrndRoiXAllCells = bkgrndRoiXAllCells;
dffCtx.bkgrndRoiYAllCells = bkgrndRoiYAllCells;
dffCtx.voltMapping = voltMapping;
dffCtx.excludeTrials = excludeTrials;
dffCtx.laserArtifactMcSecondSweepForDff = laserArtifactMcSecondSweepForDff;
dffCtx.laserArtifactGateColFirst = laserArtifactGateColFirst;
dffCtx.laserArtifactGateColLast = laserArtifactGateColLast;
dffCtx.laserArtifactThreshMode = laserArtifactThreshMode;
dffCtx.laserArtifactThreshParam = laserArtifactThreshParam;
dffCtx.startTime = startTime;
dffCtx.imagingFreq = imagingFreq;
dffCtx.preStimWindow = preStimWindow;
dffCtx.postStimWindow = postStimWindow;
dffCtx.ipi = ipi;
dffCtx.nPulses = nPulses;
dffCtx.UpOrDown = UpOrDown;
dffCtx.nHolos = nHolos;
dffCtx.zeroDummySequence = zds;
dffCtx.commonEarlyF0BaselineMs   = commonEarlyF0BaselineMs;
dffCtx.commonEarlyF0RollingWinMs = commonEarlyF0RollingWinMs;

trialDffPack = cell(nTrials, 1);
if useParallelDffTrials
    dffTrialProgressQ = parallel.pool.DataQueue;
    afterEach(dffTrialProgressQ, @(ttDone) fprintf(1, 'dF/F (per-trial / stim windows) completed trial %d of %d\n', ttDone, nTrials));
    parfor tt = 1:nTrials
        trialDffPack{tt} = VoltImg_mapping_parallel_dffOneTrial_commonEarlyF0(tt, dffCtx);
        send(dffTrialProgressQ, tt);
    end
else
    for tt = 1:nTrials
        disp(['Trial number: ', num2str(tt)]);
        trialDffPack{tt} = VoltImg_mapping_parallel_dffOneTrial_commonEarlyF0(tt, dffCtx);
    end
end

maxNumFrames = 0;
for tt = 1:nTrials
    p = trialDffPack{tt};
    if ~isstruct(p)
        error('VoltImg:trialDffPackBadElement', ...
            'trialDffPack{%d} must be a struct (class=%s). Re-run the per-trial dF/F loop; if you use parfor, check the pool and that trialDffPack is not overwritten.', ...
            tt, class(p));
    end
    if ~isfield(p, 'numFrames')
        error('VoltImg:trialDffPackNoNumFrames', 'trialDffPack{%d} has no numFrames field.', tt);
    end
    maxNumFrames = max(maxNumFrames, p.numFrames);
end
if maxNumFrames < 1
    error('VoltImg:maxNumFramesZero', ...
        'maxNumFrames is 0: every trial pack reported numFrames==0 (missing MC TIFFs under mcTiffFolder, or stacks failed to open).');
end
% Per-trial columns: ephys-excluded trials (p.isExcluded) never enter the else branch below, so globalF0/dF/dFF stay NaN for those tt.
for nn = 1:nCells
    analysisStruct.(roiMeanFCellNames{nn}) = nan(maxNumFrames, nTrials);
    analysisStruct.(bkgrndMeanFCellNames{nn}) = nan(maxNumFrames, nTrials);
    analysisStruct.(roiMeanFCorrectedCellNames{nn}) = nan(maxNumFrames, nTrials);
    analysisStruct.(subScalarCellNames{nn}) = nan(nTrials, 1);
    % Full-trial dF/dFF and one F0 per trial (same F0 for all holos in common-early-F0 pipeline).
    analysisStruct.(globalF0CellNames{nn}) = nan(nTrials, 1);
    analysisStruct.(dFCellNames{nn}) = nan(maxNumFrames, nTrials);
    analysisStruct.(dFFCellNames{nn}) = nan(maxNumFrames, nTrials);
end

for tt = 1:nTrials
    p = trialDffPack{tt};
    cc = p.cc;
    for nn = 1:nCells
        fineRoiXAllCells{nn}{tt} = p.fineRoiXAllCells{nn};
        fineRoiYAllCells{nn}{tt} = p.fineRoiYAllCells{nn};
        bkgrndRoiXAllCells_trial{nn}{tt} = p.bkgrndRoiXAllCells_trial{nn};
        bkgrndRoiYAllCells_trial{nn}{tt} = p.bkgrndRoiYAllCells_trial{nn};

        nf = p.numFrames;
        if p.isExcluded
            analysisStruct.(roiMeanFCellNames{nn})(1:nf, tt) = nan(nf, 1);
            analysisStruct.(bkgrndMeanFCellNames{nn})(1:nf, tt) = nan(nf, 1);
            analysisStruct.(roiMeanFCorrectedCellNames{nn})(1:nf, tt) = nan(nf, 1);
            analysisStruct.(subScalarCellNames{nn})(tt, 1) = nan;
        else
            analysisStruct.(roiMeanFCellNames{nn})(1:nf, tt) = p.roiMeanFcell{nn};
            analysisStruct.(bkgrndMeanFCellNames{nn})(1:nf, tt) = p.bkgrndMeanFcell{nn};
            analysisStruct.(roiMeanFCorrectedCellNames{nn})(1:nf, tt) = p.roiMeanFCorrectedCell{nn};
            analysisStruct.(subScalarCellNames{nn})(tt, 1) = p.alphaScalarVec(nn);
            f0Trial = p.f0Cell{nn}(1);
            analysisStruct.(globalF0CellNames{nn})(tt, 1) = f0Trial;
            roiC = p.roiMeanFCorrectedCell{nn};
            dFcol = roiC(:) - f0Trial;
            if isnan(f0Trial) || f0Trial == 0
                dFFcol = nan(nf, 1);
            else
                dFFcol = dFcol / f0Trial;
            end
            if UpOrDown == '2'
                dFFcol = -dFFcol;
            end
            analysisStruct.(dFCellNames{nn})(1:nf, tt) = dFcol;
            analysisStruct.(dFFCellNames{nn})(1:nf, tt) = dFFcol;
        end

        for hh = 1:p.nHolosThis
            slot = p.holoSeqThisTrial(hh);
            if p.isExcluded
                analysisStruct.(F0CellNames{nn}){cc}{slot, 1} = [analysisStruct.(F0CellNames{nn}){cc}{slot, 1}, NaN];
                analysisStruct.(holoSortedImagingCellNames{nn}){cc}{slot, 1} = [analysisStruct.(holoSortedImagingCellNames{nn}){cc}{slot, 1}, p.dFFCell{nn}{hh}];
                analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc}{slot, 1} = [analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc}{slot, 1}, p.filtDFFCell{nn}{hh}];
            else
                analysisStruct.(F0CellNames{nn}){cc}{slot, 1} = [analysisStruct.(F0CellNames{nn}){cc}{slot, 1}, p.f0Cell{nn}(hh)];
                analysisStruct.(holoSortedImagingCellNames{nn}){cc}{slot, 1} = [analysisStruct.(holoSortedImagingCellNames{nn}){cc}{slot, 1}, p.dFFCell{nn}{hh}];
                analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc}{slot, 1} = [analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc}{slot, 1}, p.filtDFFCell{nn}{hh}];
            end
        end
    end
end


%% Calculate mean response (and CI) for each hologram across trials and per condition
holoSortedMeanCellNames = cell(nCells, 1);
filtHoloSortedMeanCellNames = cell(nCells, 1);
for nn = 1:nCells
    holoSortedMeanCellNames{nn}        = ['holoSortedImagingMean_', 'cell', num2str(nn)];
    filtHoloSortedMeanCellNames{nn}    = ['filtHoloSortedImagingMean_', 'cell', num2str(nn)];
    analysisStruct.(holoSortedMeanCellNames{nn})        = cell(nConds, 1);
    analysisStruct.(filtHoloSortedMeanCellNames{nn})    = cell(nConds, 1);
    
    for cc = 1:nConds
        analysisStruct.(holoSortedMeanCellNames{nn}){cc}        = cell(nHolos(cc), 1);
        analysisStruct.(filtHoloSortedMeanCellNames{nn}){cc}    = cell(nHolos(cc), 1);
    end
end

for nn = 1:nCells
    for cc = 1:nConds
        for hh = 1:nHolos(cc)
            analysisStruct.(holoSortedMeanCellNames{nn}){cc}{hh}        = nanmean(analysisStruct.(holoSortedImagingCellNames{nn}){cc}{hh}, 2);
            analysisStruct.(filtHoloSortedMeanCellNames{nn}){cc}{hh}    = nanmean(analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc}{hh}, 2);        
        end
    end
end

CIDffAllCondsCellNames = cell(nCells, 1);
filtCIDffAllCondsCellNames = cell(nCells, 1);
for nn = 1:nCells
    CIDffAllCondsCellNames{nn}        = ['CIDffAllConds_', 'cell', num2str(nn)];
    filtCIDffAllCondsCellNames{nn}    = ['filtCIDffAllConds_', 'cell', num2str(nn)];
    analysisStruct.(CIDffAllCondsCellNames{nn})        = cell(nConds, 1);
    analysisStruct.(filtCIDffAllCondsCellNames{nn})    = cell(nConds, 1);
end

for nn = 1:nCells
    for cc = 1:nConds 
        for hh = 1:nHolos(cc)
            confidence_level = 0.95;
            means     = nanmean(analysisStruct.(holoSortedImagingCellNames{nn}){cc}{hh, 1}, 2);
            filtMeans = nanmean(analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc}{hh, 1}, 2);       
            std_errors     = std(analysisStruct.(holoSortedImagingCellNames{nn}){cc}{hh, 1}, 0, 2, "omitnan") / sqrt(size(analysisStruct.(holoSortedImagingCellNames{nn}){cc}{hh, 1}, 2));
            filtStd_errors = std(analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc}{hh, 1}, 0, 2, "omitnan") / sqrt(size(analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc}{hh, 1}, 2));
      
            t_score     = tinv((1 + confidence_level) / 2, size(analysisStruct.(holoSortedImagingCellNames{nn}){cc}{hh, 1}, 2) - 1);
            filtT_score = tinv((1 + confidence_level) / 2, size(analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc}{hh, 1}, 2) - 1);
            margin_of_error     = t_score * std_errors;
            filtMargin_of_error = filtT_score * filtStd_errors;
            lower_bounds     = means - margin_of_error;
            filtLower_bounds = filtMeans - filtMargin_of_error;
            upper_bounds     = means + margin_of_error;
            filtUpper_bounds = filtMeans + filtMargin_of_error;
            if UpOrDown == '2'
                analysisStruct.(CIDffAllCondsCellNames{nn}){cc}{hh, 1}        = [lower_bounds, upper_bounds];
                analysisStruct.(filtCIDffAllCondsCellNames{nn}){cc}{hh, 1}    = [filtLower_bounds, filtUpper_bounds];
            elseif UpOrDown =='1'
                analysisStruct.(CIDffAllCondsCellNames{nn}){cc}{hh, 1}        = [-lower_bounds, -upper_bounds];
                analysisStruct.(filtCIDffAllCondsCellNames{nn}){cc}{hh, 1}    = [-filtLower_bounds, -filtUpper_bounds];
            end
        end
    end
end

voltMapping.nPulseCoordsImaging = nPulseCoordsImaging;
voltMapping.holoSeqIndex        = holoSeqIndex;
voltMapping.preStimWindow       = preStimWindow;
voltMapping.postStimWindow      = postStimWindow;
voltMapping.fineRoiXAllCells    = fineRoiXAllCells;
voltMapping.fineRoiYAllCells    = fineRoiYAllCells;
voltMapping.bkgrndRoiXAllCells_trial = bkgrndRoiXAllCells_trial;
voltMapping.bkgrndRoiYAllCells_trial = bkgrndRoiYAllCells_trial;

%%
VoltImg_mapping_analysis_MultiCell_trialExcluder;

%% Reorganize voltMapping structs by cells
voltMapping.holoSortedImagingCellNames                  = holoSortedImagingCellNames;
voltMapping.filtHoloSortedImagingCellNames              = filtHoloSortedImagingCellNames;
voltMapping.holoSortedMeanCellNames                     = holoSortedMeanCellNames;
voltMapping.filtHoloSortedMeanCellNames                 = filtHoloSortedMeanCellNames;
voltMapping.CIDffAllCondsCellNames                      = CIDffAllCondsCellNames;
voltMapping.filtCIDffAllCondsCellNames                  = filtCIDffAllCondsCellNames;
voltMapping.stdImagingAllTrialsCellNames                = stdImagingAllTrialsCellNames;
voltMapping.stdFiltImagingAllTrialsCellNames            = stdFiltImagingAllTrialsCellNames;
voltMapping.exclHoloSortedImagingAllTrialsCellNames     = exclHoloSortedImagingAllTrialsCellNames;
voltMapping.exclFiltHoloSortedImagingAllTrialsCellNames = exclFiltHoloSortedImagingAllTrialsCellNames;
voltMapping.exclHoloSortedImagingMeanCellNames          = exclHoloSortedImagingMeanCellNames;
voltMapping.exclFiltHoloSortedImagingMeanCellNames      = exclFiltHoloSortedImagingMeanCellNames;

cellID = cell(nCells, 1);
for nn = 1:nCells
    cellID{nn} = ['Cell', num2str(nn)];
    structname = ['voltMapping.', 'Cell', num2str(nn)];
    eval([structname ' = struct();']);

    % Original Analysis data
    eval([structname '.' 'holoSortedImagingAllTrials' ' = analysisStruct.(holoSortedImagingCellNames{nn});']);
    eval([structname '.' 'filtHoloSortedImagingAllTrials' ' = analysisStruct.(filtHoloSortedImagingCellNames{nn});']);
    eval([structname '.' 'holoSortedImagingMean' ' = analysisStruct.(holoSortedMeanCellNames{nn});']);
    eval([structname '.' 'filtHoloSortedImagingMean' ' = analysisStruct.(filtHoloSortedMeanCellNames{nn});']);
    eval([structname '.' 'CIDffAllConds' ' = analysisStruct.(CIDffAllCondsCellNames{nn});']);
    eval([structname '.' 'filtCIDffAllConds' ' = analysisStruct.(filtCIDffAllCondsCellNames{nn});']);
    % Full-trial imaging traces (frames x trials); holoSorted* above are per-stim windows only.
    eval([structname '.' 'roiMeanF_allTrials' ' = analysisStruct.(roiMeanFCellNames{nn});']);
    eval([structname '.' 'roiMeanFCorrected_allTrials' ' = analysisStruct.(roiMeanFCorrectedCellNames{nn});']);
    eval([structname '.' 'globalF0_perTrial' ' = analysisStruct.(globalF0CellNames{nn});']);
    eval([structname '.' 'dF_fullTrial' ' = analysisStruct.(dFCellNames{nn});']);
    eval([structname '.' 'dFF_fullTrial' ' = analysisStruct.(dFFCellNames{nn});']);
    eval([structname '.' 'F0AllTrials_byHolo' ' = analysisStruct.(F0CellNames{nn});']);
    
    % Exclusion criteria applied data
    eval([structname '.' 'stdImagingAllTrialsCellNames' ' = analysisStruct.(stdImagingAllTrialsCellNames{nn});']);
    eval([structname '.' 'stdFiltImagingAllTrialsCellNames' ' = analysisStruct.(stdFiltImagingAllTrialsCellNames{nn});']);    
    eval([structname '.' 'exclHoloSortedImagingAllTrials' ' = analysisStruct.(exclHoloSortedImagingAllTrialsCellNames{nn});']);    
    eval([structname '.' 'exclFiltHoloSortedImagingAllTrials' ' = analysisStruct.(exclFiltHoloSortedImagingAllTrialsCellNames{nn});']);
    eval([structname '.' 'exclHoloSortedImagingMean' ' = analysisStruct.(exclHoloSortedImagingMeanCellNames{nn});']);
    eval([structname '.' 'exclFiltHoloSortedImagingMean' ' = analysisStruct.(exclFiltHoloSortedImagingMeanCellNames{nn});']);
    eval([structname '.' 'exclCIDffAllConds' ' = analysisStruct.(exclCIDffAllCondsCellNames{nn});']);
    eval([structname '.' 'exclFiltCIDffAllConds' ' = analysisStruct.(exclFiltCIDffAllCondsCellNames{nn});']);
    
    % Background subtracted Analysis data
%     eval([structname '.' 'holoSortedImagingAllTrials_bkgrndsubtrct' ' = analysisStruct.(holoSortedImagingCellNames_bkgrndsubtrct{nn});']);
%     eval([structname '.' 'filtHoloSortedImagingAllTrials_bkgrndsubtrct' ' = analysisStruct.(filtHoloSortedImagingCellNames_bkgrndsubtrct{nn});']);
%     eval([structname '.' 'holoSortedImagingMean_bkgrndsubtrct' ' = analysisStruct.(holoSortedMeanCellNames_bkgrndsubtrct{nn});']);
%     eval([structname '.' 'filtHoloSortedImagingMean_bkgrndsubtrct' ' = analysisStruct.(filtHoloSortedMeanCellNames_bkgrndsubtrct{nn});']);
%     eval([structname '.' 'CIDffAllConds_bkgrndsubtrct' ' = analysisStruct.(CIDffAllCondsCellNames_bkgrndsubtrct{nn});']);
%     eval([structname '.' 'filtCIDffAllConds_bkgrndsubtrct' ' = analysisStruct.(filtCIDffAllCondsCellNames_bkgrndsubtrct{nn});']);
    
end
voltMapping.cellID = cellID;

%% Align Exclusion-applied imaging with ephys traces
nn = double(input('which cell number? '));

exclCIDffAllConds = voltMapping.(cellID{nn}).exclCIDffAllConds;
exclHoloSortedImagingMean = voltMapping.(cellID{nn}).exclHoloSortedImagingMean;
exclFiltCIDffAllConds = voltMapping.(cellID{nn}).exclFiltCIDffAllConds;
exclFiltHoloSortedImagingMean = voltMapping.(cellID{nn}).exclFiltHoloSortedImagingMean;

for cc = 1:nConds
    for hh = 1:nHolos(cc)
        figure(cc*1000+hh);
        set(gcf, 'Position',  [100, 100, 500, 375])
        clf
        % subplot(1,2,1)
        % hold on
        % 
        % fill([linspace(0, size(exclCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclCIDffAllConds{cc}{hh}, 1)), fliplr(linspace(0, size(exclCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclCIDffAllConds{cc}{hh}, 1)))],...
        %     [exclCIDffAllConds{cc}{hh}(:, 1)'*100, fliplr(exclCIDffAllConds{cc}{hh}(:, 2)'*100)], [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
        % % plot CI lowerbound
        % plot(linspace(0, size(exclCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclCIDffAllConds{cc}{hh}, 1)), exclCIDffAllConds{cc}{hh}(:, 1)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);
        % % plot CI upperbound
        % plot(linspace(0, size(exclCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclCIDffAllConds{cc}{hh}, 1)), exclCIDffAllConds{cc}{hh}(:, 2)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);
        % 
        % %         % plot ephys and voltage traces
        % %         [ax, hl, h2] = plotyy(linspace(0, size(holoSortedImagingMean{cc}(:, hh), 1)/imagingFreq, size(holoSortedImagingMean{cc}(:, hh), 1)), holoSortedImagingMean{cc}(:, hh), linspace(0, size(holoSortedDataMean{cc}(:, hh), 1)/Fs, size(holoSortedDataMean{cc}(:, hh), 1)), holoSortedDataMean{cc}(:, hh));
        % % %         set(ax, 'XAxisLocation', 'origin');
        % % %         ylim([min(holoSortedImagingMean{cc}(:, hh)), max(holoSortedImagingMean20ƒ{cc}(:, hh))]);
        % % %         ylim([-1, 2])
        % %         ax(2).YLim = ax(1).YLim;
        % %         set(ax(1), 'ytick', floor(min(holoSortedImagingMean{cc}(:, hh))):1:ceil(max(holoSortedImagingMean{cc}(:, hh))));
        % %         set(ax, {'ycolor'},{'g';'k'});
        % %         set(ax, 'FontSize', 14);
        % %         ylabel(ax(1), 'dF/F');
        % %         ylabel(ax(2), 'dV');
        % %         xlabel(ax(2), 'Time (s)');
        % %         set(hl, 'Color', 'g');
        % %         set(h2, 'Color', [0.7 0.7 0.7]);
        % %         set(hl, 'LineWidth', 2);
        % 
        % % plot voltage imaging trace only
        % % yyaxis right
        % plot(linspace(0, size(exclHoloSortedImagingMean{cc}{hh}, 1)/imagingFreq, size(exclHoloSortedImagingMean{cc}{hh}, 1)), exclHoloSortedImagingMean{cc}{hh}*100, '-', 'linewidth', 2, 'color', 'g');
        % ylabel('dF/F (%)');
        % xlabel('Time (s)');
        % %         ylim([-1 max(holoSortedImagingMean{cc}{hh})*100])
        % %         ylim([min(holoSortedImagingMean{cc}{hh})*100 max(holoSortedImagingMean{cc}{hh})*100])
        % ylim([-0.5 1.5])
        % ax.YColor = [0 1 0];
        % % xticks([0:0.05:size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);
        % 
        % % plot ephys trace only
        % %         yyaxis left
        % %     %     axes('Position',[.70 .12 .2 .2]);
        % %     %     box on
        % %         plot(linspace(0, size(holoSortedDataMean{cc}(:, hh), 1)/Fs, size(holoSortedDataMean{cc}(:, hh), 1)), holoSortedDataMean{cc}(:, hh), 'linewidth', 1.5, 'color', [0 0 0]);
        % % gca;
        % % set(gca,'xtick',[], 'fontsize', 18);
        % % %         ylim([min(holoSortedDataMean{cc}(:, hh)), max(holoSortedDataMean{cc}(:, hh))]);
        % %         ylim([-0.5 2]);
        % %         ylabel('dV');
        % %         ax.YColor = [0, 0, 0];
        % % axis off
        % %         xlim([0 size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);
        % %         xticks([0:0.05:size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);
        % 
        % for nn = 1:length(nPulseCoords)
        %     xline(nPulseCoords(nn)/Fs, '-', 'LineWidth', 10, 'color', [1 0 0]);
        % end
        % % plot([0.01; 0.035], [0.5; 0.5], '-k', 'LineWidth', 2);
        % % plot([0.01; 0.01], [0.5; 1], '-k', 'LineWidth', 2);
        % hold off
        % subplot(1,2,1)
        % hold on
        % 
        % fill([linspace(0, size(exclCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclCIDffAllConds{cc}{hh}, 1)), fliplr(linspace(0, size(exclCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclCIDffAllConds{cc}{hh}, 1)))],...
        %     [exclCIDffAllConds{cc}{hh}(:, 1)'*100, fliplr(exclCIDffAllConds{cc}{hh}(:, 2)'*100)], [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
        % % plot CI lowerbound
        % plot(linspace(0, size(exclCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclCIDffAllConds{cc}{hh}, 1)), exclCIDffAllConds{cc}{hh}(:, 1)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);
        % % plot CI upperbound
        % plot(linspace(0, size(exclCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclCIDffAllConds{cc}{hh}, 1)), exclCIDffAllConds{cc}{hh}(:, 2)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);
        % 
        % %         % plot ephys and voltage traces
        % %         [ax, hl, h2] = plotyy(linspace(0, size(holoSortedImagingMean{cc}(:, hh), 1)/imagingFreq, size(holoSortedImagingMean{cc}(:, hh), 1)), holoSortedImagingMean{cc}(:, hh), linspace(0, size(holoSortedDataMean{cc}(:, hh), 1)/Fs, size(holoSortedDataMean{cc}(:, hh), 1)), holoSortedDataMean{cc}(:, hh));
        % % %         set(ax, 'XAxisLocation', 'origin');
        % % %         ylim([min(holoSortedImagingMean{cc}(:, hh)), max(holoSortedImagingMean20ƒ{cc}(:, hh))]);
        % % %         ylim([-1, 2])
        % %         ax(2).YLim = ax(1).YLim;
        % %         set(ax(1), 'ytick', floor(min(holoSortedImagingMean{cc}(:, hh))):1:ceil(max(holoSortedImagingMean{cc}(:, hh))));
        % %         set(ax, {'ycolor'},{'g';'k'});
        % %         set(ax, 'FontSize', 14);
        % %         ylabel(ax(1), 'dF/F');
        % %         ylabel(ax(2), 'dV');
        % %         xlabel(ax(2), 'Time (s)');
        % %         set(hl, 'Color', 'g');
        % %         set(h2, 'Color', [0.7 0.7 0.7]);
        % %         set(hl, 'LineWidth', 2);
        % 
        % % plot voltage imaging trace only
        % % yyaxis right
        % plot(linspace(0, size(exclHoloSortedImagingMean{cc}{hh}, 1)/imagingFreq, size(exclHoloSortedImagingMean{cc}{hh}, 1)), exclHoloSortedImagingMean{cc}{hh}*100, '-', 'linewidth', 2, 'color', 'g');
        % ylabel('dF/F (%)');
        % xlabel('Time (s)');
        % %         ylim([-1 max(holoSortedImagingMean{cc}{hh})*100])
        % %         ylim([min(holoSortedImagingMean{cc}{hh})*100 max(holoSortedImagingMean{cc}{hh})*100])
        % ylim([-0.5 1.5])
        % ax.YColor = [0 1 0];
        % % xticks([0:0.05:size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);
        % 
        % % plot ephys trace only
        % %         yyaxis left
        % %     %     axes('Position',[.70 .12 .2 .2]);
        % %     %     box on
        % %         plot(linspace(0, size(holoSortedDataMean{cc}(:, hh), 1)/Fs, size(holoSortedDataMean{cc}(:, hh), 1)), holoSortedDataMean{cc}(:, hh), 'linewidth', 1.5, 'color', [0 0 0]);
        % % gca;
        % % set(gca,'xtick',[], 'fontsize', 18);
        % % %         ylim([min(holoSortedDataMean{cc}(:, hh)), max(holoSortedDataMean{cc}(:, hh))]);
        % %         ylim([-0.5 2]);
        % %         ylabel('dV');
        % %         ax.YColor = [0, 0, 0];
        % % axis off
        % %         xlim([0 size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);
        % %         xticks([0:0.05:size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);
        % 
        % for nn = 1:length(nPulseCoords)
        %     xline(nPulseCoords(nn)/Fs, '-', 'LineWidth', 10, 'color', [1 0 0]);
        % end
        % % plot([0.01; 0.035], [0.5; 0.5], '-k', 'LineWidth', 2);
        % % plot([0.01; 0.01], [0.5; 1], '-k', 'LineWidth', 2);
        % hold off
        % subplot(1,2,1)
        % hold on
        % 
        % fill([linspace(0, size(exclCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclCIDffAllConds{cc}{hh}, 1)), fliplr(linspace(0, size(exclCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclCIDffAllConds{cc}{hh}, 1)))],...
        %     [exclCIDffAllConds{cc}{hh}(:, 1)'*100, fliplr(exclCIDffAllConds{cc}{hh}(:, 2)'*100)], [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
        % % plot CI lowerbound
        % plot(linspace(0, size(exclCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclCIDffAllConds{cc}{hh}, 1)), exclCIDffAllConds{cc}{hh}(:, 1)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);
        % % plot CI upperbound
        % plot(linspace(0, size(exclCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclCIDffAllConds{cc}{hh}, 1)), exclCIDffAllConds{cc}{hh}(:, 2)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);
        % 
        % %         % plot ephys and voltage traces
        % %         [ax, hl, h2] = plotyy(linspace(0, size(holoSortedImagingMean{cc}(:, hh), 1)/imagingFreq, size(holoSortedImagingMean{cc}(:, hh), 1)), holoSortedImagingMean{cc}(:, hh), linspace(0, size(holoSortedDataMean{cc}(:, hh), 1)/Fs, size(holoSortedDataMean{cc}(:, hh), 1)), holoSortedDataMean{cc}(:, hh));
        % % %         set(ax, 'XAxisLocation', 'origin');
        % % %         ylim([min(holoSortedImagingMean{cc}(:, hh)), max(holoSortedImagingMean20ƒ{cc}(:, hh))]);
        % % %         ylim([-1, 2])
        % %         ax(2).YLim = ax(1).YLim;
        % %         set(ax(1), 'ytick', floor(min(holoSortedImagingMean{cc}(:, hh))):1:ceil(max(holoSortedImagingMean{cc}(:, hh))));
        % %         set(ax, {'ycolor'},{'g';'k'});
        % %         set(ax, 'FontSize', 14);
        % %         ylabel(ax(1), 'dF/F');
        % %         ylabel(ax(2), 'dV');
        % %         xlabel(ax(2), 'Time (s)');
        % %         set(hl, 'Color', 'g');
        % %         set(h2, 'Color', [0.7 0.7 0.7]);
        % %         set(hl, 'LineWidth', 2);
        % 
        % % plot voltage imaging trace only
        % % yyaxis right
        % plot(linspace(0, size(exclHoloSortedImagingMean{cc}{hh}, 1)/imagingFreq, size(exclHoloSortedImagingMean{cc}{hh}, 1)), exclHoloSortedImagingMean{cc}{hh}*100, '-', 'linewidth', 2, 'color', 'g');
        % ylabel('dF/F (%)');
        % xlabel('Time (s)');
        % %         ylim([-1 max(holoSortedImagingMean{cc}{hh})*100])
        % %         ylim([min(holoSortedImagingMean{cc}{hh})*100 max(holoSortedImagingMean{cc}{hh})*100])
        % ylim([-0.5 1.5])
        % ax.YColor = [0 1 0];
        % % xticks([0:0.05:size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);
        % 
        % % plot ephys trace only
        % %         yyaxis left
        % %     %     axes('Position',[.70 .12 .2 .2]);
        % %     %     box on
        % %         plot(linspace(0, size(holoSortedDataMean{cc}(:, hh), 1)/Fs, size(holoSortedDataMean{cc}(:, hh), 1)), holoSortedDataMean{cc}(:, hh), 'linewidth', 1.5, 'color', [0 0 0]);
        % % gca;
        % % set(gca,'xtick',[], 'fontsize', 18);
        % % %         ylim([min(holoSortedDataMean{cc}(:, hh)), max(holoSortedDataMean{cc}(:, hh))]);
        % %         ylim([-0.5 2]);
        % %         ylabel('dV');
        % %         ax.YColor = [0, 0, 0];
        % % axis off
        % %         xlim([0 size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);
        % %         xticks([0:0.05:size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);
        % 
        % for nn = 1:length(nPulseCoords)
        %     xline(nPulseCoords(nn)/Fs, '-', 'LineWidth', 10, 'color', [1 0 0]);
        % end
        % % plot([0.01; 0.035], [0.5; 0.5], '-k', 'LineWidth', 2);
        % % plot([0.01; 0.01], [0.5; 1], '-k', 'LineWidth', 2);
        % hold off
        % subplot(1,2,1)
        % hold on
        % 
        % fill([linspace(0, size(exclCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclCIDffAllConds{cc}{hh}, 1)), fliplr(linspace(0, size(exclCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclCIDffAllConds{cc}{hh}, 1)))],...
        %     [exclCIDffAllConds{cc}{hh}(:, 1)'*100, fliplr(exclCIDffAllConds{cc}{hh}(:, 2)'*100)], [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
        % % plot CI lowerbound
        % plot(linspace(0, size(exclCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclCIDffAllConds{cc}{hh}, 1)), exclCIDffAllConds{cc}{hh}(:, 1)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);
        % % plot CI upperbound
        % plot(linspace(0, size(exclCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclCIDffAllConds{cc}{hh}, 1)), exclCIDffAllConds{cc}{hh}(:, 2)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);
        % 
        % %         % plot ephys and voltage traces
        % %         [ax, hl, h2] = plotyy(linspace(0, size(holoSortedImagingMean{cc}(:, hh), 1)/imagingFreq, size(holoSortedImagingMean{cc}(:, hh), 1)), holoSortedImagingMean{cc}(:, hh), linspace(0, size(holoSortedDataMean{cc}(:, hh), 1)/Fs, size(holoSortedDataMean{cc}(:, hh), 1)), holoSortedDataMean{cc}(:, hh));
        % % %         set(ax, 'XAxisLocation', 'origin');
        % % %         ylim([min(holoSortedImagingMean{cc}(:, hh)), max(holoSortedImagingMean20ƒ{cc}(:, hh))]);
        % % %         ylim([-1, 2])
        % %         ax(2).YLim = ax(1).YLim;
        % %         set(ax(1), 'ytick', floor(min(holoSortedImagingMean{cc}(:, hh))):1:ceil(max(holoSortedImagingMean{cc}(:, hh))));
        % %         set(ax, {'ycolor'},{'g';'k'});
        % %         set(ax, 'FontSize', 14);
        % %         ylabel(ax(1), 'dF/F');
        % %         ylabel(ax(2), 'dV');
        % %         xlabel(ax(2), 'Time (s)');
        % %         set(hl, 'Color', 'g');
        % %         set(h2, 'Color', [0.7 0.7 0.7]);
        % %         set(hl, 'LineWidth', 2);
        % 
        % % plot voltage imaging trace only
        % % yyaxis right
        % plot(linspace(0, size(exclHoloSortedImagingMean{cc}{hh}, 1)/imagingFreq, size(exclHoloSortedImagingMean{cc}{hh}, 1)), exclHoloSortedImagingMean{cc}{hh}*100, '-', 'linewidth', 2, 'color', 'g');
        % ylabel('dF/F (%)');
        % xlabel('Time (s)');
        % %         ylim([-1 max(holoSortedImagingMean{cc}{hh})*100])
        % %         ylim([min(holoSortedImagingMean{cc}{hh})*100 max(holoSortedImagingMean{cc}{hh})*100])
        % ylim([-0.5 1.5])
        % ax.YColor = [0 1 0];
        % % xticks([0:0.05:size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);
        % 
        % % plot ephys trace only
        % %         yyaxis left
        % %     %     axes('Position',[.70 .12 .2 .2]);
        % %     %     box on
        % %         plot(linspace(0, size(holoSortedDataMean{cc}(:, hh), 1)/Fs, size(holoSortedDataMean{cc}(:, hh), 1)), holoSortedDataMean{cc}(:, hh), 'linewidth', 1.5, 'color', [0 0 0]);
        % % gca;
        % % set(gca,'xtick',[], 'fontsize', 18);
        % % %         ylim([min(holoSortedDataMean{cc}(:, hh)), max(holoSortedDataMean{cc}(:, hh))]);
        % %         ylim([-0.5 2]);
        % %         ylabel('dV');
        % %         ax.YColor = [0, 0, 0];
        % % axis off
        % %         xlim([0 size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);
        % %         xticks([0:0.05:size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);
        % 
        % for nn = 1:length(nPulseCoords)
        %     xline(nPulseCoords(nn)/Fs, '-', 'LineWidth', 10, 'color', [1 0 0]);
        % end
        % % plot([0.01; 0.035], [0.5; 0.5], '-k', 'LineWidth', 2);
        % % plot([0.01; 0.01], [0.5; 1], '-k', 'LineWidth', 2);
        % hold off
        % subplot(1,2,1)
        % hold on
        % 
        % fill([linspace(0, size(exclCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclCIDffAllConds{cc}{hh}, 1)), fliplr(linspace(0, size(exclCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclCIDffAllConds{cc}{hh}, 1)))],...
        %     [exclCIDffAllConds{cc}{hh}(:, 1)'*100, fliplr(exclCIDffAllConds{cc}{hh}(:, 2)'*100)], [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
        % % plot CI lowerbound
        % plot(linspace(0, size(exclCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclCIDffAllConds{cc}{hh}, 1)), exclCIDffAllConds{cc}{hh}(:, 1)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);
        % % plot CI upperbound
        % plot(linspace(0, size(exclCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclCIDffAllConds{cc}{hh}, 1)), exclCIDffAllConds{cc}{hh}(:, 2)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);
        % 
        % %         % plot ephys and voltage traces
        % %         [ax, hl, h2] = plotyy(linspace(0, size(holoSortedImagingMean{cc}(:, hh), 1)/imagingFreq, size(holoSortedImagingMean{cc}(:, hh), 1)), holoSortedImagingMean{cc}(:, hh), linspace(0, size(holoSortedDataMean{cc}(:, hh), 1)/Fs, size(holoSortedDataMean{cc}(:, hh), 1)), holoSortedDataMean{cc}(:, hh));
        % % %         set(ax, 'XAxisLocation', 'origin');
        % % %         ylim([min(holoSortedImagingMean{cc}(:, hh)), max(holoSortedImagingMean20ƒ{cc}(:, hh))]);
        % % %         ylim([-1, 2])
        % %         ax(2).YLim = ax(1).YLim;
        % %         set(ax(1), 'ytick', floor(min(holoSortedImagingMean{cc}(:, hh))):1:ceil(max(holoSortedImagingMean{cc}(:, hh))));
        % %         set(ax, {'ycolor'},{'g';'k'});
        % %         set(ax, 'FontSize', 14);
        % %         ylabel(ax(1), 'dF/F');
        % %         ylabel(ax(2), 'dV');
        % %         xlabel(ax(2), 'Time (s)');
        % %         set(hl, 'Color', 'g');
        % %         set(h2, 'Color', [0.7 0.7 0.7]);
        % %         set(hl, 'LineWidth', 2);
        % 
        % % plot voltage imaging trace only
        % % yyaxis right
        % plot(linspace(0, size(exclHoloSortedImagingMean{cc}{hh}, 1)/imagingFreq, size(exclHoloSortedImagingMean{cc}{hh}, 1)), exclHoloSortedImagingMean{cc}{hh}*100, '-', 'linewidth', 2, 'color', 'g');
        % ylabel('dF/F (%)');
        % xlabel('Time (s)');
        % %         ylim([-1 max(holoSortedImagingMean{cc}{hh})*100])
        % %         ylim([min(holoSortedImagingMean{cc}{hh})*100 max(holoSortedImagingMean{cc}{hh})*100])
        % ylim([-0.5 1.5])
        % ax.YColor = [0 1 0];
        % % xticks([0:0.05:size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);
        % 
        % % plot ephys trace only
        % %         yyaxis left
        % %     %     axes('Position',[.70 .12 .2 .2]);
        % %     %     box on
        % %         plot(linspace(0, size(holoSortedDataMean{cc}(:, hh), 1)/Fs, size(holoSortedDataMean{cc}(:, hh), 1)), holoSortedDataMean{cc}(:, hh), 'linewidth', 1.5, 'color', [0 0 0]);
        % % gca;
        % % set(gca,'xtick',[], 'fontsize', 18);
        % % %         ylim([min(holoSortedDataMean{cc}(:, hh)), max(holoSortedDataMean{cc}(:, hh))]);
        % %         ylim([-0.5 2]);
        % %         ylabel('dV');
        % %         ax.YColor = [0, 0, 0];
        % % axis off
        % %         xlim([0 size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);
        % %         xticks([0:0.05:size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);
        % 
        % for nn = 1:length(nPulseCoords)
        %     xline(nPulseCoords(nn)/Fs, '-', 'LineWidth', 10, 'color', [1 0 0]);
        % end
        % % plot([0.01; 0.035], [0.5; 0.5], '-k', 'LineWidth', 2);
        % % plot([0.01; 0.01], [0.5; 1], '-k', 'LineWidth', 2);
        % hold off
        % subplot(1,2,1)
        % hold on
        % 
        % fill([linspace(0, size(exclCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclCIDffAllConds{cc}{hh}, 1)), fliplr(linspace(0, size(exclCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclCIDffAllConds{cc}{hh}, 1)))],...
        %     [exclCIDffAllConds{cc}{hh}(:, 1)'*100, fliplr(exclCIDffAllConds{cc}{hh}(:, 2)'*100)], [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
        % % plot CI lowerbound
        % plot(linspace(0, size(exclCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclCIDffAllConds{cc}{hh}, 1)), exclCIDffAllConds{cc}{hh}(:, 1)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);
        % % plot CI upperbound
        % plot(linspace(0, size(exclCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclCIDffAllConds{cc}{hh}, 1)), exclCIDffAllConds{cc}{hh}(:, 2)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);
        % 
        % %         % plot ephys and voltage traces
        % %         [ax, hl, h2] = plotyy(linspace(0, size(holoSortedImagingMean{cc}(:, hh), 1)/imagingFreq, size(holoSortedImagingMean{cc}(:, hh), 1)), holoSortedImagingMean{cc}(:, hh), linspace(0, size(holoSortedDataMean{cc}(:, hh), 1)/Fs, size(holoSortedDataMean{cc}(:, hh), 1)), holoSortedDataMean{cc}(:, hh));
        % % %         set(ax, 'XAxisLocation', 'origin');
        % % %         ylim([min(holoSortedImagingMean{cc}(:, hh)), max(holoSortedImagingMean20ƒ{cc}(:, hh))]);
        % % %         ylim([-1, 2])
        % %         ax(2).YLim = ax(1).YLim;
        % %         set(ax(1), 'ytick', floor(min(holoSortedImagingMean{cc}(:, hh))):1:ceil(max(holoSortedImagingMean{cc}(:, hh))));
        % %         set(ax, {'ycolor'},{'g';'k'});
        % %         set(ax, 'FontSize', 14);
        % %         ylabel(ax(1), 'dF/F');
        % %         ylabel(ax(2), 'dV');
        % %         xlabel(ax(2), 'Time (s)');
        % %         set(hl, 'Color', 'g');
        % %         set(h2, 'Color', [0.7 0.7 0.7]);
        % %         set(hl, 'LineWidth', 2);
        % 
        % % plot voltage imaging trace only
        % % yyaxis right
        % plot(linspace(0, size(exclHoloSortedImagingMean{cc}{hh}, 1)/imagingFreq, size(exclHoloSortedImagingMean{cc}{hh}, 1)), exclHoloSortedImagingMean{cc}{hh}*100, '-', 'linewidth', 2, 'color', 'g');
        % ylabel('dF/F (%)');
        % xlabel('Time (s)');
        % %         ylim([-1 max(holoSortedImagingMean{cc}{hh})*100])
        % %         ylim([min(holoSortedImagingMean{cc}{hh})*100 max(holoSortedImagingMean{cc}{hh})*100])
        % ylim([-0.5 1.5])
        % ax.YColor = [0 1 0];
        % % xticks([0:0.05:size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);
        % 
        % % plot ephys trace only
        % %         yyaxis left
        % %     %     axes('Position',[.70 .12 .2 .2]);
        % %     %     box on
        % %         plot(linspace(0, size(holoSortedDataMean{cc}(:, hh), 1)/Fs, size(holoSortedDataMean{cc}(:, hh), 1)), holoSortedDataMean{cc}(:, hh), 'linewidth', 1.5, 'color', [0 0 0]);
        % % gca;
        % % set(gca,'xtick',[], 'fontsize', 18);
        % % %         ylim([min(holoSortedDataMean{cc}(:, hh)), max(holoSortedDataMean{cc}(:, hh))]);
        % %         ylim([-0.5 2]);
        % %         ylabel('dV');
        % %         ax.YColor = [0, 0, 0];
        % % axis off
        % %         xlim([0 size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);
        % %         xticks([0:0.05:size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);
        % 
        % for nn = 1:length(nPulseCoords)
        %     xline(nPulseCoords(nn)/Fs, '-', 'LineWidth', 10, 'color', [1 0 0]);
        % end
        % % plot([0.01; 0.035], [0.5; 0.5], '-k', 'LineWidth', 2);
        % % plot([0.01; 0.01], [0.5; 1], '-k', 'LineWidth', 2);
        % hold off

        % subplot(1,2,2)
        hold on

        fill([linspace(0, size(exclFiltCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclFiltCIDffAllConds{cc}{hh}, 1)), fliplr(linspace(0, size(exclFiltCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclFiltCIDffAllConds{cc}{hh}, 1)))],...
            [exclFiltCIDffAllConds{cc}{hh}(:, 1)'*100, fliplr(exclFiltCIDffAllConds{cc}{hh}(:, 2)'*100)], [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
        % plot CI lowerbound
        plot(linspace(0, size(exclFiltCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclFiltCIDffAllConds{cc}{hh}, 1)), exclFiltCIDffAllConds{cc}{hh}(:, 1)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);
        % plot CI upperbound
        plot(linspace(0, size(exclFiltCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclFiltCIDffAllConds{cc}{hh}, 1)), exclFiltCIDffAllConds{cc}{hh}(:, 2)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);

        %         % plot ephys and voltage traces
        %         [ax, hl, h2] = plotyy(linspace(0, size(holoSortedImagingMean{cc}(:, hh), 1)/imagingFreq, size(holoSortedImagingMean{cc}(:, hh), 1)), holoSortedImagingMean{cc}(:, hh), linspace(0, size(holoSortedDataMean{cc}(:, hh), 1)/Fs, size(holoSortedDataMean{cc}(:, hh), 1)), holoSortedDataMean{cc}(:, hh));
        % %         set(ax, 'XAxisLocation', 'origin');
        % %         ylim([min(holoSortedImagingMean{cc}(:, hh)), max(holoSortedImagingMean{cc}(:, hh))]);
        % %         ylim([-1, 2])
        %         ax(2).YLim = ax(1).YLim;
        %         set(ax(1), 'ytick', floor(min(holoSortedImagingMean{cc}(:, hh))):1:ceil(max(holoSortedImagingMean{cc}(:, hh))));
        %         set(ax, {'ycolor'},{'g';'k'});
        %         set(ax, 'FontSize', 14);
        %         ylabel(ax(1), 'dF/F');
        %         ylabel(ax(2), 'dV');
        %         xlabel(ax(2), 'Time (s)');
        %         set(hl, 'Color', 'g');
        %         set(h2, 'Color', [0.7 0.7 0.7]);
        %         set(hl, 'LineWidth', 2);

        % plot voltage imaging trace only
        % yyaxis right
        ylabel('dF/F (%)');
        xlabel('Time (s)');
        %         ylim([-1 max(holoSortedImagingMean{cc}{hh})*100])
        %         ylim([min(holoSortedImagingMean{cc}{hh})*100 max(holoSortedImagingMean{cc}{hh})*100])
        ylim([-1 2])
        
        if cc > 2
            ylim([-1 4])
        end

        ax.YColor = [0 1 0];
        % xticks([0:0.05:size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);

        % plot ephys trace only
        %         yyaxis left
        %     %     axes('Position',[.70 .12 .2 .2]);
        %     %     box on
        %         plot(linspace(0, size(holoSortedDataMean{cc}(:, hh), 1)/Fs, size(holoSortedDataMean{cc}(:, hh), 1)), holoSortedDataMean{cc}(:, hh), 'linewidth', 1.5, 'color', [0 0 0]);
        % gca;
        % set(gca,'xtick',[], 'fontsize', 18);
        % %         ylim([min(holoSortedDataMean{cc}(:, hh)), max(holoSortedDataMean{cc}(:, hh))]);
        %         ylim([-0.5 2]);
        %         ylabel('dV');
        %         ax.YColor = [0, 0, 0];
        % axis off
        %         xlim([0 size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);
        %         xticks([0:0.05:size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);

        pdMs = voltMapping.outParams.pulseDur(:);
        pulseDurMsHere = pdMs(min(cc, numel(pdMs)));
        if pulseDurMsHere <= 0 && ~isempty(pulseDurs), pulseDurMsHere = pulseDurs(1); end
        pulseDurSecHere = pulseDurMsHere / 1000;
        ylStim = ylim;
        for nn = 1:length(nPulseCoords)
            tOn = nPulseCoords(nn) / Fs;
            patch([tOn, tOn + pulseDurSecHere, tOn + pulseDurSecHere, tOn], ...
                [ylStim(1), ylStim(1), ylStim(2), ylStim(2)], [1, 0, 0], ...
                'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HitTest', 'off');
        end
        plot(linspace(0, size(exclFiltHoloSortedImagingMean{cc}{hh}, 1)/imagingFreq, size(exclFiltHoloSortedImagingMean{cc}{hh}, 1)), exclFiltHoloSortedImagingMean{cc}{hh}*100, '-', 'linewidth', 2, 'color', 'g');
        % plot([0.01; 0.035], [0.5; 0.5], '-k', 'LineWidth', 2);
        % plot([0.01; 0.01], [0.5; 1], '-k', 'LineWidth', 2);
        hold off
        pause
    end
end

%% Four-panel comparison: mean + 95% CI for all pipelines (1×4 subplots)
cellIdx = double(input('which cell number (4-panel mean + CI)? '));

holoSortedImagingMean         = voltMapping.(cellID{cellIdx}).holoSortedImagingMean;
filtHoloSortedImagingMean     = voltMapping.(cellID{cellIdx}).filtHoloSortedImagingMean;
exclHoloSortedImagingMean     = voltMapping.(cellID{cellIdx}).exclHoloSortedImagingMean;
exclFiltHoloSortedImagingMean = voltMapping.(cellID{cellIdx}).exclFiltHoloSortedImagingMean;
CIDffAllConds                 = voltMapping.(cellID{cellIdx}).CIDffAllConds;
filtCIDffAllConds             = voltMapping.(cellID{cellIdx}).filtCIDffAllConds;
exclCIDffAllConds             = voltMapping.(cellID{cellIdx}).exclCIDffAllConds;
exclFiltCIDffAllConds         = voltMapping.(cellID{cellIdx}).exclFiltCIDffAllConds;

meanSeries = {holoSortedImagingMean, filtHoloSortedImagingMean, exclHoloSortedImagingMean, exclFiltHoloSortedImagingMean};
ciSeries   = {CIDffAllConds, filtCIDffAllConds, exclCIDffAllConds, exclFiltCIDffAllConds};
yScales    = [10, 10, 10, 10];
subTitles  = {'Mean', 'Filt mean', 'Excl mean', 'Excl filt mean'};

FsEphys = voltMapping.daqParams.Fs;
plotEphysWithImaging = exist('ePhysAvail', 'var') && ePhysAvail == 1 ...
    && isfield(voltMapping, 'ephys') && isfield(voltMapping.ephys, 'holoSortedDataMean');

for cc = 1:nConds
    for hh = 1:nHolos(cc)
        figure(20000 + cc*1000 + hh);
        % Extra horizontal room so yyaxis-right ylabels are not clipped
        set(gcf, 'Position', [100, 100, 2400, 420]);
        clf

        nT = numel(holoSortedImagingMean{cc}{hh});
        tAxis = linspace(0, nT/imagingFreq, nT);
        pdMs = voltMapping.outParams.pulseDur(:);
        pulseDurMsHere = pdMs(min(cc, numel(pdMs)));
        if pulseDurMsHere <= 0 && ~isempty(pulseDurs), pulseDurMsHere = pulseDurs(1); end
        pulseDurSecHere = pulseDurMsHere / 1000;

        ephysOnImgGrid = [];
        if plotEphysWithImaging && hh <= size(voltMapping.ephys.holoSortedDataMean{cc}, 2)
            ephysMeanCol = voltMapping.ephys.holoSortedDataMean{cc}(:, hh);
            nEp = numel(ephysMeanCol);
            if nEp > 1
                tEp = linspace(0, nEp/FsEphys, nEp);
                ephysOnImgGrid = interp1(tEp(:), double(ephysMeanCol(:)), tAxis(:), 'linear', 'extrap');
            elseif nEp == 1
                ephysOnImgGrid = repmat(double(ephysMeanCol), numel(tAxis), 1);
            end
        end

        ephysRgb = [0.9, 0.9, 0.9];
        for sp = 1:4
            ax = subplot(1, 4, sp);
            m = meanSeries{sp}{cc}{hh}*10;
            ci = ciSeries{sp}{cc}{hh}*10;
            ys = yScales(sp);

            yyaxis(ax, 'left')
            hold(ax, 'on')
            fill(ax, [tAxis, fliplr(tAxis)], [ci(:, 1)'*ys, fliplr(ci(:, 2)'*ys)], ...
                [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
            ylStim = ax.YLim;
            for pulseIdx = 1:length(nPulseCoords)
                tOn = nPulseCoords(pulseIdx) / FsEphys;
                patch(ax, [tOn, tOn + pulseDurSecHere, tOn + pulseDurSecHere, tOn], ...
                    [ylStim(1), ylStim(1), ylStim(2), ylStim(2)], [1, 0, 0], ...
                    'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HitTest', 'off');
            end
            plot(ax, tAxis, ci(:, 1)*ys, '--', 'linewidth', 1, 'color', [0.7, 0.7, 0.7]);
            plot(ax, tAxis, ci(:, 2)*ys, '--', 'linewidth', 1, 'color', [0.7, 0.7, 0.7]);
            imgMeanLine = m*ys;
            plot(ax, tAxis, imgMeanLine, '-', 'linewidth', 2, 'color', 'g');
            ylabel(ax, 'dF/F (%)');

            if plotEphysWithImaging && ~isempty(ephysOnImgGrid)
                imgPk = max(abs(imgMeanLine(:)), [], 'omitnan');
                epPk = max(abs(ephysOnImgGrid(:)), [], 'omitnan');
                yyaxis(ax, 'right')
                hold(ax, 'on')
                scaleEp = nan;
                if isfinite(imgPk) && isfinite(epPk) && epPk > 0 && imgPk > 0
                    scaleEp = imgPk / epPk;
                    ephysScaledRow = ephysOnImgGrid(:).' * scaleEp;
                    plot(ax, tAxis, ephysScaledRow, '-', 'linewidth', 2, 'color', ephysRgb);
                else
                    plot(ax, tAxis, ephysOnImgGrid(:).', '-', 'linewidth', 2, 'color', ephysRgb);
                end
                if isfinite(scaleEp) && scaleEp > 0
                    yRightLbl = sprintf('Rel. Vm (ephys); value / %.3g = mV', scaleEp);
                else
                    yRightLbl = 'Vm ephys (mV, rel. to baseline)';
                end
                ylabel(ax, yRightLbl, 'Color', ephysRgb, 'Interpreter', 'none');
                ax.YAxis(2).Color = ephysRgb;
                ax.YAxis(2).Label.Visible = 'on';
                ax.YAxis(2).Label.FontWeight = 'normal';
            end

            yyaxis(ax, 'left')
            hold(ax, 'off')
            xlabel(ax, 'Time (s)');
            title(ax, subTitles{sp});
        end
        if plotEphysWithImaging && ~isempty(ephysOnImgGrid)
            sgtitle(sprintf('Cell %d — cond %d, holo %d (blue ephys: right axis; peak-matched to green)', cellIdx, cc, hh));
        else
            sgtitle(sprintf('Cell %d — cond %d, holo %d (4-panel)', cellIdx, cc, hh));
        end
        pause
    end
end

%% Publication (light mode): filt mean + 95% CI + ephys, shared y-scale per condition, scale bars
% One figure per (condition, hologram); advance with pause like the four-panel section. For each
% condition, all holos share identical xlim, left ylim (dF/F percent), and right ylim (rel. Vm, mV).
% Pre-stim mean subtracted from filt + ephys so baselines meet at y=0; symmetric ylims align zeros.
% Scale bars: 1%% dF/F, 0.5 mV Vm (both vertical bars on the left), 10 ms horizontal.
if ~exist('cellIdx', 'var')
    cellIdxPub = double(input('which cell number (publication filt + scale bars)? '));
else
    cellIdxPub = cellIdx;
end
filtPubMean     = voltMapping.(cellID{cellIdxPub}).filtHoloSortedImagingMean;
filtPubCI       = voltMapping.(cellID{cellIdxPub}).filtCIDffAllConds;
FsEphysPub      = voltMapping.daqParams.Fs;
plotEphysPub    = exist('ePhysAvail', 'var') && ePhysAvail == 1 ...
    && isfield(voltMapping, 'ephys') && isfield(voltMapping.ephys, 'holoSortedDataMean');
scaleBarPct     = 1;    % dF/F scale bar height (percent)
scaleBar_mV     = 1;    % Vm scale bar height (mV)
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

for ccPub = 1:nConds
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
        yBarBotL = yLimL(1) + 0.34 * ylRng;
        plot(axh, [xBarPct, xBarPct], yBarBotL + [0, scaleBarPct], '-', 'Color', 'k', ...
            'LineWidth', 1.2, 'Clipping', 'off', 'HandleVisibility', 'off');
        text(axh, xBarPct - 0.030 * xlR, yBarBotL + 0.5 * scaleBarPct, sprintf('1%s dF/F', '%'), ...
            'FontName', pubFont, 'FontSize', pubFontSize, 'Color', 'k', ...
            'Rotation', 90, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', ...
            'Interpreter', 'none', 'Clipping', 'off');

        tBarW = min(timeBarSec, 0.85 * xlR);
        yBarTime = yLimL(1) + 0.29 * ylRng;
        plot(axh, xBarPct + [0, tBarW], [yBarTime, yBarTime], '-', 'Color', 'k', ...
            'LineWidth', 1.2, 'Clipping', 'off', 'HandleVisibility', 'off');
        text(axh, xBarPct + 0.5 * tBarW, yBarTime - 0.018 * ylRng, sprintf('%.0f ms', tBarW * 1000), ...
            'FontName', pubFont, 'FontSize', pubFontSize - 0.5, 'Color', 'k', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'Interpreter', 'none', 'Clipping', 'off');

        if ~isempty(ephysLine) && ~isempty(yLimR)
            ylRngV = yLimR(2) - yLimR(1);
            yBarBotV = yLimR(1) + 0.34 * ylRngV;
            yyaxis(axh, 'right')
            ylim(axh, yLimR)
            xlim(axh, xLimAll)
            plot(axh, [xBarVmL, xBarVmL], yBarBotV + [0, scaleBar_mV], '-', 'Color', 'k', ...
                'LineWidth', 1.2, 'Clipping', 'off', 'HandleVisibility', 'off');
            ymidVm = yBarBotV + 0.5 * scaleBar_mV;
            vFrac = (ymidVm - yLimR(1)) / max(yLimR(2) - yLimR(1), eps);
            yLblVm = yLimL(1) + vFrac * ylRng;
            yyaxis(axh, 'left')
            text(axh, xBarVmL - 0.030 * xlR, yLblVm, sprintf('%.1f mV', scaleBar_mV), ...
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

%% Power curve: peak mean dF/F (filt vs excl-filt) vs stim power
if ~exist('cellIdx', 'var')
    cellIdx = double(input('which cell number (power curve)? '));
end
filtHoloSortedImagingMean     = voltMapping.(cellID{cellIdx}).filtHoloSortedImagingMean;
exclFiltHoloSortedImagingMean = voltMapping.(cellID{cellIdx}).exclFiltHoloSortedImagingMean;
powerVec = voltMapping.outParams.power(:);
maxHoloAcrossConds = max(nHolos);
% Plateau: mean peak dF/F at the highest nPlateauPowerPts tested powers (asymptote estimate).
nPlateauPowerPts = 3;
holoColors = lines(max(maxHoloAcrossConds, 2));

figure(31001); clf;
set(gcf, 'Position', [120, 120, 600, 450]);
hold on
for hh = 1:maxHoloAcrossConds
    pwr = []; pk = [];
    for cc = 1:nConds
        if hh <= nHolos(cc)
            tr = filtHoloSortedImagingMean{cc}{hh};
            pwr(end+1) = powerVec(cc);
            pk(end+1)  = max(tr(:), [], 'omitnan');
        end
    end
    if ~isempty(pwr)
        [pwrSorted, ord] = sort(pwr);
        pkSorted = pk(ord);
        kPl = min(nPlateauPowerPts, numel(pkSorted));
        plateauPct = mean(pkSorted(end - kPl + 1:end), 'omitnan') * 100;
        c = holoColors(hh, :);
        plot(pwrSorted, pkSorted * 100, '-o', 'LineWidth', 1.5, 'Color', c, 'MarkerFaceColor', c, ...
            'DisplayName', sprintf('holo %d (plateau ~ %.2f%%)', hh, plateauPct));
        yline(plateauPct, '--', 'Color', c, 'LineWidth', 1, 'Alpha', 0.55, 'HandleVisibility', 'off');
    end
end
hold off
xlabel('Power');
ylabel('Peak mean dF/F (%)');
title(sprintf('Cell %d — filtHoloSortedImagingMean (plateau = mean at top %d powers)', cellIdx, nPlateauPowerPts));
legend('Location', 'best');
grid on

figure(31002); clf;
set(gcf, 'Position', [140, 90, 600, 450]);
hold on
for hh = 1:maxHoloAcrossConds
    pwr = []; pk = [];
    for cc = 1:nConds
        if hh <= nHolos(cc)
            tr = exclFiltHoloSortedImagingMean{cc}{hh};
            pwr(end+1) = powerVec(cc);
            pk(end+1)  = max(tr(:), [], 'omitnan');
        end
    end
    if ~isempty(pwr)
        [pwrSorted, ord] = sort(pwr);
        pkSorted = pk(ord);
        kPl = min(nPlateauPowerPts, numel(pkSorted));
        plateauPct = mean(pkSorted(end - kPl + 1:end), 'omitnan') * 100;
        c = holoColors(hh, :);
        plot(pwrSorted, pkSorted * 100, '-o', 'LineWidth', 1.5, 'Color', c, 'MarkerFaceColor', c, ...
            'DisplayName', sprintf('holo %d (plateau ~ %.2f%%)', hh, plateauPct));
        yline(plateauPct, '--', 'Color', c, 'LineWidth', 1, 'Alpha', 0.55, 'HandleVisibility', 'off');
    end
end
hold off
xlabel('Power');
ylabel('Peak mean dF/F (%)');
title(sprintf('Cell %d — exclFiltHoloSortedImagingMean (plateau = mean at top %d powers)', cellIdx, nPlateauPowerPts));
legend('Location', 'best');
grid on

%% Optional sanity-check visualization: per-trial ROI vs neuropil ring overlap
doNeuropilSanityPlot = true;   % set false to skip; samples trials every sanityTrialStep
sanityTrialStep = 300;          % trials 1, 1+step, 1+2*step, ... up to nTrials
if doNeuropilSanityPlot
    sanityTrialIndices = 1:sanityTrialStep:nTrials;
    nSanityPanels = numel(sanityTrialIndices);
    if nSanityPanels < 1
        disp('Neuropil sanity plot: no trials in range.');
    else
        imgH = size(meanFluorMaxDvStack, 1);
        imgW = size(meanFluorMaxDvStack, 2);
        % Single column: larger, easier to inspect than a tight grid
        ncols = 1;
        nrows = nSanityPanels;
        panelW = 700;
        panelH = 480;
        cellColors = lines(max(nCells, 2));

        figure(99); clf;
        set(gcf, 'Position', [80, 40, panelW, panelH * nSanityPanels + 140]);

        for ip = 1:nSanityPanels
            sanityTrialIdx = sanityTrialIndices(ip);

            allTrialRoiMaskDbg = false(imgH, imgW);
            for cc = 1:nCells
                xTmp = fineRoiXAllCells{cc}{sanityTrialIdx};
                yTmp = fineRoiYAllCells{cc}{sanityTrialIdx};
                if ~isempty(xTmp)
                    allTrialRoiMaskDbg(sub2ind([imgH, imgW], xTmp, yTmp)) = true;
                end
            end

            overlapPerCell = zeros(nCells, 1);
            for cc = 1:nCells
                neuropilMaskThis = false(imgH, imgW);
                if ~isempty(bkgrndRoiXAllCells_trial{cc}{sanityTrialIdx})
                    neuropilMaskThis(sub2ind([imgH, imgW], ...
                        bkgrndRoiXAllCells_trial{cc}{sanityTrialIdx}, ...
                        bkgrndRoiYAllCells_trial{cc}{sanityTrialIdx})) = true;
                end
                overlapPerCell(cc) = nnz(neuropilMaskThis & allTrialRoiMaskDbg);
            end

            disp(['Sanity check trial ', num2str(sanityTrialIdx), ': overlap(np, ANY ROI) per cell = ', ...
                mat2str(overlapPerCell'), ' | max = ', num2str(max(overlapPerCell))]);

            subplot(nrows, ncols, ip);
            imagesc(meanFluorMaxDvStack); axis image; colormap(gray); hold on;
            [rAll, cAll] = find(allTrialRoiMaskDbg);
            plot(cAll, rAll, '.', 'Color', [1, 1, 0], 'MarkerSize', 3);
            for cc = 1:nCells
                xN = bkgrndRoiXAllCells_trial{cc}{sanityTrialIdx};
                yN = bkgrndRoiYAllCells_trial{cc}{sanityTrialIdx};
                if ~isempty(xN)
                    plot(yN, xN, '.', 'Color', cellColors(cc, :), 'MarkerSize', 5);
                end
            end
            title({['trial ', num2str(sanityTrialIdx)], ...
                ['max overlap(any ROI)=', num2str(max(overlapPerCell))]}, 'FontSize', 9);
            if ip == 1
                legend({'All ROIs (trial)', 'Neuropil: color = cell #'}, ...
                    'TextColor', 'k', 'Location', 'best', 'FontSize', 8);
            end
            hold off;
        end
        sgtitle({['Neuropil sanity: all ', num2str(nCells), ' cells; neuropil color = cell index'], ...
            ['every ', num2str(sanityTrialStep), ' trials (1:', num2str(sanityTrialStep), ':', num2str(nTrials), ')']});
    end
end

%% Save Analysis Results
expID = num2str(mouseID);
voltMapping.mouseID = ['voltMapping_Analysis_', expID, '_MultiCellAnalysis_MCfineROI_laserRowArtifact_parallel'];
% Persist analysis flags inside voltMapping for robust reloads.
if exist('UpOrDown', 'var')
    voltMapping.UpOrDown = UpOrDown;
end
if exist('ePhysAvail', 'var')
    voltMapping.ePhysAvail = ePhysAvail;
end
directory = '/Users/masatosadahiro/Documents/Data/Voltage Imaging/Voltage Imaging/voltMapping/Analysis Results/Analysis_newMCanddFFcalc';
fileName = [num2str(voltMapping.mouseID), '.mat'];

% Persist everything needed to reload traces, re-index analysisStruct, re-run MC/dF/F
% logic, and reproduce figures (paths, filters, ROI name keys, trial excluder outputs).
if ~exist('confidence_level', 'var') || isempty(confidence_level)
    confidence_level = 0.95;
end
analysisSaveTimestamp = datetime('now');
varsToSave = { ...
    'voltMapping', 'analysisStruct', ...
    'F0CellNames', 'roiMeanFCellNames', 'bkgrndMeanFCellNames', 'subScalarCellNames', ...
    'roiMeanFCorrectedCellNames', 'globalF0CellNames', 'dFCellNames', 'dFFCellNames', ...
    'holoSortedImagingCellNames', 'filtHoloSortedImagingCellNames', ...
    'holoSortedMeanCellNames', 'filtHoloSortedMeanCellNames', ...
    'CIDffAllCondsCellNames', 'filtCIDffAllCondsCellNames', ...
    'stdImagingAllTrialsCellNames', 'stdFiltImagingAllTrialsCellNames', ...
    'exclHoloSortedImagingAllTrialsCellNames', 'exclFiltHoloSortedImagingAllTrialsCellNames', ...
    'exclHoloSortedImagingMeanCellNames', 'exclFiltHoloSortedImagingMeanCellNames', ...
    'exclCIDffAllCondsCellNames', 'exclFiltCIDffAllCondsCellNames', ...
    'mouseID', 'ePhysAvail', ...
    'nCells', 'nTrials', 'nConds', 'nHolos', 'imagingFreq', 'Fs', 'ipi', 'nPulses', ...
    'preStimWindow', 'postStimWindow', 'startTime', 'UpOrDown', 'excludeTrials', ...
    'powers', 'nextHoloDelay', 'pulseDurs', 'trialTime', ...
    'ImgfolderContents', 'savePath', 'saveDirectory', 'directory', 'fileName', ...
    'normcorrePath', ...
    'cutOffFreq', 'blp', 'alp', 'cutOffFreqIm', 'bIm', 'aIm', ...
    'vThreshold', 'confidence_level', ...
    'mappingInputs', 'allBkgrndRois', 'nPlateauPowerPts', ...
    'useLaserRowArtifactFilter', 'laserArtifactGateColFirst', 'laserArtifactGateColLast', ...
    'laserArtifactThreshMode', 'laserArtifactThreshParam', 'laserArtifactMcMode', ...
    'mcUseGateColumnsOnly', 'laserArtifactMcSecondSweepForDff', ...
    'cellID', 'analysisSaveTimestamp'};

if exist('zeroDummySequence', 'var')
    varsToSave{end+1} = 'zeroDummySequence';
end

varsToSave = varsToSave(cellfun(@(n) exist(n, 'var') == 1, varsToSave));
save(fullfile(directory, fileName), varsToSave{:}, '-v7.3');

TimeNow = datetime('now', 'Format', 'dd-MMM-yyyy HH:mm:ss');
disp(['finished saving at: ' char(TimeNow)])

%% Load Analysis Results
% Restores the same workspace names as "Save Analysis Results" (varsToSave)
% plus fields of voltMapping / voltMapping.ephys for downstream sections.
%
% Usage:
%   A) load(fullfile(directory,fileName)) yourself, then run this section, or
%   B) run this section with no voltMapping in workspace: you will be asked
%      to pick the .mat (or it loads directory/fileName if those exist).

destWs = 'caller'; % invoking workspace (command window / script base)

if ~exist('voltMapping', 'var') || isempty(voltMapping)
    matPath = '';
    if exist('directory', 'var') && exist('fileName', 'var')
        cand = fullfile(directory, fileName);
        if exist(cand, 'file') == 2
            matPath = cand;
        end
    end
    if isempty(matPath)
        [fn, pth] = uigetfile('*.mat', 'Select saved voltMapping analysis .mat');
        if isequal(fn, 0)
            return;
        end
        matPath = fullfile(pth, fn);
    end
    S = load(matPath);
    topNames = fieldnames(S);
    for ii = 1:numel(topNames)
        assignin(destWs, topNames{ii}, S.(topNames{ii}));
    end
    clear S topNames ii matPath
    if exist('fn', 'var')
        clear fn pth
    end
end

if ~exist('voltMapping', 'var') || isempty(voltMapping)
    error('LoadAnalysisResults:NoVoltMapping', 'voltMapping missing after load.');
end

% Unpack voltMapping fields (daqParams, outParams, per-cell Cell#, etc.)
vmNames = fieldnames(voltMapping);
for i = 1:numel(vmNames)
    assignin(destWs, vmNames{i}, voltMapping.(vmNames{i}));
end

% Alias for comments / snippets that still say ExpStruct
assignin(destWs, 'ExpStruct', voltMapping);
if isfield(voltMapping, 'ExpStruct2')
    assignin(destWs, 'ExpStruct2', voltMapping.ExpStruct2);
end

% Unpack nested ephys (also appears as top-level fields after this loop)
if isfield(voltMapping, 'ephys') && isstruct(voltMapping.ephys)
    epNames = fieldnames(voltMapping.ephys);
    for i = 1:numel(epNames)
        assignin(destWs, epNames{i}, voltMapping.ephys.(epNames{i}));
    end
end

% Canonical scalars / vectors (aligned with header setup and voltMapping.* copies)
Fs = voltMapping.daqParams.Fs;
if isfield(voltMapping, 'imagingFreq') && ~isempty(voltMapping.imagingFreq)
    imagingFreq = voltMapping.imagingFreq;
elseif isfield(voltMapping, 'sampleFreq')
    imagingFreq = voltMapping.sampleFreq;
else
    error('LoadAnalysisResults:NoImagingFreq', 'voltMapping has neither imagingFreq nor sampleFreq.');
end
trialTime = voltMapping.daqParams.maxSweepLengthSec;
nTrials = length(voltMapping.trialCond);
outParams = voltMapping.outParams;
powers = outParams.power;
if isfield(voltMapping, 'nConds') && ~isempty(voltMapping.nConds)
    nConds = voltMapping.nConds;
else
    nConds = length(outParams.sequence);
end
if isfield(voltMapping, 'nHolos') && ~isempty(voltMapping.nHolos)
    nHolos = voltMapping.nHolos;
else
    nHolos = voltMapping.holoStimParams.nHolos;
    nHolos(1) = max(nHolos);
end
if isfield(voltMapping, 'pulseDurs') && ~isempty(voltMapping.pulseDurs)
    pulseDurs = voltMapping.pulseDurs;
else
    pulseDurs = unique(outParams.pulseDur);
    pulseDurs = nonzeros(pulseDurs);
end
if isfield(voltMapping, 'nPulses') && ~isempty(voltMapping.nPulses)
    nPulses = voltMapping.nPulses;
else
    nPulses = unique(outParams.nPulses);
    nPulses = nonzeros(nPulses);
end
if isfield(voltMapping, 'ipi') && ~isempty(voltMapping.ipi)
    ipi = voltMapping.ipi;
else
    ipi = unique(outParams.ipi);
    ipi = nonzeros(ipi);
end
if isfield(voltMapping, 'nextHoloDelay') && ~isempty(voltMapping.nextHoloDelay)
    nextHoloDelay = voltMapping.nextHoloDelay;
else
    nextHoloDelay = unique(voltMapping.holoStimParams.nextHoloDelay);
    nextHoloDelay = nonzeros(nextHoloDelay);
end
startTime = (voltMapping.holoStimParams.startTime) / 1000;
if isfield(voltMapping, 'imagesIndex')
    imagesIndex = voltMapping.imagesIndex;
end
if isfield(voltMapping, 'UpOrDown')
    UpOrDown = voltMapping.UpOrDown;
end
if isfield(voltMapping, 'ephysFilePath')
    ephysFilePath = voltMapping.ephysFilePath;
end
if isfield(voltMapping, 'ImgsFilePath')
    ImgsFilePath = voltMapping.ImgsFilePath;
end
if isfield(voltMapping, 'cellID')
    cellID = voltMapping.cellID;
end
if isfield(voltMapping, 'ePhysAvail')
    ePhysAvail = voltMapping.ePhysAvail;
end

disp('Load Analysis Results: voltMapping, analysisStruct, name lists, filters, and paths restored; structs unpacked to workspace.');
