%% Voltage Imaging Mapping Analysis 021226 newDFFcalc and MCfineROI (global neuropil)
% Same as VoltImg_mapping_analysis_MultiCell_newDFF_021226_MCfineROI.m except:
% Neuropil/background subtraction uses the global ring from the mean motion-corrected
% stack (fibermetric ring per cell), not trial-specific rings and not ROI-exclusion logic.
% Per-trial fine ROIs for soma F are unchanged.
%
% Pipeline:
% 1) Per-trial NoRMCorre motion correction; 2) maxDvStack (50%% random eligible trials, first 4s frames per plane);
% 3) rough ROIs;
% 4) global fine ROI + global neuropil ring per cell on mean image;
% 5) per-trial fine ROIs; F with global neuropil subtraction for dF/F.

%%
clear all
close all

%% Load files and setup
% Step 1: Read the ephys file
% ephysFilePath = char(uigetdir('E:\Voltage Imaging\VoltMapping\Ephys Data')); % Select and set root folder where all experiments with cells you want to analyze are located
% ephysFilePath = char(uigetdir('/Volumes/ExData2/Voltage Imaging/VoltMapping/Ephys Data')); % Select and set root folder where all experiments with cells you want to analyze are located
ephysFilePath = char(uigetdir('/Volumes/ExData3/voltMapping/ephys data'));

ephysFileDir = dir(ephysFilePath);
load([ephysFileDir(end).folder, '/', ephysFileDir(end).name]);
disp(ephysFileDir(end).name);
 
% Step 2: Identify the folder containing the imaging files correspdonding to the ephys file
% ImgsFilePath = char(uigetdir('E:\Voltage Imaging\VoltMapping\Imaging Data')); % Select and set root folder where all experiments with cells you want to analyze are located
% ImgsFilePath = char(uigetdir('/Volumes/ExData2/Voltage Imaging/VoltMapping/Imaging Data')); % Select and set root folder where all experiments with cells you want to analyze are located
ImgsFilePath = char(uigetdir('/Volumes/ExData3/voltMapping/imaging data'));

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

% Step 4: Identify GEVI type
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

% Step 5: Save directory (later overwritten near save step)
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
    
    sweepThisTrial = filtfilt(blp, alp, mappingInputsBaselined{tt});
    
    for hh = 1:nHolos(voltMapping.trialCond(tt, 1))
        if isempty(voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)})
            voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)} = voltMapping.outParams.firstStimTimes{1, 2};
        end

        thisHoloSweep = sweepThisTrial((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*Fs:((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*Fs+((ipi*nPulses+preStimWindow+postStimWindow)/1000)*Fs));
        thisHoloSweep = thisHoloSweep - mean(thisHoloSweep(1:nPulseCoords(1)-100));
   
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

%% Motion-correct all imaging trials with NoRMCorre and build maxDvStack
input('Running NoRMCorre on all imaging trials, saving motion-corrected stacks, and building maxDvStack (continue or ctrl+c to stop!)');

% Setup save directory for motion-corrected images
expID = num2str(mouseID);
voltMapping.mouseID = ['voltMapping_Analysis_', expID, '_MultiCellAnalysis_MCfineROI'];
savePath = '/Volumes/phoenixinthesky/Masato/Voltage Imaging Data_Phoenix/voltMapping/MC Imaging Data';
saveDirectory = fullfile(savePath, num2str(voltMapping.mouseID));
if ~exist(saveDirectory, 'dir')
    mkdir(saveDirectory);
end

% Create subfolder for motion-corrected TIFFs
mcTiffFolder = fullfile(saveDirectory, 'Motion_Corrected_Tiffs');
if ~exist(mcTiffFolder, 'dir')
    mkdir(mcTiffFolder);
end

% maxDvStack sampling: 50%% of eligible (non-excluded) trials at random; first 4s of frames per mean plane
nImgTrials = length(imagesIndex);
eligibleTrialTT = [];
for ttx = 1:nImgTrials
    if ~ismember(imagesIndex(ttx), excludeTrials)
        eligibleTrialTT(end+1) = ttx; %#ok<AGROW>
    end
end
[maxDvTrialMask, maxDvFrameCap] = VoltImg_mapping_maxDvStackSamplingPlan(nImgTrials, imagingFreq, eligibleTrialTT);

% Preallocate maxDvStack using first TIFF dimensions
firstImgPath = fullfile(ImgfolderContents(imagesIndex(1)).folder, ImgfolderContents(imagesIndex(1)).name);
infoFirst = imfinfo(firstImgPath);
maxDvStack = nan(infoFirst(1).Height, infoFirst(1).Width, length(imagesIndex));

% Pass 1: build one global template from trial mean images so all trials
% are registered into the same reference frame.
globalTemplateAccum = zeros(infoFirst(1).Height, infoFirst(1).Width, 'single');
nTemplateTrials = 0;
for tt = 1:length(imagesIndex)
    disp(['Global template registering trial ', num2str(tt), ' of ', num2str(length(imagesIndex))]);
    if ismember(imagesIndex(tt), excludeTrials)
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

    globalTemplateAccum = globalTemplateAccum + mean(single(imageStack), 3);
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

for tt = 1:length(imagesIndex)
    disp(['Motion correcting trial ', num2str(tt), ' of ', num2str(length(imagesIndex))]);

    % Optionally skip excluded ephys trials
    if ismember(imagesIndex(tt), excludeTrials)
        continue
    end

    % Load raw multi-page TIFF (same "new loading method" as original)
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
    [H,W] = size(firstFrame);

    imageStack = zeros(H, W, nKeep, 'like', firstFrame);
    imageStack(:,:,1) = firstFrame;

    for ki = 2:nKeep
        t.setDirectory(dirList(ki));
        imageStack(:,:,ki) = t.read();
    end
    t.close();

    % NoRMCorre rigid motion correction on this trial (anchored to one
    % global template so all stacks share the same coordinate frame).
    Y = single(imageStack);
    [d1,d2,~] = size(Y);
    options_rigid = NoRMCorreSetParms('d1',d1,'d2',d2,'bin_width',35,'max_shift',4,'us_fac',50,'init_batch',1);
    [M_mc, ~] = normcorre(Y, options_rigid, globalTemplate);
    imageStack_mc = M_mc; % single

    % Save motion-corrected stack as multi-page TIFF in mcTiffFolder
    [~, rawName, ~] = fileparts(ImgfolderContents(imagesIndex(tt)).name);
    mcName = [rawName, '_mc.tif'];
    mcPath = fullfile(mcTiffFolder, mcName);

    % Convert to uint16 for saving (adjust if your data type differs)
    if ~isa(imageStack_mc, 'uint16')
        % simple scaling: assume original-like range fits into uint16
        mcMin = min(imageStack_mc(:));
        mcMax = max(imageStack_mc(:));
        if mcMax > mcMin
            imageStack_mc_uint16 = uint16( (imageStack_mc - mcMin) ./ (mcMax - mcMin) * double(intmax('uint16')) );
        else
            imageStack_mc_uint16 = uint16(zeros(size(imageStack_mc)));
        end
    else
        imageStack_mc_uint16 = imageStack_mc;
    end

    tOut = Tiff(mcPath, 'w');
    tagstruct.ImageLength = d1;
    tagstruct.ImageWidth  = d2;
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 16;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Compression = Tiff.Compression.LZW;

    for k = 1:size(imageStack_mc_uint16, 3)
        tOut.setTag(tagstruct);
        tOut.write(imageStack_mc_uint16(:,:,k));
        if k < size(imageStack_mc_uint16, 3)
            tOut.writeDirectory();
        end
    end
    tOut.close();

    nCap = min(size(imageStack_mc, 3), maxDvFrameCap);
    imageStackMean_mc = mean(imageStack_mc(:, :, 1:nCap), 3);
    if ~maxDvTrialMask(tt)
        maxDvStack(:, :, tt) = nan(size(imageStackMean_mc));
    else
        maxDvStack(:, :, tt) = imageStackMean_mc;
    end

end

meanMaxDvStack = mean(maxDvStack, 3, 'omitnan'); % Grand mean of motion-corrected trial means
meanFluorMaxDvStack = meanMaxDvStack;

voltMapping.mcTiffFolder           = mcTiffFolder;
voltMapping.maxDvStack         = maxDvStack;
voltMapping.meanFluorMaxDvStack = meanFluorMaxDvStack;
voltMapping.globalMcTemplate = globalTemplate;
voltMapping.maxDvTrialMask = maxDvTrialMask;
voltMapping.maxDvFrameCap = maxDvFrameCap;

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
        roiMeanMaxDvStack(roughRoiXAllCells{nn}(rr), roughRoiYAllCells{nn}(rr)) = mean(maxDvStack(roughRoiXAllCells{nn}(rr), roughRoiYAllCells{nn}(rr),:), 3);
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

%% Preallocate containers for per-trial fine ROIs
fineRoiXAllCells = cell(nCells, 1);
fineRoiYAllCells = cell(nCells, 1);
for nn = 1:nCells
    fineRoiXAllCells{nn} = cell(length(imagesIndex), 1);
    fineRoiYAllCells{nn} = cell(length(imagesIndex), 1);
end

%% F, F0, dF, dF/F0 Calculation using motion-corrected stacks and per-trial fine ROIs
input('This step uses motion-corrected movies, computes trial-specific fine ROIs, and calculates dF/F (ctrl+c to stop!)');

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
    dFFCellNames{nn}                   = ['dFF', 'cell', num2str(nn)];
    holoSortedImagingCellNames{nn}     = ['holoSortedImagingAllTrials_', 'cell', num2str(nn)];
    filtHoloSortedImagingCellNames{nn} = ['filtHoloSortedImagingAllTrials_', 'cell', num2str(nn)];

    analysisStruct.(roiMeanFCellNames{nn})              = [];
    analysisStruct.(bkgrndMeanFCellNames{nn})           = [];
    analysisStruct.(subScalarCellNames{nn})             = [];
    analysisStruct.(roiMeanFCorrectedCellNames{nn})     = [];
    analysisStruct.(globalF0CellNames{nn})              = [];
    analysisStruct.(dFCellNames{nn})                    = [];
    analysisStruct.(dFFCellNames{nn})                   = [];
    analysisStruct.(F0CellNames{nn})                    = cell(nConds, 1);
    analysisStruct.(holoSortedImagingCellNames{nn})     = cell(nConds, 1);
    analysisStruct.(filtHoloSortedImagingCellNames{nn}) = cell(nConds, 1);
    
    for cc = 1:nConds
        analysisStruct.(F0CellNames{nn}){cc}                    = cell(nHolos(cc), 1);
        analysisStruct.(holoSortedImagingCellNames{nn}){cc}     = cell(nHolos(cc), 1);
        analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc} = cell(nHolos(cc), 1);
    end
end

counter = 0;
for tt = 1:nTrials
    counter = counter+1;
    disp(['Trial number: ', num2str(counter)]);

    % Load motion-corrected TIFF for this trial
    rawName = ImgfolderContents(imagesIndex(tt)).name;
    [~, baseName, ~] = fileparts(rawName);
    mcName = [baseName, '_mc.tif'];
    currImgPath = fullfile(mcTiffFolder, mcName);

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
    [H,W] = size(firstFrame);
    
    imageStack = zeros(H, W, numFrames, 'like', firstFrame);
    imageStack(:,:,1) = firstFrame;
    
    for k = 2:numFrames
        t.setDirectory(k);
        imageStack(:,:,k) = t.read();
    end
    t.close();

    % Mean image for this trial (motion-corrected)
    meanImgThisTrial = mean(single(imageStack), 3);

    % Start time of first hologram stimulation in image sample
    startTimeImaging = floor(startTime*imagingFreq);    

    %--------------------------------------------------------------
    % Per-trial fine ROIs for all cells
    %--------------------------------------------------------------
    counter2 = 0;
    for nn = 1:nCells
        counter2 = counter2+1;
        disp(['Computing trial-specific ROI for Cell: ', num2str(counter2)]);

        %--------------------------------------------------------------
        % Compute per-trial fine ROI for this cell, restricted to rough ROI
        % using same gaussian + fibermetric logic as above.
        %--------------------------------------------------------------
        if isempty(fineRoiXAllCells{nn}{tt})
            roiMeanThisTrial = zeros(size(meanImgThisTrial));
            for rr = 1:length(roughRoiXAllCells{nn})
                roiMeanThisTrial(roughRoiXAllCells{nn}(rr), roughRoiYAllCells{nn}(rr)) = ...
                    meanImgThisTrial(roughRoiXAllCells{nn}(rr), roughRoiYAllCells{nn}(rr));
            end

            % roiMaxStackDouble = im2double(roiMeanThisTrial);
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
            
            % roiMaxStackRidge = fibermetric(roiMaxStackNorm, 'StructureSensitivity', 2);
            
            % Apply fibermetric only within the rough ROI pixels, not the whole image
            roiMaxStackRidge = zeros(size(roiMaxStackNorm));
            % Create a mask for the rough ROI pixels
            roiMask = false(size(roiMaxStackNorm));
            for rr = 1:length(roughRoiXAllCells{nn})
                roiMask(roughRoiXAllCells{nn}(rr), roughRoiYAllCells{nn}(rr)) = true;
            end
            % Extract a tight bounding box around the rough ROI to speed processing
            [rows, cols] = find(roiMask);
            rmin = max(min(rows)-2, 1); rmax = min(max(rows)+2, size(roiMask,1));
            cmin = max(min(cols)-2, 1); cmax = min(max(cols)+2, size(roiMask,2));
            subImg = roiMaxStackNorm(rmin:rmax, cmin:cmax);
            subMask = roiMask(rmin:rmax, cmin:cmax);
            % Apply fibermetric to the sub-image, then zero outside subMask
            subRidge = fibermetric(subImg, 'StructureSensitivity', 2);
            subRidge(~subMask) = 0;
            % Place result back into full-size image
            roiMaxStackRidge(rmin:rmax, cmin:cmax) = subRidge;

            valsR = nonzeros(roiMaxStackRidge);
            if ~isempty(valsR)
                thr  = prctile(valsR, 60);
            else
                thr = 0;
            end
            roiMaxStackRidgeReduced = roiMaxStackRidge;
            roiMaxStackRidgeReduced(roiMaxStackRidgeReduced < thr) = 0;
            roiMaxStackRidgeReduced(roiMaxStackRidgeReduced > 0) = 1;

            [roiFineX, roiFineY] = find(roiMaxStackRidgeReduced);
            if isempty(roiFineX)
                % fall back to rough ROI if fibermetric fails
                roiFineX = roughRoiXAllCells{nn};
                roiFineY = roughRoiYAllCells{nn};
            end

            fineRoiXAllCells{nn}{tt} = roiFineX;
            fineRoiYAllCells{nn}{tt} = roiFineY;
        end
    end

    %--------------------------------------------------------------
    % Pass 2: extract F traces; neuropil from global ring only (same pixels every trial)
    %--------------------------------------------------------------
    counter2 = 0;
    for nn = 1:nCells
        counter2 = counter2+1;
        disp(['Getting F traces and F0 for Cell: ', num2str(counter2)]);

        % Use fine, trial-specific ROI pixels for F extraction
        rawWholeRoiF = imageStack(fineRoiXAllCells{nn}{tt}, fineRoiYAllCells{nn}{tt}, :);
        roiMeanF = [];
        for ff = 1:size(rawWholeRoiF, 3)
            roiMeanF(ff, 1) = mean(mean(rawWholeRoiF(:, :, ff)));
        end
        
        % Global neuropil ring (mean max-dV fibermetric soma + dilated annulus); not redefined per trial.
        bkgrndRoiXUse = bkgrndRoiXAllCells{nn};
        bkgrndRoiYUse = bkgrndRoiYAllCells{nn};
        bkgrndRoiXAllCells_trial{nn}{tt} = bkgrndRoiXUse;
        bkgrndRoiYAllCells_trial{nn}{tt} = bkgrndRoiYUse;

        if isempty(bkgrndRoiXUse)
            bkgrndMeanF = zeros(numFrames, 1);
        else
            rawWholeBkgrndF = imageStack(bkgrndRoiXUse, bkgrndRoiYUse, :);
            bkgrndMeanF = [];
            for ff = 1:size(rawWholeBkgrndF, 3)
                bkgrndMeanF(ff, 1) = mean(mean(rawWholeBkgrndF(:, :, ff)));
            end
        end
        
        % Neuropil correction subscalar (same logic as original)
        baselineIndices = 1:startTimeImaging;
        if ~isempty(bkgrndRoiXUse)
            bFit = robustfit(bkgrndMeanF(baselineIndices), roiMeanF(baselineIndices));
            alphaScalar = bFit(2);
            alphaScalar = min(max(alphaScalar, 0), 1);
            if alphaScalar > 0.8
                alphaScalar = 0.8;
            end
        end
        alphaScalar = 0.90;

        roiMeanFCorrected = roiMeanF - alphaScalar*bkgrndMeanF;

        % Break apart the imageStack for this trial into stim windows and rearrange
        if isempty((voltMapping.outParams.sequenceThisTrial{tt}))
            voltMapping.outParams.sequenceThisTrial{tt} = zeroDummySequence;
        end

        cutOffFreqIm = 40;
        [bIm, aIm] = butter(4, cutOffFreqIm/(imagingFreq/2));        

        holoSeqThisTrial = (unique(voltMapping.outParams.sequenceThisTrial{tt}, 'stable') - min(unique(voltMapping.outParams.sequenceThisTrial{tt}, 'stable')) + 1)';

        for hh = 1:nHolos(voltMapping.trialCond(tt, 1))
            if isempty(voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)})
                voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)} = voltMapping.outParams.firstStimTimes{1, 2};
            end
            
            roiFCorrectedThisHolo = roiMeanFCorrected(...
                floor((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*imagingFreq):(ceil((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*imagingFreq)+ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)));

            roiFCorrectedThisHoloPreStim = roiFCorrectedThisHolo(1:(preStimWindow/1000*imagingFreq)-1);
            f0ThisHolo = mean(roiFCorrectedThisHoloPreStim);
            
            dFThisHolo = roiFCorrectedThisHolo - f0ThisHolo;
            dFFThisHolo = dFThisHolo/f0ThisHolo;

            if UpOrDown == '2'
                dFFThisHolo = -dFFThisHolo;
            elseif UpOrDown =='1'
                dFFThisHolo = dFFThisHolo;
            end
            
            filtdffThisHolo = filter(bIm, aIm, dFFThisHolo);

            if ismember(tt, excludeTrials)
                analysisStruct.(F0CellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}                    = [analysisStruct.(F0CellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, NaN];
                analysisStruct.(holoSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}     = [analysisStruct.(holoSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, NaN(ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)+2, 1)];
                analysisStruct.(filtHoloSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(filtHoloSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, NaN(ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)+2, 1)];
            else
                analysisStruct.(F0CellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}                    = [analysisStruct.(F0CellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, f0ThisHolo];
                analysisStruct.(holoSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}     = [analysisStruct.(holoSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, dFFThisHolo];
                analysisStruct.(filtHoloSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(filtHoloSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, filtdffThisHolo];
            end        
        end

        if ismember(tt, excludeTrials)
            analysisStruct.(roiMeanFCellNames{nn})(:, tt)              = NaN(numFrames, 1);
            analysisStruct.(bkgrndMeanFCellNames{nn})(:, tt)           = NaN(numFrames, 1);
            analysisStruct.(subScalarCellNames{nn})(tt, 1)             = NaN;
            analysisStruct.(roiMeanFCorrectedCellNames{nn})(:, tt)     = NaN(numFrames, 1);
        else
            analysisStruct.(roiMeanFCellNames{nn})(:, tt)              = roiMeanF;
            analysisStruct.(bkgrndMeanFCellNames{nn})(:, tt)           = bkgrndMeanF;
            analysisStruct.(subScalarCellNames{nn})(tt, 1)             = alphaScalar;
            analysisStruct.(roiMeanFCorrectedCellNames{nn})(:, tt)     = roiMeanFCorrected;
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
        plot(linspace(0, size(exclFiltHoloSortedImagingMean{cc}{hh}, 1)/imagingFreq, size(exclFiltHoloSortedImagingMean{cc}{hh}, 1)), exclFiltHoloSortedImagingMean{cc}{hh}*100, '-', 'linewidth', 2, 'color', 'g');
        ylabel('dF/F (%)');
        xlabel('Time (s)');
        %         ylim([-1 max(holoSortedImagingMean{cc}{hh})*100])
        %         ylim([min(holoSortedImagingMean{cc}{hh})*100 max(holoSortedImagingMean{cc}{hh})*100])
        % ylim([- 2])
        
        % if cc > 2
            ylim([-5 30])
        % end

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

        for nn = 1:length(nPulseCoords)
            xline(nPulseCoords(nn)/Fs, '-', 'LineWidth', 5, 'color', [1 0 0 0.1]);
        end
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
yScales    = [100, 10, 100, 10];
subTitles  = {'Mean', 'Filt mean', 'Excl mean', 'Excl filt mean'};

for cc = 1:nConds
    for hh = 1:nHolos(cc)
        figure(20000 + cc*1000 + hh);
        set(gcf, 'Position', [100, 100, 2000, 400]);
        clf

        nT = numel(holoSortedImagingMean{cc}{hh});
        tAxis = linspace(0, nT/imagingFreq, nT);

        for sp = 1:4
            subplot(1, 4, sp);
            hold on
            m = meanSeries{sp}{cc}{hh};
            ci = ciSeries{sp}{cc}{hh};
            ys = yScales(sp);
            fill([tAxis, fliplr(tAxis)], [ci(:, 1)'*ys, fliplr(ci(:, 2)'*ys)], [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
            plot(tAxis, ci(:, 1)*ys, '--', 'linewidth', 1, 'color', [0.7, 0.7, 0.7]);
            plot(tAxis, ci(:, 2)*ys, '--', 'linewidth', 1, 'color', [0.7, 0.7, 0.7]);
            plot(tAxis, m*ys, '-', 'linewidth', 2, 'color', 'g');
            ylabel('dF/F (%)');
            xlabel('Time (s)');
            title(subTitles{sp});
            for pulseIdx = 1:length(nPulseCoords)
                xline(nPulseCoords(pulseIdx)/Fs, '-', 'LineWidth', 5, 'color', [1, 0, 0, 0.1]);
            end
            hold off
        end
        sgtitle(sprintf('Cell %d — cond %d, holo %d (4-panel)', cellIdx, cc, hh));
        pause
    end
end

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
        sgtitle({['Neuropil sanity (global ring, identical every trial): ', num2str(nCells), ' cells; color = cell index'], ...
            ['every ', num2str(sanityTrialStep), ' trials (1:', num2str(sanityTrialStep), ':', num2str(nTrials), ')']});
    end
end

%% Save Analysis Results
expID = num2str(mouseID);
voltMapping.mouseID = ['voltMapping_Analysis_', expID, '_MultiCellAnalysis_MCfineROI'];
directory = '/Users/masatosadahiro/Documents/Data/Voltage Imaging/Voltage Imaging/voltMapping/Analysis Results/Analysis_newMCanddFFcalc';
fileName = [num2str(voltMapping.mouseID), '.mat'];

% Persist imaging traces/stats plus field-name cell arrays (needed to index analysisStruct after load).
varsToSave = [{'voltMapping', 'analysisStruct', ...
    'F0CellNames', 'roiMeanFCellNames', 'bkgrndMeanFCellNames', 'subScalarCellNames', ...
    'roiMeanFCorrectedCellNames', 'globalF0CellNames', 'dFCellNames', 'dFFCellNames', ...
    'holoSortedImagingCellNames', 'filtHoloSortedImagingCellNames', ...
    'holoSortedMeanCellNames', 'filtHoloSortedMeanCellNames', ...
    'CIDffAllCondsCellNames', 'filtCIDffAllCondsCellNames', ...
    'mouseID', 'ePhysAvail'}, ...
    {'nCells', 'nTrials', 'nConds', 'nHolos', 'imagingFreq', 'Fs', 'ipi', 'nPulses', ...
    'preStimWindow', 'postStimWindow', 'startTime', 'UpOrDown', 'excludeTrials'}];
if exist('zeroDummySequence', 'var')
    varsToSave{end+1} = 'zeroDummySequence';
end
save(fullfile(directory, fileName), varsToSave{:}, '-v7.3');

TimeNow = datetime('now', 'Format', 'dd-MMM-yyyy HH:mm:ss');
disp(['finished saving at: ' char(TimeNow)])

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
nHolos = voltMapping.nHolos; % number of holograms in grid
pulseDurs = unique(voltMapping.pulseDur);
nPulses = unique(voltMapping.nPulses);
ipi = voltMapping.ipi;
nextHoloDelay = voltMapping.nextHoloDelay;
% SpotCoordinates = SortedData.holoRequest.targets;
startTime = (voltMapping.holoStimParams.startTime)/1000;
imagesIndex = voltMapping.imagesIndex;
UpOrDown = voltMapping.UpOrDown;
ephysFilePath = voltMapping.ephysFilePath;
ImgsFilePath = voltMapping.ImgsFilePath;
cellID = voltMapping.cellID;
