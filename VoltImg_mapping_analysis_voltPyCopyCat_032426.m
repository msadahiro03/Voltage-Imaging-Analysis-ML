%% Voltage Imaging Mapping Analysis — VolPy-style pipeline (no spike pursuit)
% 1) Per-trial NoRMCorre (unchanged); motion-corrected TIFFs saved as before.
% 2) Summary images for segmentation: z-scored mean + local correlation image
%    on a high-passed concatenated movie (VolPy, PLOS Comp Biol 2021).
% 3) Segmentation: load neuron masks from .mat (e.g. Mask R-CNN / CaImAn VolPy
%    output) or draw freehand on the summary RGB; manual checkpoint before traces.
% 4) Per-trial traces: context region, high-pass (bleaching), SVD of local
%    background + ridge subtraction, ridge spatial filter, optional temporal
%    smoothing, low-pass subthreshold band (~20 Hz default).
% 5) Hologram windows / dF/F / trial exclusion: same structure as the MCfineROI
%    script, operating on the VolPy-style denoised whole-trial trace.

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
rawImgNChannels = input('1 if raw TIFF is single-channel (use every page); 2 if two-color interleaved (green = TIFF pages 1,3,5,...): ');

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
voltMapping.mouseID = ['voltMapping_Analysis_', expID, '_VolPyStyle'];
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

% Preallocate maxDvStack using first TIFF dimensions
firstImgPath = fullfile(ImgfolderContents(imagesIndex(1)).folder, ImgfolderContents(imagesIndex(1)).name);
infoFirst = imfinfo(firstImgPath);
maxDvStack = zeros(infoFirst(1).Height, infoFirst(1).Width, length(imagesIndex));

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

maskCounter = 0;
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

        maskCounter = maskCounter + 1;
        imageStackMean_mc = mean(imageStack_mc, 3);
        maxDvStack(:, :, maskCounter) = imageStackMean_mc;

end

% Drop unused preallocated planes if any
if maskCounter < size(maxDvStack,3)
    maxDvStack(:,:,maskCounter+1:end) = [];
end

meanMaxDvStack = mean(maxDvStack, 3); % Grand mean of motion-corrected trial means
meanFluorMaxDvStack = meanMaxDvStack;

voltMapping.mcTiffFolder        = mcTiffFolder;
voltMapping.maxDvStack          = maxDvStack;
voltMapping.meanFluorMaxDvStack = meanFluorMaxDvStack;
voltMapping.globalMcTemplate    = globalTemplate;

%% VolPy-style summary images (mean z-score + local correlation on high-pass movie)
input('Building concatenated movie for VolPy summary images (continue or ctrl+c to stop!)');

[H0, W0] = size(meanFluorMaxDvStack);
% VolPy summary: ~30k frames = (frames per trial) × (X evenly spaced trials).
% All MC TIFFs assumed equal length; first *_mc.tif in folder sets frame count.
targetFramesPerTrialVolpy = 200;
maxFramesVolpySummary = 30000;

mcDirList = dir(fullfile(mcTiffFolder, '*_mc.tif'));
if isempty(mcDirList)
    error('VolPy summary: no *_mc.tif files found in mcTiffFolder.');
end
[~, sortMc] = sort({mcDirList.name});
mcDirList = mcDirList(sortMc);
firstMcPath = fullfile(mcTiffFolder, mcDirList(1).name);
tifProbe = Tiff(firstMcPath, 'r');
nProbe = 1;
while true
    try
        tifProbe.setDirectory(nProbe);
        nProbe = nProbe + 1;
    catch
        nProbe = nProbe - 1;
        break
    end
end
tifProbe.close();
nFramesPerMcTiff = nProbe;
if nFramesPerMcTiff < 1
    error('VolPy summary: could not read frames from %s', firstMcPath);
end

if nFramesPerMcTiff >= targetFramesPerTrialVolpy
    framesPerTrialVolpyActual = targetFramesPerTrialVolpy;
else
    framesPerTrialVolpyActual = nFramesPerMcTiff;
end

nTrialsVolpyTarget = floor(maxFramesVolpySummary / framesPerTrialVolpyActual);

% Pool: all imaging trials not excluded (same order as imagesIndex)
poolTT = [];
for tt = 1:length(imagesIndex)
    if ~ismember(imagesIndex(tt), excludeTrials)
        poolTT(end + 1) = tt; 
    end
end
N = numel(poolTT);
if N == 0
    error('VolPy summary: no trials after exclusions.');
end

nPick = min(nTrialsVolpyTarget, N);
idxPick = round(linspace(1, N, nPick));
idxPick = unique(idxPick, 'stable');
selectedTT = poolTT(idxPick);

disp(['VolPy summary: first MC = ', mcDirList(1).name, ' → ', num2str(nFramesPerMcTiff), ' frames/trial; ', ...
    'using ', num2str(framesPerTrialVolpyActual), ' frames × ', num2str(numel(selectedTT)), ...
    ' evenly spaced trials (of ', num2str(N), ' in pool).']);

Ycat = [];
for ii = 1:numel(selectedTT)
    tt = selectedTT(ii);
    disp(['Loading summary chunk ', num2str(ii), ' / ', num2str(numel(selectedTT)), ...
        ' (trial index ', num2str(tt), ')']);
    rawName = ImgfolderContents(imagesIndex(tt)).name;
    [~, baseName, ~] = fileparts(rawName);
    mcPathTry = fullfile(mcTiffFolder, [baseName, '_mc.tif']);
    if exist(mcPathTry, 'file') ~= 2
        error('VolPy summary: missing MC file for trial index %d: %s', tt, mcPathTry);
    end
    tif = Tiff(mcPathTry, 'r');
    tif.setDirectory(1);
    fr0 = tif.read();
    Ychunk = zeros(size(fr0, 1), size(fr0, 2), framesPerTrialVolpyActual, 'like', single(fr0));
    Ychunk(:, :, 1) = single(fr0);
    for kf = 2:framesPerTrialVolpyActual
        tif.setDirectory(kf);
        Ychunk(:, :, kf) = single(tif.read());
    end
    tif.close();
    if UpOrDown == '2'
        Ychunk = -Ychunk;
    end
    Ycat = cat(3, Ycat, Ychunk);
end

if isempty(Ycat)
    error('No frames collected for VolPy summary images (check mcTiffFolder / exclusions).');
end

[meanImgZ, corrImgZ, ~] = voltImg_volpySummaryImages(Ycat, imagingFreq);
volpyRgbBase = fullfile(saveDirectory, [num2str(mouseID), '_volpy_summary']);
voltImg_volpySaveMaskRcnnRgb(meanImgZ, corrImgZ, volpyRgbBase);

voltMapping.volpy.meanImgZ = meanImgZ;
voltMapping.volpy.corrImgZ = corrImgZ;
voltMapping.volpy.maxFramesUsed = size(Ycat, 3);
voltMapping.volpy.targetFramesPerTrialVolpy = targetFramesPerTrialVolpy;
voltMapping.volpy.framesPerTrialVolpyActual = framesPerTrialVolpyActual;
voltMapping.volpy.nFramesPerMcTiffFromProbe = nFramesPerMcTiff;
voltMapping.volpy.firstMcTiffProbed = firstMcPath;
voltMapping.volpy.maxFramesVolpySummary = maxFramesVolpySummary;
voltMapping.volpy.nTrialsVolpyTarget = nTrialsVolpyTarget;
voltMapping.volpy.nTrialsVolpySelected = numel(selectedTT);
voltMapping.volpy.nTrialsInPoolVolpy = N;
voltMapping.volpy.selectedTrialIndices = imagesIndex(selectedTT);
voltMapping.volpy.rgbExportPng = [volpyRgbBase, '_volpy_maskrcnn_rgb.png'];
% Mask R-CNN (VolPy): run inference in Python/CaImAn on the saved RGB PNG, then
% save a .mat with cellMasks{{1},{2},...} each HxW logical matching this FOV.
% See https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008806
% and https://github.com/flatironinstitute/CaImAn/wiki/Training-Mask-R-CNN

%% Segmentation: load masks (Mask R-CNN / external) or draw manually + checkpoint
input('Segmentation: load masks or draw on summary RGB (continue or ctrl+c to stop!)');

disp('--- VolPy-style segmentation ---');
disp('1 = Load binary masks from a .MAT file');
disp('    Expect cellMasks: 1xN or Nx1 cell of HxW logical, OR masks: HxWxN logical');
disp('2 = Draw freehand ROIs on mean fluorescence (MC trial means, winter colormap)');
segChoice = input('Enter 1 or 2: ', 's');

roiMaskLogical = {};
if segChoice == '1'
    [maskMatFile, maskMatDir] = uigetfile('*.mat', 'Select MAT with cellMasks or masks');
    if isequal(maskMatFile, 0)
        error('Mask file selection canceled.');
    end
    S = load(fullfile(maskMatDir, maskMatFile));
    if isfield(S, 'cellMasks')
        cm = S.cellMasks;
        if ~iscell(cm)
            error('cellMasks must be a cell array of HxW logical masks.');
        end
        nCells = numel(cm);
        for nn = 1:nCells
            roiMaskLogical{nn} = logical(cm{nn});
            if ~isequal(size(roiMaskLogical{nn}), [H0, W0])
                error('Mask %d size [%d,%d] does not match FOV [%d,%d].', nn, ...
                    size(roiMaskLogical{nn}, 1), size(roiMaskLogical{nn}, 2), H0, W0);
            end
        end
    elseif isfield(S, 'masks')
        M = S.masks;
        szM = size(M);
        if numel(szM) == 2
            nCells = 1;
            roiMaskLogical{1} = logical(M);
        else
            nCells = szM(3);
            for nn = 1:nCells
                roiMaskLogical{nn} = logical(M(:, :, nn));
            end
        end
        for nn = 1:nCells
            if ~isequal(size(roiMaskLogical{nn}), [H0, W0])
                error('Mask %d size does not match FOV.', nn);
            end
        end
    else
        error('MAT file must contain cellMasks or masks.');
    end
else
    % Preview (same spirit as original rough-ROI section): inspect FOV before choosing count
    figure(9);
    set(gcf, 'Position', [100, 100, 1800, 900]);
    clf
    imagesc(meanFluorMaxDvStack);
    axis equal;
    axis image;
    colormap(gca, winter);
    colorbar;
    set(gca, 'FontSize', 12);
    title('Preview: mean fluorescence (motion-corrected trial means) — then enter how many ROIs to draw');
    nCells = input('How many neurons to draw?: ', 's');
    nCells = str2double(nCells);
    % Familiar mean-fluorescence view (same as rough ROI elsewhere); not the z-scored VolPy RGB.
    for nn = 1:nCells
        f1 = figure(10);
        set(f1, 'Position', [100, 100, 1600, 800]);
        clf(f1);
        ax = axes('Parent', f1);
        imagesc(ax, meanFluorMaxDvStack);
        axis(ax, 'image');
        colormap(ax, winter);
        colorbar(ax);
        set(ax, 'FontSize', 12);
        title(ax, sprintf(['Cell %d / %d: draw freehand ROI on mean fluorescence ', ...
            '(motion-corrected trial means; finish ROI then close or proceed)'], nn, nCells));
        roiHandSelect = drawfreehand(ax);
        roiMaskLogical{nn} = logical(createMask(roiHandSelect));
        if ~isequal(size(roiMaskLogical{nn}), [H0, W0])
            error('Manual ROI %d: mask size [%d,%d] does not match FOV [%d,%d].', nn, ...
                size(roiMaskLogical{nn}, 1), size(roiMaskLogical{nn}, 2), H0, W0);
        end
        close(f1);
    end
end

% Ring masks for sanity plots (same spirit as original neuropil ring)
innerBuffer = 2;
ringWidth = 3;
minArea = 50;
roiXAllCells_global = cell(nCells, 1);
roiYAllCells_global = cell(nCells, 1);
bkgrndRoiXAllCells = cell(nCells, 1);
bkgrndRoiYAllCells = cell(nCells, 1);
roughRoiXAllCells = cell(nCells, 1);
roughRoiYAllCells = cell(nCells, 1);
allRois = zeros(H0, W0);
for nn = 1:nCells
    R = roiMaskLogical{nn};
    [rx, ry] = find(R);
    roiXAllCells_global{nn} = rx;
    roiYAllCells_global{nn} = ry;
    roughRoiXAllCells{nn} = rx;
    roughRoiYAllCells{nn} = ry;
    allRois = allRois + double(R) * nn;

    innerSelect = imdilate(R, strel('disk', innerBuffer));
    outerSelect = imdilate(R, strel('disk', innerBuffer + ringWidth));
    backgroundRing = outerSelect & ~innerSelect;
    meanD = im2double(meanFluorMaxDvStack);
    valsBk = meanD(backgroundRing);
    if ~isempty(valsBk)
        brightCut = prctile(valsBk, 95);
        ringClean = backgroundRing & (meanD <= brightCut);
    else
        ringClean = backgroundRing;
    end
    ringClean = bwareaopen(ringClean, 7);
    if nnz(ringClean) < minArea
        ringClean = backgroundRing;
    end
    [bkgrndRoiXAllCells{nn}, bkgrndRoiYAllCells{nn}] = find(ringClean);
end

figure(25); clf; set(gcf, 'Position', [80, 80, 1600, 700]);
tiledlayout(1, 2);
nexttile; imagesc(meanImgZ); axis image; colormap(turbo); colorbar;
title('Mean (z-score)');
hold on;
for nn = 1:nCells
    boundary = bwboundaries(roiMaskLogical{nn}, 'noholes');
    if ~isempty(boundary)
        b = boundary{1};
        plot(b(:, 2), b(:, 1), 'w-', 'LineWidth', 1.2);
    end
end
hold off
nexttile; imagesc(corrImgZ); axis image; colormap(turbo); colorbar;
title('Local correlation (z-score)');
hold on;
for nn = 1:nCells
    boundary = bwboundaries(roiMaskLogical{nn}, 'noholes');
    if ~isempty(boundary)
        b = boundary{1};
        plot(b(:, 2), b(:, 1), 'w-', 'LineWidth', 1.2);
    end
end
hold off
sgtitle('Manual checkpoint: confirm ROIs (save new MAT and re-run if needed)');

input('Press Enter after inspecting figure 25 (ROIs overlaid on summary images).');

voltMapping.nCells = nCells;
voltMapping.volpy.roiMaskLogical = roiMaskLogical;
voltMapping.roughRoiXAllCells = roughRoiXAllCells;
voltMapping.roughRoiYAllCells = roughRoiYAllCells;
voltMapping.roiXAllCells_global = roiXAllCells_global;
voltMapping.roiYAllCells_global = roiYAllCells_global;
voltMapping.bkgrndRoiXAllCells = bkgrndRoiXAllCells;
voltMapping.bkgrndRoiYAllCells = bkgrndRoiYAllCells;
voltMapping.allRois = allRois;

% VolPy trace parameters (edit as needed)
volpyTraceParams = struct();
volpyTraceParams.fps = imagingFreq;
volpyTraceParams.contextDilatePx = 15;
volpyTraceParams.bgMinDistPx = 12;
volpyTraceParams.nPc = 8;
volpyTraceParams.lambdaB = 0.01;
volpyTraceParams.lambdaW = 0.01;
volpyTraceParams.fcHpBleach = 1/3;
volpyTraceParams.fcLpSubthresh = 20;
volpyTraceParams.reversePolarity = (UpOrDown == '2');
volpyTraceParams.temporalSmoothFrames = 0;
voltMapping.volpy.traceParams = volpyTraceParams;

%% Preallocate per-trial ROI indices (global masks replicated for each trial index)
fineRoiXAllCells = cell(nCells, 1);
fineRoiYAllCells = cell(nCells, 1);
bkgrndRoiXAllCells_trial = cell(nCells, 1);
bkgrndRoiYAllCells_trial = cell(nCells, 1);
for nn = 1:nCells
    fineRoiXAllCells{nn} = cell(length(imagesIndex), 1);
    fineRoiYAllCells{nn} = cell(length(imagesIndex), 1);
    bkgrndRoiXAllCells_trial{nn} = cell(length(imagesIndex), 1);
    bkgrndRoiYAllCells_trial{nn} = cell(length(imagesIndex), 1);
    [gx, gy] = find(roiMaskLogical{nn});
    for tt = 1:length(imagesIndex)
        fineRoiXAllCells{nn}{tt} = gx;
        fineRoiYAllCells{nn}{tt} = gy;
    end
end

%% VolPy-style traces + dF/F + hologram sorting (motion-corrected stacks, global ROIs)
input('VolPy-style trace extraction and dF/F / hologram windows (ctrl+c to stop!)');

% Preallocate data for each cell
holoSortedImagingCellNames = cell(nCells, 1);
filtHoloSortedImagingCellNames = cell(nCells, 1);
holoSortedVolpySubCellNames = cell(nCells, 1);
filtHoloSortedVolpySubCellNames = cell(nCells, 1);
volpySubthreshFullCellNames = cell(nCells, 1);
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
    holoSortedVolpySubCellNames{nn}    = ['holoSortedVolpySubAllTrials_', 'cell', num2str(nn)];
    filtHoloSortedVolpySubCellNames{nn} = ['filtHoloSortedVolpySubAllTrials_', 'cell', num2str(nn)];
    volpySubthreshFullCellNames{nn}   = ['volpySubthreshFull_', 'cell', num2str(nn)];

    analysisStruct.(roiMeanFCellNames{nn})              = [];
    analysisStruct.(bkgrndMeanFCellNames{nn})           = [];
    analysisStruct.(subScalarCellNames{nn})             = [];
    analysisStruct.(roiMeanFCorrectedCellNames{nn})     = [];
    analysisStruct.(globalF0CellNames{nn})              = [];
    analysisStruct.(dFCellNames{nn})                    = [];
    analysisStruct.(dFFCellNames{nn})                   = [];
    analysisStruct.(volpySubthreshFullCellNames{nn})    = [];
    analysisStruct.(F0CellNames{nn})                    = cell(nConds, 1);
    analysisStruct.(holoSortedImagingCellNames{nn})     = cell(nConds, 1);
    analysisStruct.(filtHoloSortedImagingCellNames{nn}) = cell(nConds, 1);
    analysisStruct.(holoSortedVolpySubCellNames{nn})    = cell(nConds, 1);
    analysisStruct.(filtHoloSortedVolpySubCellNames{nn}) = cell(nConds, 1);
    
    for cc = 1:nConds
        analysisStruct.(F0CellNames{nn}){cc}                    = cell(nHolos(cc), 1);
        analysisStruct.(holoSortedImagingCellNames{nn}){cc}     = cell(nHolos(cc), 1);
        analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc} = cell(nHolos(cc), 1);
        analysisStruct.(holoSortedVolpySubCellNames{nn}){cc}    = cell(nHolos(cc), 1);
        analysisStruct.(filtHoloSortedVolpySubCellNames{nn}){cc} = cell(nHolos(cc), 1);
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

    meanImgThisTrial = mean(single(imageStack), 3);
    startTimeImaging = floor(startTime * imagingFreq);
    meanImgThisTrialDouble = im2double(meanImgThisTrial);

    allTrialRoiMask = false(size(meanImgThisTrial));
    for nn = 1:nCells
        allTrialRoiMask = allTrialRoiMask | roiMaskLogical{nn};
    end

    Ysingle = single(imageStack);
    if UpOrDown == '2'
        Yvolpy = -Ysingle;
    else
        Yvolpy = Ysingle;
    end

    for nn = 1:nCells
        disp(['VolPy-style traces, cell ', num2str(nn), ' / ', num2str(nCells)]);

        rawWholeRoiF = Ysingle(fineRoiXAllCells{nn}{tt}, fineRoiYAllCells{nn}{tt}, :);
        roiMeanF = squeeze(mean(mean(rawWholeRoiF, 1, 'omitnan'), 2, 'omitnan'));
        roiMeanF = roiMeanF(:);

        innerBuffer = 2;
        ringWidth = 3;
        minArea = 50;
        roiMaskThisCell = roiMaskLogical{nn};
        innerSelect = imdilate(roiMaskThisCell, strel('disk', innerBuffer));
        outerSelect = imdilate(roiMaskThisCell, strel('disk', innerBuffer + ringWidth));
        backgroundRing = outerSelect & ~innerSelect & ~allTrialRoiMask;
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
            if ~isempty(bkgrndRoiXAllCells{nn})
                globalInd = sub2ind(size(ringGlobalMask), bkgrndRoiXAllCells{nn}, bkgrndRoiYAllCells{nn});
                ringGlobalMask(globalInd) = true;
            end
            ringClean = ringGlobalMask & ~allTrialRoiMask;
        end
        [bkgrndRoiXTrial, bkgrndRoiYTrial] = find(ringClean);
        bkgrndRoiXAllCells_trial{nn}{tt} = bkgrndRoiXTrial;
        bkgrndRoiYAllCells_trial{nn}{tt} = bkgrndRoiYTrial;
        rawWholeBkgrndF = Ysingle(bkgrndRoiXTrial, bkgrndRoiYTrial, :);
        bkgrndMeanF = squeeze(mean(mean(rawWholeBkgrndF, 1, 'omitnan'), 2, 'omitnan'));
        bkgrndMeanF = bkgrndMeanF(:);

        baselineIndices = 1:min(max(1, startTimeImaging), numel(bkgrndMeanF));
        if numel(baselineIndices) < 3
            baselineIndices = 1:min(30, numel(bkgrndMeanF));
        end
        bFit = robustfit(bkgrndMeanF(baselineIndices), roiMeanF(baselineIndices));
        alphaScalar = bFit(2);
        alphaScalar = min(max(alphaScalar, 0), 1);
        if alphaScalar > 0.8
            alphaScalar = 0.8;
        end

        V = voltImg_volpyTraceExtractNoSpikes(Yvolpy, roiMaskLogical{nn}, ...
            'fps', volpyTraceParams.fps, ...
            'contextDilatePx', volpyTraceParams.contextDilatePx, ...
            'bgMinDistPx', volpyTraceParams.bgMinDistPx, ...
            'nPc', volpyTraceParams.nPc, ...
            'lambdaB', volpyTraceParams.lambdaB, ...
            'lambdaW', volpyTraceParams.lambdaW, ...
            'fcHpBleach', volpyTraceParams.fcHpBleach, ...
            'fcLpSubthresh', volpyTraceParams.fcLpSubthresh, ...
            'reversePolarity', false, ...
            'temporalSmoothFrames', volpyTraceParams.temporalSmoothFrames);

        roiMeanFCorrected = double(V.tSpatial(:));
        volpySubFull = double(V.tSubthresh(:));

        if isempty((voltMapping.outParams.sequenceThisTrial{tt}))
            voltMapping.outParams.sequenceThisTrial{tt} = zeroDummySequence;
        end

        cutOffFreqIm = 40;
        [bIm, aIm] = butter(4, cutOffFreqIm / (imagingFreq / 2));

        holoSeqThisTrial = (unique(voltMapping.outParams.sequenceThisTrial{tt}, 'stable') - min(unique(voltMapping.outParams.sequenceThisTrial{tt}, 'stable')) + 1)';

        for hh = 1:nHolos(voltMapping.trialCond(tt, 1))
            if isempty(voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)})
                voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)} = voltMapping.outParams.firstStimTimes{1, 2};
            end

            i0 = floor((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh) - preStimWindow / 1000) * imagingFreq);
            i1 = ceil((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh) - preStimWindow / 1000) * imagingFreq) + ceil((ipi * nPulses + (preStimWindow + postStimWindow)) / 1000 * imagingFreq);
            i0 = max(1, i0);
            i1 = min(numFrames, i1);
            if i1 < i0
                i1 = i0;
            end
            roiFCorrectedThisHolo = roiMeanFCorrected(i0:i1);
            volpySubThisHolo = volpySubFull(i0:i1);

            nPre = max(1, round((preStimWindow / 1000 * imagingFreq) - 1));
            roiFCorrectedThisHoloPreStim = roiFCorrectedThisHolo(1:min(nPre, numel(roiFCorrectedThisHolo)));
            f0ThisHolo = mean(roiFCorrectedThisHoloPreStim, 'omitnan');
            dFThisHolo = roiFCorrectedThisHolo - f0ThisHolo;
            dFFThisHolo = dFThisHolo / max(abs(f0ThisHolo), eps);

            volpySubPre = volpySubThisHolo(1:min(nPre, numel(volpySubThisHolo)));
            f0Sub = mean(volpySubPre, 'omitnan');
            dFFSubThisHolo = (volpySubThisHolo - f0Sub) / max(abs(f0Sub), eps);
            % 
            % if UpOrDown == '2'
            %     dFFThisHolo = -dFFThisHolo;
            %     dFFSubThisHolo = -dFFSubThisHolo;
            % end

            filtdffThisHolo = filter(bIm, aIm, dFFThisHolo);
            filtdffSubThisHolo = filter(bIm, aIm, dFFSubThisHolo);

            nanPad = NaN(ceil((ipi * nPulses + (preStimWindow + postStimWindow)) / 1000 * imagingFreq) + 2, 1);

            if ismember(tt, excludeTrials)
                analysisStruct.(F0CellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(F0CellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, NaN];
                analysisStruct.(holoSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(holoSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, nanPad];
                analysisStruct.(filtHoloSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(filtHoloSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, nanPad];
                analysisStruct.(holoSortedVolpySubCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(holoSortedVolpySubCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, nanPad];
                analysisStruct.(filtHoloSortedVolpySubCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(filtHoloSortedVolpySubCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, nanPad];
            else
                analysisStruct.(F0CellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(F0CellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, f0ThisHolo];
                analysisStruct.(holoSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(holoSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, dFFThisHolo];
                analysisStruct.(filtHoloSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(filtHoloSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, filtdffThisHolo];
                analysisStruct.(holoSortedVolpySubCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(holoSortedVolpySubCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, dFFSubThisHolo];
                analysisStruct.(filtHoloSortedVolpySubCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(filtHoloSortedVolpySubCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, filtdffSubThisHolo];
            end
        end

        if ismember(tt, excludeTrials)
            analysisStruct.(roiMeanFCellNames{nn})(:, tt) = NaN(numFrames, 1);
            analysisStruct.(bkgrndMeanFCellNames{nn})(:, tt) = NaN(numFrames, 1);
            analysisStruct.(subScalarCellNames{nn})(tt, 1) = NaN;
            analysisStruct.(roiMeanFCorrectedCellNames{nn})(:, tt) = NaN(numFrames, 1);
            analysisStruct.(volpySubthreshFullCellNames{nn})(:, tt) = NaN(numFrames, 1);
        else
            analysisStruct.(roiMeanFCellNames{nn})(:, tt) = roiMeanF;
            analysisStruct.(bkgrndMeanFCellNames{nn})(:, tt) = bkgrndMeanF;
            analysisStruct.(subScalarCellNames{nn})(tt, 1) = alphaScalar;
            analysisStruct.(roiMeanFCorrectedCellNames{nn})(:, tt) = roiMeanFCorrected;
            analysisStruct.(volpySubthreshFullCellNames{nn})(:, tt) = volpySubFull;
        end
    end
end

%% Calculate mean response (and CI) for each hologram across trials and per condition
holoSortedMeanCellNames = cell(nCells, 1);
filtHoloSortedMeanCellNames = cell(nCells, 1);
holoSortedVolpySubMeanCellNames = cell(nCells, 1);
filtHoloSortedVolpySubMeanCellNames = cell(nCells, 1);
for nn = 1:nCells
    holoSortedMeanCellNames{nn}        = ['holoSortedImagingMean_', 'cell', num2str(nn)];
    filtHoloSortedMeanCellNames{nn}    = ['filtHoloSortedImagingMean_', 'cell', num2str(nn)];
    holoSortedVolpySubMeanCellNames{nn} = ['holoSortedVolpySubMean_', 'cell', num2str(nn)];
    filtHoloSortedVolpySubMeanCellNames{nn} = ['filtHoloSortedVolpySubMean_', 'cell', num2str(nn)];
    analysisStruct.(holoSortedMeanCellNames{nn})        = cell(nConds, 1);
    analysisStruct.(filtHoloSortedMeanCellNames{nn})    = cell(nConds, 1);
    analysisStruct.(holoSortedVolpySubMeanCellNames{nn}) = cell(nConds, 1);
    analysisStruct.(filtHoloSortedVolpySubMeanCellNames{nn}) = cell(nConds, 1);
    
    for cc = 1:nConds
        analysisStruct.(holoSortedMeanCellNames{nn}){cc}        = cell(nHolos(cc), 1);
        analysisStruct.(filtHoloSortedMeanCellNames{nn}){cc}    = cell(nHolos(cc), 1);
        analysisStruct.(holoSortedVolpySubMeanCellNames{nn}){cc} = cell(nHolos(cc), 1);
        analysisStruct.(filtHoloSortedVolpySubMeanCellNames{nn}){cc} = cell(nHolos(cc), 1);
    end
end

for nn = 1:nCells
    for cc = 1:nConds
        for hh = 1:nHolos(cc)
            analysisStruct.(holoSortedMeanCellNames{nn}){cc}{hh}        = nanmean(analysisStruct.(holoSortedImagingCellNames{nn}){cc}{hh}, 2);
            analysisStruct.(filtHoloSortedMeanCellNames{nn}){cc}{hh}    = nanmean(analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc}{hh}, 2);
            analysisStruct.(holoSortedVolpySubMeanCellNames{nn}){cc}{hh} = nanmean(analysisStruct.(holoSortedVolpySubCellNames{nn}){cc}{hh}, 2);
            analysisStruct.(filtHoloSortedVolpySubMeanCellNames{nn}){cc}{hh} = nanmean(analysisStruct.(filtHoloSortedVolpySubCellNames{nn}){cc}{hh}, 2);
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
voltMapping.holoSortedVolpySubCellNames                 = holoSortedVolpySubCellNames;
voltMapping.filtHoloSortedVolpySubCellNames             = filtHoloSortedVolpySubCellNames;
voltMapping.volpySubthreshFullCellNames                 = volpySubthreshFullCellNames;
voltMapping.holoSortedMeanCellNames                     = holoSortedMeanCellNames;
voltMapping.filtHoloSortedMeanCellNames                 = filtHoloSortedMeanCellNames;
voltMapping.holoSortedVolpySubMeanCellNames             = holoSortedVolpySubMeanCellNames;
voltMapping.filtHoloSortedVolpySubMeanCellNames         = filtHoloSortedVolpySubMeanCellNames;
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
    eval([structname '.' 'holoSortedVolpySubAllTrials' ' = analysisStruct.(holoSortedVolpySubCellNames{nn});']);
    eval([structname '.' 'filtHoloSortedVolpySubAllTrials' ' = analysisStruct.(filtHoloSortedVolpySubCellNames{nn});']);
    eval([structname '.' 'volpySubthreshFull' ' = analysisStruct.(volpySubthreshFullCellNames{nn});']);
    eval([structname '.' 'holoSortedImagingMean' ' = analysisStruct.(holoSortedMeanCellNames{nn});']);
    eval([structname '.' 'filtHoloSortedImagingMean' ' = analysisStruct.(filtHoloSortedMeanCellNames{nn});']);
    eval([structname '.' 'holoSortedVolpySubMean' ' = analysisStruct.(holoSortedVolpySubMeanCellNames{nn});']);
    eval([structname '.' 'filtHoloSortedVolpySubMean' ' = analysisStruct.(filtHoloSortedVolpySubMeanCellNames{nn});']);
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
            [exclFiltCIDffAllConds{cc}{hh}(:, 1)'*10, fliplr(exclFiltCIDffAllConds{cc}{hh}(:, 2)'*10)], [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
        % plot CI lowerbound
        plot(linspace(0, size(exclFiltCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclFiltCIDffAllConds{cc}{hh}, 1)), exclFiltCIDffAllConds{cc}{hh}(:, 1)*10, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);
        % plot CI upperbound
        plot(linspace(0, size(exclFiltCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclFiltCIDffAllConds{cc}{hh}, 1)), exclFiltCIDffAllConds{cc}{hh}(:, 2)*10, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);

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
        plot(linspace(0, size(exclFiltHoloSortedImagingMean{cc}{hh}, 1)/imagingFreq, size(exclFiltHoloSortedImagingMean{cc}{hh}, 1)), exclFiltHoloSortedImagingMean{cc}{hh}*10, '-', 'linewidth', 2, 'color', 'g');
        ylabel('dF/F (%)');
        xlabel('Time (s)');
        %         ylim([-1 max(holoSortedImagingMean{cc}{hh})*100])
        %         ylim([min(holoSortedImagingMean{cc}{hh})*100 max(holoSortedImagingMean{cc}{hh})*100])
        % ylim([-1 2])
        
        % if cc > 2
            ylim([-5 30])
        % end

        % ax.YColor = [0 1 0];
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
sanityTrialStep = 250;          % trials 1, 1+step, 1+2*step, ... up to nTrials
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
voltMapping.mouseID = ['voltMapping_Analysis_', expID, '_VolPyStyle'];
directory = '/Users/masatosadahiro/Documents/Data/Voltage Imaging/Voltage Imaging/voltMapping/Analysis Results/Analysis_newMCanddFFcalc';
fileName = [num2str(voltMapping.mouseID), '.mat'];

% Persist imaging traces/stats plus field-name cell arrays (needed to index analysisStruct after load).
varsToSave = [{'voltMapping', 'analysisStruct', ...
    'F0CellNames', 'roiMeanFCellNames', 'bkgrndMeanFCellNames', 'subScalarCellNames', ...
    'roiMeanFCorrectedCellNames', 'globalF0CellNames', 'dFCellNames', 'dFFCellNames', ...
    'holoSortedImagingCellNames', 'filtHoloSortedImagingCellNames', ...
    'holoSortedVolpySubCellNames', 'filtHoloSortedVolpySubCellNames', 'volpySubthreshFullCellNames', ...
    'holoSortedMeanCellNames', 'filtHoloSortedMeanCellNames', ...
    'holoSortedVolpySubMeanCellNames', 'filtHoloSortedVolpySubMeanCellNames', ...
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
