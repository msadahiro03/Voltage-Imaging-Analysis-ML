%% Voltage Imaging Mapping Analysis 070125
% This analysis is used for 2P optogenetic mapping experiments where the readout is either ephys or imaging.

%% Notes dump
% ipi, nHolos, among other variables need to be fixed to accommodate for
% 0mV trials, just because 0mV trials typically do not have stim profies
% because 0mV simply means no holograms are generated.

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
    % Check if the entry is a regular file (not a directory) and its name doesn't start with a period
    if ~ImgfolderContents(ii).isdir && ~startsWith(ImgfolderContents(ii).name, '.') && endsWith(ImgfolderContents(ii).name, fileType)
        % Add the file name to the cell array
        fileNames{ii, 1} = ImgfolderContents(ii).name;
    end
end
imagesIndex = find(~cellfun(@isempty, fileNames));

% Step 3: Setup struct
voltMapping = ExpStruct;

if exist('ExpStruct2', 'var') % If the experiment has a second patch electrode
    voltMapping.ExpStruct2 = ExpStruct2;
end

% Stimulation properties
imagingFreq = voltMapping.sampleFreq;
Fs = voltMapping.daqParams.Fs;
trialTime = voltMapping.daqParams.maxSweepLengthSec;
nTrials = length(voltMapping.trialCond); % ALTERNATIVELY "length(find(cellfun(@isempty, ExpStruct.inputs)==0))". Instead of "length(ExpStruct.inputs)" this puts out true number of trials successfully recorded
powers = voltMapping.outParams.power; % ALTERNATIVELY "unique(ExpStruct.trialCond)" what powers were used
nConds = length(voltMapping.outParams.sequence); % total number of powers used
nHolos = voltMapping.holoStimParams.nHolos; % number of holograms in grid
    nHolos(1) = max(nHolos); % Hack for 0 holos conditions because of 0mW trials involved

% Stimulation properties - assume all of these properties are the same across all conditions/powers
pulseDurs = unique(voltMapping.outParams.pulseDur);  % Assume same pulse duration or stim rate for all conditions
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

% % Step 5: Calculate frameclock offsets
% frameClock_inputs = [];
% risingEdgeIndices = [];
% risingEdgeTimes = [];
% if isfield(ExpStruct.dvStepParams, 'frameClock_inputs')
%     frameClock_inputs = ExpStruct.dvStepParams.frameClock_inputs;
%     for tt = 1:size(frameClock_inputs, 2)
%         trialTime = (0:length(frameClock_inputs(:, tt))-1)*(1/Fs);
%         riseThreshold = 0; % rising edge threshold, any number well above baseline will do
%         if any(frameClock_inputs(:, tt) == 1)
%             risingEdgeIndices(tt) = find(frameClock_inputs(:, tt) > riseThreshold, 1); % find the index for first rising edge
%             risingEdgeTimes(tt) = trialTime(risingEdgeIndices(tt));
%         else
%             risingEdgeIndices(tt) = 0;
%             risingEdgeTimes(tt) = 0;
%         end
%     end
% else
%     risingEdgeTimes = repmat(0.0011, 1, size(vsTest_inputs, 2) - size(excludeTrials, 2));
% end

% Step 6: Save directory
mouseID = ExpStruct.mouseID;
directory = 'D:\Data\Voltage Imaging\voltMapping\Analysis Results\';
fileName = ['voltMapping ', num2str(ExpStruct.mouseID), '.mat'];

voltMapping.imagesIndex   = imagesIndex;
voltMapping.imagingFreq   = imagingFreq;
voltMapping.UpOrDown      = UpOrDown;
voltMapping.ephysFilePath = ephysFilePath;
voltMapping.ImgsFilePath  = ImgsFilePath;

%% Trial exclusion and baselining
% Excludes all trials where the baseline Vm suddenly jumps (losing cell etc.)
vThreshold = -60; %set threshold Vm here
mappingInputs = ExpStruct.inputs;

% for tt = 1:length(ExpStruct.inputs)
%     ExpStruct.inputs{tt} = ExpStruct.inputs{tt}*2;
% end

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
% Stimulation window is the same as ones set for breaking apart ephys data

% Filtering parameters (if and when data gets filtered below)
cutOffFreq = 480;   % Cutoff frequency for Butterworth filter
[blp, alp] = butter(4, cutOffFreq/(Fs/2), 'low');  % 4th order Butterworth filter

holoSeqIndex = cell(nConds, 1);
holoSortedDataAllTrials = cell(nConds, 1);
for cc = 1:nConds
    holoSortedDataAllTrials{cc} = cell(nHolos(cc), 1);
end

% Set the stimulation window
postStimWindow = 50; % time(ms) to add to stim window after final pulse + ipi
preStimWindow = nextHoloDelay - postStimWindow; % time(ms) window before first pulse 

nPulseCoords = []; % indices of when the pulses happen across the stimulation window (in ephys samples)
for pp = 1:nPulses
    nPulseCoords = [nPulseCoords, ((pp-1)*ipi/1000*Fs) + preStimWindow/1000*Fs];
end

nPulseCoordsImaging = []; % indices of when the pulses happen across the stimulation window (in imaging frames)
for pp = 1:nPulses
    nPulseCoordsImaging = [nPulseCoordsImaging, ((pp-1)*ipi/1000*imagingFreq) + preStimWindow/1000*imagingFreq];
end

% Create dummy sequence for 0mV trials (if any)
if any(voltMapping.outParams.power == 0)
    zeroDummySequence = voltMapping.outParams.sequence{1, 2};
end

condSortedInputs = cell(nConds, 1);
for tt = 1:nTrials
    
    if isempty(voltMapping.outParams.sequenceThisTrial{tt}) % Hack for 0mV trials (replacing empty holo sequence with another from a random trial)
        voltMapping.outParams.sequenceThisTrial{tt} = zeroDummySequence;
    end

    holoSeqThisTrial = unique(voltMapping.outParams.sequenceThisTrial{tt}, 'stable') - min(unique(voltMapping.outParams.sequenceThisTrial{tt}, 'stable')) + 1;
    holoSeqIndex{voltMapping.trialCond(tt, 1)} = [holoSeqIndex{voltMapping.trialCond(tt, 1)}, holoSeqThisTrial']; % Compile hologram sequences for every trial

    sortingIndex = holoSeqIndex{voltMapping.trialCond(tt, 1)}(:, end) - min(holoSeqIndex{voltMapping.trialCond(tt, 1)})+1;
    % if ismember(tt, excludeTrials) 
    %     holoSeqIndex{voltMapping.trialCond(tt, 1)} = [holoSeqIndex{voltMapping.trialCond(tt, 1)}, NaN(nHolos(voltMapping.trialCond(tt, 1)), 1)]; % Find hologram sequences for every trial
    %     condSortedInputs{voltMapping.trialCond(tt, 1)} = [condSortedInputs{voltMapping.trialCond(tt, 1)}, NaN(length(mappingInputsBaselined{tt}), 1)];
    %     for hh = 1:nHolos(voltMapping.trialCond(tt, 1))
    %         holoSortedDataAllTrials{voltMapping.trialCond(tt, 1)}{sortingIndex(hh), 1} = [holoSortedDataAllTrials{voltMapping.trialCond(tt, 1)}{sortingIndex(hh), 1}, NaN((ipi*nPulses+(preStimWindow+postStimWindow))/1000*Fs+1, 1)];         
    %     end
    %     continue
    % end
    
    condSortedInputs{voltMapping.trialCond(tt, 1)} = [condSortedInputs{voltMapping.trialCond(tt, 1)}, mappingInputsBaselined{tt}];
    
%     sweepThisTrial =  mappingInputsBaselined{tt}; %Run this line instead of one below if data needs filtering
    sweepThisTrial = filtfilt(blp, alp, mappingInputsBaselined{tt}); %If ephys data needs to be filtered
    
    for hh = 1:nHolos(voltMapping.trialCond(tt, 1)) % NEED TO USE sorting index but...
        if isempty(voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}) % Hack for 0mV trials
            voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)} = voltMapping.outParams.firstStimTimes{1, 2};
        end

        thisHoloSweep = sweepThisTrial((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*Fs:((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*Fs+((ipi*nPulses+preStimWindow+postStimWindow)/1000)*Fs)); % Break whole sweep into stimulation-timed windows
        thisHoloSweep = thisHoloSweep - mean(thisHoloSweep(1:nPulseCoords(1)-100)); % Baseline the stimulation-timed window
   
        % if max(thisHoloSweep(nPulseCoords:end)) > 0.00% IMPORTANT: comment out this if condition (or lower the threshold to 0) if we want to use all trials instead of just those that actually show some depolarization
            holoSortedDataAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, end), 1} = [holoSortedDataAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, end), 1}, thisHoloSweep];
        % else
            % holoSortedDataAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, end), 1} = [holoSortedDataAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, end), 1}, NaN(length(thisHoloSweep), 1)];
        % end
    end
end

% % Show whole trial traces
% for tt = 1:20
%     figure(0+tt); clf;
%         plot(linspace(0, length(mappingInputsBaselined{tt})/Fs, length(mappingInputsBaselined{tt})), mappingInputsBaselined{tt});
% %     hold on
% %         plot(linspace(0, length(mappingInputsBaselined{tt})/Fs, length(mappingInputsBaselined{tt})), ExpStruct.outParams.nextHoloStims{1, 1} );
%     xlabel('time(s)')
%     ylabel('dV(mV)')
% end

% Show average traces for each hologram
holoSortedDataMean = cell(nConds, 1);
for cc = 1:nConds
    for hh = 1:nHolos(cc)
        holoSortedDataMean{cc}(:, hh) = nanmean(holoSortedDataAllTrials{cc}{hh}, 2);

        % Baseline the mean holo traces
%         holoSortedDataMean{cc}(:, hh) = holoSortedDataMean{cc}(:, hh) - mean(holoSortedDataMean{cc}((ipi+preStimWindow)/1000*Fs-100:(ipi+preStimWindow)/1000*Fs, hh));

    %     figure(cc*100+hh);
    %     clf
    %     for nn = 1:length(nPulseCoords)
    %         xline(nPulseCoords(nn)/Fs, '--', 'LineWidth', 1.5, 'color', [.8 .8 .8]);
    %     end
    %     hold on
    %     plot(linspace(0, length(holoSortedDataMean{cc}(:, hh))/Fs, length(holoSortedDataMean{cc}(:, hh))), holoSortedDataMean{cc}(:, hh), 'LineWidth', 1.5);
    % %     axis off
    %     hold off
    %     ylim([-0.2, 2]);
    %     xlim([0, 0.15]);
    %     xlabel('Time(s)')
    %     ylabel('mV')
        % pause
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

% % Combine all hologram traces, calculate one grand mean
% holoCombinedDataAllTrials = [];
% for cc = 1
%     for hh = 1:nHolos(1)
%         holoCombinedDataAllTrials = [holoCombinedDataAllTrials, holoSortedDataAllTrials{cc}{hh}(:, :)];
%     end
% end
% holoComboDataGrandMean = nanmean(holoCombinedDataAllTrials, 2);

voltMapping.ephys.holoSeqIndex            = holoSeqIndex; 
voltMapping.ephys.holoSortedDataAllTrials = holoSortedDataAllTrials;
voltMapping.ephys.nPulseCoords            = nPulseCoords;
voltMapping.ephys.CIephysAllConds         = CIephysAllConds;
voltMapping.ephys.holoSortedDataMean      = holoSortedDataMean;
voltMapping.ephys.preStimWindow           = preStimWindow;
voltMapping.ephys.postStimWindow          = postStimWindow;
voltMapping.ephys.condSortedInputs        = condSortedInputs;

%% Calculate ROI mask
input('Averaging a set of trials for viewing and selecting ROIs (continue or ctrl+c to stop!)');
% Step 1: Calculate mask for ROI 
if length(imagesIndex) < 100
    randTrialsForMask = randperm(length(imagesIndex), length(imagesIndex));
else
    randTrialsForMask = randperm(length(imagesIndex), length(imagesIndex)); % Select random n trials/sweeps to be used to generate an average z-stacked image;
end

% Preallocate maxDvStack, which takes average z-stack image for each trial used
currImgPath = [ImgfolderContents(imagesIndex(1)).folder, '/', ImgfolderContents(imagesIndex(1)).name];
info = imfinfo(currImgPath);
maxDvStack = zeros(info(1).Height, info(1).Width, length(randTrialsForMask));

counter = 0;
for tt = randTrialsForMask
    counter = counter+1;
    disp(['Trial number: ', num2str(counter)]);

    if ismember(imagesIndex(tt), excludeTrials)
        continue
    end

    % %%%%%%%%%%%%%%%% Old Loading Method
    % tic
    % currImgPath = [ImgfolderContents(imagesIndex(tt)).folder, '/', ImgfolderContents(imagesIndex(tt)).name];
    % info = imfinfo(currImgPath);
    % numFrames = numel(info);
    % if numFrames > imagingFreq*4 % Just take frames from first 4 seconds (if actual sweeps are longer than 4s), because not entire sweeps are needed
    %     numFrames = floor(imagingFreq*4);
    % end
    % 
    % % Preallocate the image stack
    % imageStack = zeros(info(1).Height, info(1).Width, numFrames);
    % 
    % % Read each frame and store in the stack
    % for frameIndex = 1:numFrames
    %     imageStack(: ,:, frameIndex) = imread(currImgPath, 'Index', frameIndex, 'Info', info);
    % end
    % toc
    % %%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%% New Loading Method
    tic
    currImgPath = fullfile(ImgfolderContents(imagesIndex(tt)).folder, ...
                           ImgfolderContents(imagesIndex(tt)).name);
    
    t = Tiff(currImgPath, 'r');
    
    % Get frame count without imfinfo overhead (often faster):
    % Note: some TIFFs require directory walking anyway.
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
    toc
    %%%%%%%%%%%%%%%%

    imageStackMean = zeros(info(1).Height, info(1).Width);
    imageStackMean = mean(imageStack(:, :, :), 3);
    %     for ff = 1:size(imageStack, 3); figure(101010); imagesc(imageStack(:, :, ff)); caxis([min(min(min(imageStack(:,:,:)))), max(max(max(imageStack(:,:,:))))]); axis equal; axis image; pause; end
    maxDvStack(:, :, tt) = imageStackMean;
end

meanMaxDvStack = mean(maxDvStack(:, :, :), 3); % The average of all the per trial average z-stacks
meanFluorMaxDvStack = meanMaxDvStack; % Grand average z-stack image to use for hand selecting ROI
% maxDvStackBrighter = meanFluorMaxDvStack + 255; % Brighter z-stack image if needed
% maxDvStackBrighter(maxDvStackBrighter > 255) = 255; % Clamping for brightness effect

figure(9); set(gcf, 'Position',  [100, 100, 1800, 900]); clf; 
colormap('winter'); imagesc(meanFluorMaxDvStack); axis equal; axis image; colorbar; set(gca, 'fontsize', 12); % caxis([-7 -4]);
nCells = input('How many neurons to analyze?: ', 's');
nCells = str2double(nCells);

% Hand select cell or area of interest, by freehand drawing
roughRoiXAllCells = cell(nCells, 1);
roughRoiYAllCells = cell(nCells, 1);
roiXAllCells = cell(nCells, 1);
roiYAllCells = cell(nCells, 1);
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

    % roiXAllCells{nn} = roughRoiX;
    % roiYAllCells{nn} = roughRoiY;
    close(f1);
end

for nn = 1:nCells
    figure(10+(nn));
    set(gcf, 'Position',  [100, 100, 1800, 600]);
    clf
    subplot(3,1,1)
    colormap('winter'); imagesc(meanFluorMaxDvStack); axis equal; axis image; set(gca, 'fontsize', 12);

    % Calculate mean fluorescence of the area of interest and then form into final ROI mask
    roiMeanMaxDvStack = zeros(size(maxDvStack, 1), size(maxDvStack, 2));
    for rr = 1:length(roughRoiXAllCells{nn})
        roiMeanMaxDvStack(roughRoiXAllCells{nn}(rr), roughRoiYAllCells{nn}(rr)) = mean(maxDvStack(roughRoiXAllCells{nn}(rr), roughRoiYAllCells{nn}(rr),:), 3);
    end
    
    subplot(3,1,2);
    imagesc(roiMeanMaxDvStack); axis equal; axis image; set(gca, 'fontsize', 12);
    
    %%%%%%%%%%%%%%%%%%
    % Old ROI selection method
    % stdFluor = std(nonzeros(roiMeanMaxDvStack));
    % meanFluor = mean(nonzeros(roiMeanMaxDvStack));
    % 
    % % Designate cutoff fluorescence for pixels to be selected for ROI
    % cutOffFluor = meanFluor; %stdFluor*1 + meanFluor; % currently cutoff is 1 standard devs from mean fluorescence
    % 
    % roiStack = roiMeanMaxDvStack;
    % roiStack(roiStack <= cutOffFluor) = 0;
    % roiStack(roiStack > 0) = 1;
    % 
    % subplot(3,1,3);
    % imagesc(roiStack); axis equal; axis image; set(gca, 'fontsize', 12);
    % 
    % [roiX, roiY] = find(roiStack);
    % roiXAllCells{nn} = roiX;
    % roiYAllCells{nn} = roiY;

    %%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%
    % New ROI selection method
    roiMaxStackDouble   = im2double(roiMeanMaxDvStack);
    % 1) ROI mask from your "selected pixels"
    roiPixels = roiMaxStackDouble > 0;
    
    % 2) Denoise / smooth ONLY for feature extraction
    % (Gaussian is usually fine; you can try medfilt2 if salt/pepper)
    % Is = imgaussfilt(I, 1.5);   % sigma in pixels (try 0.6 to 1.5)
    roiMaxStackFilt = imgaussfilt(roiMaxStackDouble, 0.7); 
    
    % Optional: normalize intensity inside ROI to reduce bias
    vals = roiMaxStackFilt(roiPixels);
    if ~isempty(vals)
        lo = prctile(vals, 10);
        hi = prctile(vals, 99);
        roiMaxStackNorm = (roiMaxStackFilt - lo) / max(hi - lo, eps);
        roiMaxStackNorm = min(max(roiMaxStackNorm, 0), 1);
    else
        roiMaxStackNorm = roiMaxStackFilt;
    end
    
    % 3) Ridge enhancement (membrane as bright, thin fiber-like structure)
    % fibermetric expects grayscale image. It returns higher values for line-like features.
    % StructureSensitivity ~ expected half-width (in px) of the ridge.
    roiMaxStackRidge = fibermetric(roiMaxStackNorm, 'StructureSensitivity', 2);  % try 1 to 4
    stdFluor = std(nonzeros(roiMaxStackRidge));
    meanFluor = mean(nonzeros(roiMaxStackRidge));
    
    roiMaxStackRidgeReduced = roiMaxStackRidge; % USE this -RidgeReduced (cleaner version of roiMaxStackRdige) instead if needed (change to roiMaxStackRidgeReduced to select background ROI)
    roiMaxStackRidgeReduced(roiMaxStackRidgeReduced <= stdFluor*1-meanFluor) = 0;
    roiMaxStackRidge(roiMaxStackRidge > 0) = 1;
    roiMaxStackRidgeReduced(roiMaxStackRidgeReduced > 0) = 1;

    subplot(3,1,3);
    imagesc(roiMaxStackRidge); axis equal; axis image; set(gca, 'fontsize', 12);
    hold on

    [roiX, roiY] = find(roiMaxStackRidge);
    roiXAllCells{nn} = roiX;
    roiYAllCells{nn} = roiY;

    % [roiX, roiY] = find(roiMaxStackRidgeReduced);
    % roiXAllCells{nn} = roiX;
    % roiYAllCells{nn} = roiY;

    %%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%
    % Background ROI
    % Selection of background ROI for this cell
    % Parameters (pixels)
    innerBuffer = 2;   % exclude this many px around membrane (spillover guard)
    ringWidth   = 2;   % thickness of neuropil ring
    minArea     = 50;  % ensure enough pixels
    
    % 1) Make an annulus around membrane
    innerSelect = imdilate(roiMaxStackRidge, strel('disk', innerBuffer));
    outerSelect = imdilate(roiMaxStackRidge, strel('disk', innerBuffer + ringWidth));
    backgroundRing  = outerSelect & ~innerSelect;
    
    % 2) Exclude other bright structures (automatic)
    % Mask out the top X% brightest pixels in the ring to avoid other cells
    vals = roiMaxStackDouble(backgroundRing);
    if ~isempty(vals)
        brightCut = prctile(vals, 95);     % tune 90–99
        ringClean = backgroundRing & (roiMaxStackDouble <= brightCut);
    else
        ringClean = backgroundRing;
    end
    
    % 3) Remove tiny islands and ensure contiguity around the cell
    ringClean = bwareaopen(ringClean, 7);
    % If ring still too small, relax exclusion
    if nnz(ringClean) < minArea
        ringClean = backgroundRing; % fall back to raw ring
    end

    [bkgrndRoiX, bkgrndRoiY] = find(ringClean);    
    bkgrndRoiXAllCells{nn} = bkgrndRoiX;
    bkgrndRoiYAllCells{nn} = bkgrndRoiY;

    [y,x] = find(ringClean);
    plot(x, y,'r.','MarkerSize',6);
    hold off
    %%%%%%%%%%%%%%%%%%
end

allRois = zeros(size(meanFluorMaxDvStack, 1), size(meanFluorMaxDvStack, 2));
centerXY = [];
for nn = 1:nCells
    thisRoi = zeros(size(meanFluorMaxDvStack, 1), size(meanFluorMaxDvStack, 2)); 
    for rr = 1:length(roiXAllCells{nn})
        thisRoi(roiXAllCells{nn}(rr), roiYAllCells{nn}(rr)) = mean(meanFluorMaxDvStack(roiXAllCells{nn}(rr), roiYAllCells{nn}(rr),:), 3);
    end
    allRois = allRois + thisRoi;
    centerXY(nn, 1) = (min(roiXAllCells{nn}) + max(roiXAllCells{nn}))/2;
    centerXY(nn, 2) = (min(roiYAllCells{nn}) + max(roiYAllCells{nn}))/2;
end
figure(20); set(gcf, 'Position',  [100, 100, 1800, 800]); clf;
subplot(2, 1, 1);
colormap('winter'); imagesc(meanFluorMaxDvStack); axis equal; axis image; colorbar; set(gca, 'fontsize', 12); % caxis([-7 -4]);
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
figure(21); set(gcf, 'Position',  [100, 100, 1800, 300]);; clf; 
colormap('winter'); axis equal; axis image; colorbar; set(gca, 'fontsize', 12);
imagesc(allRois); axis equal; axis image; colorbar;
hold on;
[y,x] = find(allBkgrndRois);
plot(x, y,'r.','MarkerSize',6);
hold on;

voltMapping.nCells              = nCells;
voltMapping.maxDvStack          = maxDvStack;
voltMapping.meanFluorMaxDvStack = meanFluorMaxDvStack;
% voltMapping.maxDvStackBrighter  = maxDvStackBrighter;
voltMapping.roiMeanMaxDvStack   = roiMeanMaxDvStack;
% voltMapping.cutOffFluor         = cutOffFluor;
voltMapping.roiXAllCells        = roughRoiXAllCells;
voltMapping.roiYAllCells        = roughRoiYAllCells;
voltMapping.roiXAllCells        = roiXAllCells;
voltMapping.roiYAllCells        = roiYAllCells;
voltMapping.bkgrndRoiXAllCells  = bkgrndRoiXAllCells;
voltMapping.bkgrndRoiYAllCells  = bkgrndRoiYAllCells;
voltMapping.allRois             = allRois;

%% Calculate df and mean df for all holos sorted by trial conditions
% Use these lines if I need to further filter the df or dff trace
input('This step splits whole trace by holograms stimmed and calculates df/f in each stim window (ctrl+c to stop!)');
cutOffFreq = 50;   % Cutoff frequency
[b, a] = butter(4, cutOffFreq/(imagingFreq/2));  % 4th order Butterworth filter

% Preallocate for data for every cell, sorted by holograms, all trials, across conditions
holoSortedImagingCellNames = cell(nCells, 1);
filtHoloSortedImagingCellNames = cell(nCells, 1);
for nn = 1:nCells
    holoSortedImagingCellNames{nn} = ['holoSortedImagingAllTrials_', 'cell', num2str(nn)];
    filtHoloSortedImagingCellNames{nn} = ['filtHoloSortedImagingAllTrials_', 'cell', num2str(nn)];
    analysisStruct.(holoSortedImagingCellNames{nn}) = cell(nConds, 1);
    analysisStruct.(filtHoloSortedImagingCellNames{nn}) = cell(nConds, 1);
    
    for cc = 1:nConds
        analysisStruct.(holoSortedImagingCellNames{nn}){cc} = cell(nHolos(cc), 1);
        analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc} = cell(nHolos(cc), 1);
    end
end

counter = 0;
for tt = 1:nTrials %size(vsTest_inputs, 2)
    windowLimits = []; firstLimit = []; windowTime = []; segmentTime = []; numSegments = []; windowWidth = [];
    counter = counter+1;
    disp(['Trial number: ', num2str(counter)]);
    
    if ismember(tt, excludeTrials)
        continue
    end
    
    % Read the multi-frame image
    currImgPath = [ImgfolderContents(imagesIndex(tt)).folder, '/', ImgfolderContents(imagesIndex(tt)).name];
    info = imfinfo(currImgPath);
    numFrames = numel(info);

    % Preallocate the image stack
    imageStack = zeros(info(1).Height, info(1).Width, numFrames);
    
    % Read each frame and store in the imageStack
    for frameIndex = 1:numFrames
        imageStack(:,:,frameIndex) = imread(currImgPath, 'Index', frameIndex, 'Info', info);

        % % Run these lines if artifact removal is necessary, comment out if not needed
        % [cleanFrame] = VoltImg_mapping_removeArtifact(imageStack(:,:,frameIndex));
        % imageStack(:,:,frameIndex) = cleanFrame; %replace the raw frame with cleaned up frame where artifact-corrupt lines are NaN'd
    end
    
    counter2 = 0;
    for nn = 1:nCells
        counter2 = counter2+1;
        disp(['Cell: ', num2str(counter2)]);
        startTimeImaging = floor(startTime*imagingFreq); % Start time of first hologram stimulation in image sample
        
        % % Baselining the membrane ROI to prestimulus at beginning of sweep: basic parameters
        % baselinePreImageStack = imageStack(roiXAllCells{nn}, roiYAllCells{nn}, 1:startTimeImaging); % image stack in ROI before start of stimulation
        % roiBaselinePreImageStack = [];
        % for ff = 1:size(baselinePreImageStack, 3)
        %     roiBaselinePreImageStack(ff, 1) = mean(mean(baselinePreImageStack(:, :, ff))); % Generate mean trace of ROI fluorescence during period before stimulation
        % end
        % 
        % % Baselining to prestimulus at beginning of sweep: parameters for moving window variance during pre-stim period
        % windowLimits = [1/imagingFreq, size(roiBaselinePreImageStack, 1)/imagingFreq]; % in sec The range within each sweep where variance is sampled by 10 sample windows, in ms. Will leave it as 500 ms total width, starting from 110ms after start (after the step pulse).
        % firstLimit = windowLimits(1)*imagingFreq; % Start of the first sampling window (in samples).
        % windowTime = ceil(size(roiBaselinePreImageStack, 1)/imagingFreq*1000)/1000; % Time of entire window (in seconds)
        % segmentTime = 0.020; % time length of each window (s)
        % numSegments = windowTime/segmentTime; % number of sample windows within the limit
        % windowWidth = floor((windowLimits(2) - windowLimits(1))/numSegments*imagingFreq); % Span of each sample window (in samples).
        % 
        % % Calculate variances/fanofactors in the prestimulus period using rolling window
        % varBaselinePre = movvar(roiBaselinePreImageStack, windowWidth); % moving window variance calculated across prestimulus 
        % fanoBaselinePre = varBaselinePre/(mean(roiBaselinePreImageStack));
        % 
        % % % Whole sweep baseline Method1: Baseline entire trace to period before first stimulation by thresholding the fanofactor
        % % % Fluorescence traces of pre-stim and post-stim periods
        % % % Choose variance/fanofactor threshold
        % % varThresholdPre = mean(varBaselinePre) - std(varBaselinePre); % variance threshold (mean-std) 
        % % fanoThresholdPre = mean(fanoBaselinePre) - std(fanoBaselinePre); % fanofactor (variance/mean)
        % % % Find variance/fanofactor beneath threshold
        % % [vLowestPre, ~] = find(varBaseline < varThreshold);
        % % [fanoLowestPre, ~] = find(fanoBaseline < fanoThreshold);
        % % % Establish the baseline fluorescence based on all the points corresponding to the lowest moments of fano
        % % roiBaselineMean = mean(roiBaselinePreImageStack(fanoLowestPre));
        % 
        % % Whole sweep baseline Method2: Calculate baseline value before first stimulation based on quantile bottom 10% of fanofactor 
        % % Step 1: Define quantile threshold (e.g., bottom 10%)
        % q = 0.10;  % Change to 0.05 for bottom 5%, etc.
        % quantileCutoff = quantile(fanoBaselinePre, q);
        % % Step 2: Select low points
        % [fanoLowestPre, ~] = find(fanoBaselinePre < quantileCutoff);
        % roiBaselineThisCellTrial = mean(roiBaselinePreImageStack(fanoLowestPre)); % baseline value to subtract
        % 
        % varBaselinePre = [];
        % fanoBaselinePre = [];
        % % Baselining the background ROI to prestimulus at beginning of sweep: basic parameters
        % bkgrndBaselinePreImageStack = imageStack(bkgrndRoiXAllCells{nn}, bkgrndRoiYAllCells{nn}, 1:startTimeImaging); % image stack in background before start of stimulation
        % bkgrndBaselinePreImageStack = [];
        % for ff = 1:size(bkgrndBaselinePreImageStack, 3)
        %     bkgrndBaselinePreImageStack(ff, 1) = mean(mean(bkgrndBaselinePreImageStack(:, :, ff))); % Generate mean trace of ROI fluorescence during period before stimulation
        % end
        % 
        % % Baselining to prestimulus at beginning of sweep: parameters for moving window variance during pre-stim period
        % windowLimits = [1/imagingFreq, size(bkgrndBaselinePreImageStack, 1)/imagingFreq]; % in sec The range within each sweep where variance is sampled by 10 sample windows, in ms. Will leave it as 500 ms total width, starting from 110ms after start (after the step pulse).
        % firstLimit = windowLimits(1)*imagingFreq; % Start of the first sampling window (in samples).
        % windowTime = ceil(size(bkgrndBaselinePreImageStack, 1)/imagingFreq*1000)/1000; % Time of entire window (in seconds)
        % segmentTime = 0.020; % time length of each window (s)
        % numSegments = windowTime/segmentTime; % number of sample windows within the limit
        % windowWidth = floor((windowLimits(2) - windowLimits(1))/numSegments*imagingFreq); % Span of each sample window (in samples).
        % 
        % % Calculate variances/fanofactors in the prestimulus period using rolling window
        % varBaselinePre = movvar(bkgrndBaselinePreImageStack, windowWidth); % moving window variance calculated across prestimulus 
        % fanoBaselinePre = varBaselinePre/(mean(bkgrndBaselinePreImageStack));
        % 
        % % % Whole sweep baseline Method1: Baseline entire trace to period before first stimulation by thresholding the fanofactor
        % % % Fluorescence traces of pre-stim and post-stim periods
        % % % Choose variance/fanofactor threshold
        % % varThresholdPre = mean(varBaselinePre) - std(varBaselinePre); % variance threshold (mean-std) 
        % % fanoThresholdPre = mean(fanoBaselinePre) - std(fanoBaselinePre); % fanofactor (variance/mean)
        % % % Find variance/fanofactor beneath threshold
        % % [vLowestPre, ~] = find(varBaseline < varThreshold);
        % % [fanoLowestPre, ~] = find(fanoBaseline < fanoThreshold);
        % % % Establish the baseline fluorescence based on all the points corresponding to the lowest moments of fano
        % % roiBaselineMean = mean(roiBaselinePreImageStack(fanoLowestPre));
        % 
        % % Whole sweep baseline Method2: Calculate baseline value before first stimulation based on quantile bottom 10% of fanofactor 
        % % Step 1: Define quantile threshold (e.g., bottom 10%)
        % q = 0.10;  % Change to 0.05 for bottom 5%, etc.
        % quantileCutoff = quantile(fanoBaselinePre, q);
        % % Step 2: Select low points
        % [fanoLowestPre, ~] = find(fanoBaselinePre < quantileCutoff);
        % bkgrndBaselineThisCellTrial = mean(bkgrndBaselinePreImageStack(fanoLowestPre)); % baseline value to subtract


    
        % Break apart the imageStack for this trial into stim windows and rearrange according to hologram sequence 
        if isempty((voltMapping.outParams.sequenceThisTrial{tt})) % Hack for 0mV trials (replacing empty holo sequence with another from a random trial)
            voltMapping.outParams.sequenceThisTrial{tt} = zeroDummySequence;
        end

        holoSeqThisTrial = (unique(voltMapping.outParams.sequenceThisTrial{tt}, 'stable') - min(unique(voltMapping.outParams.sequenceThisTrial{tt}, 'stable')) + 1)'; % Hologram sequence for this trial    
        for hh = 1:nHolos(voltMapping.trialCond(tt, 1))
            windowLimits = []; firstLimit = []; windowTime = []; segmentTime = []; numSegments = []; windowWidth = [];
    
            % Extract the frames associated with the current hologram, including pre and post stimulation windows
            if isempty(voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}) % Hack for 0mV trials
                voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)} = voltMapping.outParams.firstStimTimes{1, 2};
            end
            
            % if ceil((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*imagingFreq)+ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq) > numFrames
            %     continue
            % end
            
            roiFramesThisHolo = imageStack(roiXAllCells{nn}, roiYAllCells{nn}, ...
                floor((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*imagingFreq):(ceil((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*imagingFreq)+ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)));
            bkgrndFramesThisHolo = imageStack(bkgrndRoiXAllCells{nn}, bkgrndRoiYAllCells{nn}, ...
                floor((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*imagingFreq):(ceil((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*imagingFreq)+ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)));
            
            % Stim window baselining method: Same as whole sweep baselining but baseline value is set to moment of least variability right before specific target/hologram
            % Image stack in ROI before stim of this target
            roiFramesThisHoloPreStim = roiFramesThisHolo(:, :, 1:preStimWindow/1000*imagingFreq);
            bkgrndFramesThisHoloPreStim = bkgrndFramesThisHolo(:, :, 1:preStimWindow/1000*imagingFreq);

            roiTraceThisHoloPreStim = [];
            for ff = 1:size(roiFramesThisHoloPreStim, 3)
                roiTraceThisHoloPreStim(ff, 1) = mean(mean(roiFramesThisHoloPreStim(:, :, ff))); % Generate mean trace of ROI fluorescence during period before stimulation
            end

            bkgrndTraceThisHoloPreStim = [];
            for ff = 1:size(bkgrndFramesThisHoloPreStim, 3)
                bkgrndTraceThisHoloPreStim(ff, 1) = mean(mean(bkgrndFramesThisHoloPreStim(:, :, ff))); % Generate mean trace of ROI fluorescence during period before stimulation
            end
            
            % Baselining to prestimulus at beginning of sweep: parameters for moving window variance during pre-stim period
            % windowLimits = [1/imagingFreq, size(roiTraceThisHoloPreStim, 1)/imagingFreq]; % in sec The range within each sweep where variance is sampled by 10 sample windows, in ms. Will leave it as 500 ms total width, starting from 110ms after start (after the step pulse).
            % firstLimit = windowLimits(1)*imagingFreq; % Start of the first sampling window (in samples).
            % windowTime = ceil(size(roiTraceThisHoloPreStim, 1)/imagingFreq*1000)/1000;
            % segmentTime = 0.02; % time length of each window (s)
            % numSegments = windowTime/segmentTime; % number of sample windows within the limit
            % windowWidth = floor((windowLimits(2) - windowLimits(1))/numSegments*imagingFreq); % Span of each sample window (in samples).
            % 
            % % Calculate variances/fanofactors in the prestimulus period using rolling window
            % varBaselineThisHolo = movvar(roiTraceThisHoloPreStim, windowWidth); % moving window variance calculated across prestimulus 
            % fanoBaselineThisHolo = varBaselineThisHolo/(mean(roiTraceThisHoloPreStim));
            % 
            % % Calculate baseline value from window before this hologram stimulation based on quantile bottom 10% of fanofactor 
            % % Step 1: Define quantile threshold (e.g., bottom 10%)
            % q = 0.10;  % Change to 0.05 for bottom 5%, etc.
            % quantileCutoff = quantile(fanoBaselineThisHolo(2:end), q);
            % % Step 2: Select low points
            % [fanoLowestThisHolo, ~] = find(fanoBaselineThisHolo(2:end) < quantileCutoff);
            % roiBaselineThisHolo = mean(roiTraceThisHoloPreStim(fanoLowestThisHolo)); % baseline value to subtract
            roiBaselineThisHolo = mean(roiTraceThisHoloPreStim); % baseline value to subtract
            bkgrndBaselineThisHolo = mean(bkgrndTraceThisHoloPreStim);
    
            dfThisHolo = [];
            dffThisHolo = [];
            for ff = 1:size(roiFramesThisHolo, 3)
                % ROI pixels for this frame
                currFrameRoi = roiFramesThisHolo(:, :, ff);
                currFrameRoiMean = nanmean(currFrameRoi, 'all'); % The mean across all pixels in the ROI for the select frame
                % bkgrnd pixels for this frame
                currFrameBkgrnd = roiFramesThisHolo(:, :, ff);
                currFrameRoiMean = nanmean(currFrameRoi, 'all'); % The mean across all pixels in the ROI for the select frame               

                % Calculate df
                intensityChange = currFrameRoiMean - roiBaselineThisHolo; % essentially, intensityChange = df, and dff is intensityChnage/roiBaselineMean
                dfThisHolo = [dfThisHolo; intensityChange];
                dffThisHolo = [dffThisHolo; intensityChange/roiBaselineThisHolo];
            end
            
            if UpOrDown == '2'
                dfThisHolo = -dfThisHolo;
                dffThisHolo = -dffThisHolo;
            elseif UpOrDown =='1'
                dfThisHolo = dfThisHolo;
                dffThisHolo = dffThisHolo;
            end
            
            filtdffThisHolo = filter(b, a, dffThisHolo); % Run this if the voltage imaging trace need filtering
            
            % Re-baseline the dff data
    %         dffThisHolo = dffThisHolo - mean(dffThisHolo(1:ceil(preStimWindow/1000*imagingFreq)));
            if ismember(tt, excludeTrials) % isnan(holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, tt))
                %           holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{hh, 1}(:, tt) = NaN(ceil(ipi*nPulses/1000*imagingFreq)+1, 1);
                analysisStruct.(holoSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(holoSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, NaN(ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)+1, 1)];
                analysisStruct.(filtHoloSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(filtHoloSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, NaN(ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)+1, 1)];
            else
                %           holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, tt), 1}(:, tt) = dfThisHolo;
                analysisStruct.(holoSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(holoSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, dffThisHolo];
                analysisStruct.(filtHoloSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(filtHoloSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, filtdffThisHolo];
            end        
        end % nHolos
    end % nCells
end % nTrials

%% Calculate mean response (and CI) for each hologram across trials and per condition
holoSortedMeanCellNames = cell(nCells, 1);
filtHoloSortedMeanCellNames = cell(nCells, 1);
for nn = 1:nCells
    holoSortedMeanCellNames{nn} = ['holoSortedImagingMean_', 'cell', num2str(nn)];
    filtHoloSortedMeanCellNames{nn} = ['filtHoloSortedImagingMean_', 'cell', num2str(nn)];
    analysisStruct.(holoSortedMeanCellNames{nn}) = cell(nConds, 1);
    analysisStruct.(filtHoloSortedMeanCellNames{nn}) = cell(nConds, 1);
    
    for cc = 1:nConds
        analysisStruct.(holoSortedMeanCellNames{nn}){cc} = cell(nHolos(cc), 1);
        analysisStruct.(filtHoloSortedMeanCellNames{nn}){cc} = cell(nHolos(cc), 1);
    end
end

for nn = 1:nCells
    for cc = 1:nConds
        for hh = 1:nHolos(cc)
            analysisStruct.(holoSortedMeanCellNames{nn}){cc}{hh} = nanmean(analysisStruct.(holoSortedImagingCellNames{nn}){cc}{hh}, 2);
            analysisStruct.(filtHoloSortedMeanCellNames{nn}){cc}{hh} = nanmean(analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc}{hh}, 2);        
        end
    end
end

% % Combine all hologram traces, calculate one grand mean. Not a useful step - wrote this for the hell of it.
% holoComboImagingAllTrials = [];
% for cc = 1
%     for hh = 1:nHolos(1)
%         holoComboImagingAllTrials = [holoComboImagingAllTrials, holoSortedImagingAllTrials_ALT{cc}{hh}(:, :)];
%     end
% end
% holoComboImagingGrandMean = nanmean(holoComboImagingAllTrials, 2);

CIDffAllCondsCellNames = cell(nCells, 1);
filtCIDffAllCondsCellNames = cell(nCells, 1);
for nn = 1:nCells
    CIDffAllCondsCellNames{nn} = ['CIDffAllConds_', 'cell', num2str(nn)];
    filtCIDffAllCondsCellNames{nn} = ['filtCIDffAllConds_', 'cell', num2str(nn)];
    analysisStruct.(CIDffAllCondsCellNames{nn}) = cell(nConds, 1);
    analysisStruct.(filtCIDffAllCondsCellNames{nn}) = cell(nConds, 1);
end

for nn = 1:nCells
    for cc = 1:nConds 
        for hh = 1:nHolos(cc)
            confidence_level = 0.95;
            means = nanmean(analysisStruct.(holoSortedImagingCellNames{nn}){cc}{hh, 1}, 2);
            filtMeans = nanmean(analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc}{hh, 1}, 2);       
            std_errors = std(analysisStruct.(holoSortedImagingCellNames{nn}){cc}{hh, 1}, 0, 2, "omitnan") / sqrt(size(analysisStruct.(holoSortedImagingCellNames{nn}){cc}{hh, 1}, 2));
            filtStd_errors = std(analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc}{hh, 1}, 0, 2, "omitnan") / sqrt(size(analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc}{hh, 1}, 2));
      
            t_score = tinv((1 + confidence_level) / 2, size(analysisStruct.(holoSortedImagingCellNames{nn}){cc}{hh, 1}, 2) - 1);
            filtT_score = tinv((1 + confidence_level) / 2, size(analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc}{hh, 1}, 2) - 1);
            margin_of_error = t_score * std_errors;
            filtMargin_of_error = filtT_score * filtStd_errors;
            lower_bounds = means - margin_of_error;
            filtLower_bounds = filtMeans - filtMargin_of_error;
            upper_bounds = means + margin_of_error;
            filtUpper_bounds = filtMeans + filtMargin_of_error;
            if UpOrDown == '2'
                analysisStruct.(CIDffAllCondsCellNames{nn}){cc}{hh, 1} = [lower_bounds, upper_bounds];
                analysisStruct.(filtCIDffAllCondsCellNames{nn}){cc}{hh, 1} = [filtLower_bounds, filtUpper_bounds];
            elseif UpOrDown =='1'
                analysisStruct.(CIDffAllCondsCellNames{nn}){cc}{hh, 1} = [-lower_bounds, -upper_bounds];
                analysisStruct.(filtCIDffAllCondsCellNames{nn}){cc}{hh, 1} = [-filtLower_bounds, -filtUpper_bounds];
            end
        end
    end
end

voltMapping.nPulseCoordsImaging = nPulseCoordsImaging;
voltMapping.holoSeqIndex        = holoSeqIndex;
voltMapping.preStimWindow       = preStimWindow;
voltMapping.postStimWindow      = postStimWindow;

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

%% Align imaging with ephys traces
nn = double(input('which cell number? '));

CIDffAllConds = voltMapping.(cellID{nn}).CIDffAllConds;
holoSortedImagingMean = voltMapping.(cellID{nn}).holoSortedImagingMean;
filtCIDffAllConds = voltMapping.(cellID{nn}).filtCIDffAllConds;
filtHoloSortedImagingMean = voltMapping.(cellID{nn}).filtHoloSortedImagingMean;

for cc = 1:nConds
    for hh = 1:nHolos(cc)
        figure(cc*1000+hh);
        % set(gcf, 'Position',  [100, 100, 1600, 800])
        clf
        subplot(1,2,1)
        hold on
        fill([linspace(0, size(CIDffAllConds{cc}{hh}, 1)/imagingFreq, size(CIDffAllConds{cc}{hh}, 1)), fliplr(linspace(0, size(CIDffAllConds{cc}{hh}, 1)/imagingFreq, size(CIDffAllConds{cc}{hh}, 1)))],...
            [CIDffAllConds{cc}{hh}(:, 1)'*100, fliplr(CIDffAllConds{cc}{hh}(:, 2)'*100)], [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
        % plot CI lowerbound
        plot(linspace(0, size(CIDffAllConds{cc}{hh}, 1)/imagingFreq, size(CIDffAllConds{cc}{hh}, 1)), CIDffAllConds{cc}{hh}(:, 1)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);
        % plot CI upperbound
        plot(linspace(0, size(CIDffAllConds{cc}{hh}, 1)/imagingFreq, size(CIDffAllConds{cc}{hh}, 1)), CIDffAllConds{cc}{hh}(:, 2)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);

        %         % plot ephys and voltage traces
        %         [ax, hl, h2] = plotyy(linspace(0, size(holoSortedImagingMean{cc}(:, hh), 1)/imagingFreq, size(holoSortedImagingMean{cc}(:, hh), 1)), holoSortedImagingMean{cc}(:, hh), linspace(0, size(holoSortedDataMean{cc}(:, hh), 1)/Fs, size(holoSortedDataMean{cc}(:, hh), 1)), holoSortedDataMean{cc}(:, hh));
        % %         set(ax, 'XAxisLocation', 'origin');
        % %         ylim([min(holoSortedImagingMean{cc}(:, hh)), max(holoSortedImagingMean20ƒ{cc}(:, hh))]);
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
        plot(linspace(0, size(holoSortedImagingMean{cc}{hh}, 1)/imagingFreq, size(holoSortedImagingMean{cc}{hh}, 1)), holoSortedImagingMean{cc}{hh}*100, '-', 'linewidth', 2, 'color', 'g');
        ylabel('dF/F (%)');
        xlabel('Time (s)');
        %         ylim([-1 max(holoSortedImagingMean{cc}{hh})*100])
        %         ylim([min(holoSortedImagingMean{cc}{hh})*100 max(holoSortedImagingMean{cc}{hh})*100])
        ylim([-2 5])
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
            xline(nPulseCoords(nn)/Fs, '-', 'LineWidth', 10, 'color', [1 0 0]);
        end
        % plot([0.01; 0.035], [0.5; 0.5], '-k', 'LineWidth', 2);
        % plot([0.01; 0.01], [0.5; 1], '-k', 'LineWidth', 2);
        hold off

        subplot(1,2,2)
        hold on

        fill([linspace(0, size(filtCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(filtCIDffAllConds{cc}{hh}, 1)), fliplr(linspace(0, size(filtCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(filtCIDffAllConds{cc}{hh}, 1)))],...
            [filtCIDffAllConds{cc}{hh}(:, 1)'*100, fliplr(filtCIDffAllConds{cc}{hh}(:, 2)'*100)], [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
        % plot CI lowerbound
        plot(linspace(0, size(filtCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(filtCIDffAllConds{cc}{hh}, 1)), filtCIDffAllConds{cc}{hh}(:, 1)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);
        % plot CI upperbound
        plot(linspace(0, size(filtCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(filtCIDffAllConds{cc}{hh}, 1)), filtCIDffAllConds{cc}{hh}(:, 2)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);

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
        plot(linspace(0, size(filtHoloSortedImagingMean{cc}{hh}, 1)/imagingFreq, size(filtHoloSortedImagingMean{cc}{hh}, 1)), filtHoloSortedImagingMean{cc}{hh}*100, '-', 'linewidth', 2, 'color', 'g');
        ylabel('dF/F (%)');
        xlabel('Time (s)');
        %         ylim([-1 max(holoSortedImagingMean{cc}{hh})*100])
        %         ylim([min(holoSortedImagingMean{cc}{hh})*100 max(holoSortedImagingMean{cc}{hh})*100])
        ylim([-2 5])
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
            xline(nPulseCoords(nn)/Fs, '-', 'LineWidth', 10, 'color', [1 0 0 0]);
        end
        % plot([0.01; 0.035], [0.5; 0.5], '-k', 'LineWidth', 2);
        % plot([0.01; 0.01], [0.5; 1], '-k', 'LineWidth', 2);
        hold off
        pause
    end
end

%% Align Exclusion-applied imaging with ephys traces
nn = double(input('which cell number? '));

exclCIDffAllConds = voltMapping.(cellID{nn}).exclCIDffAllConds;
exclHoloSortedImagingMean = voltMapping.(cellID{nn}).exclHoloSortedImagingMean;
exclFiltCIDffAllConds = voltMapping.(cellID{nn}).exclFiltCIDffAllConds;
exclFiltHoloSortedImagingMean = voltMapping.(cellID{nn}).exclFiltHoloSortedImagingMean;

for cc = 1:nConds
    for hh = 1:nHolos(cc)
        figure(cc*1000+hh);
        % set(gcf, 'Position',  [100, 100, 1600, 800])
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
        ylim([-2 5])
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
   
%% Save Analysis Results
expID = num2str(mouseID);
voltMapping.mouseID = ['voltMapping_Analysis_', expID, '_MultiCellAnalysis'];
directory = '/Volumes/ExData2/Voltage Imaging/VoltMapping/Analysis Results';
directory = '/Users/masatosadahiro/Documents/Data/Voltage Imaging/Voltage Imaging/voltMapping/Analysis Results';
directory = '/Users/masatosadahiro/Documents/Data/Voltage Imaging/Voltage Imaging/voltMapping/Analysis Results/Analysis_newMCanddFFcalc'
% directory = 'E:\Voltage Imaging\VoltMapping\Analysis Results';
fileName = [num2str(voltMapping.mouseID), '.mat'];
save(fullfile(directory, fileName), 'voltMapping', '-v7.3');

TimeNow = datetime('now', 'Format', 'dd-MMM-yyyy HH:mm:ss');
disp(['finished saving at: ' char(TimeNow)])

%% Save Analysis Results (individual cells)
% for nn = 1:nCells
%     expID = num2str(mouseID);
%     voltMapping.mouseID = ['voltMapping_Analysis_', expID, '_MultiCellAnalysis', '_Cell', num2str(nn)];
%     directory = '/Volumes/ExData2/Voltage Imaging/VoltMapping/Analysis Results';
%     directory = '/Users/masatosadahiro/Documents/Data/Voltage Imaging/Voltage Imaging/voltMapping/Analysis Results';
%     fileName = [num2str(voltMapping.mouseID), '.mat'];
%     save(fullfile(directory, fileName), structname{nn}, '-v7.3');
%     TimeNow = datetime('now', 'Format', 'dd-MMM-yyyy HH:mm:ss');
%     disp(['finished saving at: ' char(TimeNow)])
% end

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
% SpotCoordinates = SortedData.holoRequest.targets;
startTime = (voltMapping.holoStimParams.startTime)/1000;
imagesIndex = voltMapping.imagesIndex;
UpOrDown = voltMapping.UpOrDown;
ephysFilePath = voltMapping.ephysFilePath;
ImgsFilePath = voltMapping.ImgsFilePath;
cellID = voltMapping.cellID;

cutOffFreq = 50;   % Cutoff frequency
[b, a] = butter(4, cutOffFreq/(imagingFreq/2));

%% Recreate imgClick and stimtarget 
figure(1010)
set(gcf, 'Position',  [300, 300, 500, 500])
clf
stimTargetFOV = flipdim(Stim_FOV_Full_115, 2);
stimTargetFOV = imrotate(stimTargetFOV, 90);
imagesc(stimTargetFOV);
hold on;
scatter(ExpStruct.holoRequest.actualtargets(:, 1), ExpStruct.holoRequest.actualtargets(:, 2), 'w', 'linewidth', 1);
set(gca,'XTick',[], 'YTick', [])
set(gca,'xdir','reverse','ydir','reverse')
camroll(-270)
axis('square');
hold off
