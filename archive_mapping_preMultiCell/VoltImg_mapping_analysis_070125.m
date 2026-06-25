%% Voltage Imaging Mapping Analysis 070125
% This analysis is used for 2P optogenetic mapping experiments where the readout is either ephys or imaging.

%%
clear all
close all

%% Load files and setup
% Step 1: Read the ephys file
% ephysFilePath = char(uigetdir('E:\Voltage Imaging\voltMapping\Ephys')); % Select and set root folder where all experiments with cells you want to analyze are located
% ephysFilePath = char(uigetdir('D:\Data\Voltage Imaging\Mapping\Data\Ephys Data')); % Select and set root folder where all experiments with cells you want to analyze are located
ephysFilePath = char(uigetdir('/Volumes/ExData2/Voltage Imaging/VoltMapping/Ephys Data')); % Select and set root folder where all experiments with cells you want to analyze are located

ephysFileDir = dir(ephysFilePath);
load([ephysFileDir(end).folder, '/', ephysFileDir(end).name]);
disp(ephysFileDir(end).name);
 
% Step 2: Identify the folder containing the imaging files correspdonding to the ephys file
% ImgsFilePath = char(uigetdir('E:\Voltage Imaging\voltMapping\Imaging')); % Select and set root folder where all experiments with cells you want to analyze are located
% ImgsFilePath = char(uigetdir('D:\Data\Voltage Imaging\Mapping\Data\Imaging Data')); % Select and set root folder where all experiments with cells you want to analyze are located
ImgsFilePath = char(uigetdir('/Volumes/ExData2/Voltage Imaging/VoltMapping/Imaging Data')); % Select and set root folder where all experiments with cells you want to analyze are located

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
pulseDurs = unique(voltMapping.outParams.pulseDur);
nPulses = unique(voltMapping.outParams.nPulses);
ipi = unique(voltMapping.outParams.ipi);
totalPulses = nHolos*nPulses;
% SpotCoordinates = SortedData.holoRequest.targets;
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

voltMapping.imagesIndex = imagesIndex;
voltMapping.imagingFreq = imagingFreq;
voltMapping.UpOrDown = UpOrDown;
voltMapping.ephysFilePath = ephysFilePath;
voltMapping.ImgsFilePath = ImgsFilePath;

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
for tt = 1:size(mappingInputs, 1)
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

voltMapping.excludeTrials = excludeTrials;
voltMapping.ephys.baseVoltAllTrials = baseVoltAllTrials;
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
preStimWindow = 25; % time(ms) window before first pulse 
postStimWindow = 70; % time(ms) to add to stim window after final pulse + ipi

nPulseCoords = []; % indices of when the pulses happen across the stimulation window (in ephys samples)
for ii = 1:nPulses
    nPulseCoords = [nPulseCoords, (ii * ipi/1000*Fs)+preStimWindow/1000*Fs];
end

nPulseCoordsImaging = []; % indices of when the pulses happen across the stimulation window (in imaging frames)
for ii = 1:nPulses
    nPulseCoordsImaging = [nPulseCoordsImaging, (ii * ipi/1000*imagingFreq) + preStimWindow/1000*imagingFreq];
end

condSortedInputs = cell(nConds, 1);
for tt = 1:nTrials
    holoSeqThisTrial = unique(voltMapping.outParams.sequenceThisTrial{tt}, 'stable') - min(unique(voltMapping.outParams.sequenceThisTrial{tt}, 'stable')) + 1;
    holoSeqIndex{voltMapping.trialCond(tt, 1)} = [holoSeqIndex{voltMapping.trialCond(tt, 1)}, holoSeqThisTrial']; % Compile hologram sequences for every trial

    sortingIndex = holoSeqIndex{voltMapping.trialCond(tt, 1)}(:, end)-min(holoSeqIndex{voltMapping.trialCond(tt, 1)})+1;
    if ismember(tt, excludeTrials) 
        holoSeqIndex{voltMapping.trialCond(tt, 1)} = [holoSeqIndex{voltMapping.trialCond(tt, 1)}, NaN(nHolos(voltMapping.trialCond(tt, 1)), 1)]; % Find hologram sequences for every trial
        condSortedInputs{voltMapping.trialCond(tt, 1)} = [condSortedInputs{voltMapping.trialCond(tt, 1)}, NaN(length(mappingInputsBaselined{tt}), 1)];
        for hh = 1:nHolos(voltMapping.trialCond(tt, 1))
            holoSortedDataAllTrials{voltMapping.trialCond(tt, 1)}{sortingIndex(hh), 1} = [holoSortedDataAllTrials{voltMapping.trialCond(tt, 1)}{sortingIndex(hh), 1}, NaN((ipi*nPulses+(preStimWindow+postStimWindow))/1000*Fs+1, 1)];         
        end
        continue
    end
    
    condSortedInputs{voltMapping.trialCond(tt, 1)} = [condSortedInputs{voltMapping.trialCond(tt, 1)}, mappingInputsBaselined{tt}];
    
%     sweepThisTrial =  mappingInputsBaselined{tt}; %Run this line instead of one below if data needs filtering
    sweepThisTrial = filtfilt(blp, alp, mappingInputsBaselined{tt}); %If ephys data needs to be filtered
    
    for hh = 1:nHolos(voltMapping.trialCond(tt, 1)) % NEED TO USE sorting index but...
        thisHoloSweep = sweepThisTrial((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*Fs:((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*Fs+((ipi*nPulses+preStimWindow+postStimWindow)/1000)*Fs)); % Break whole sweep into stimulation-timed windows
        thisHoloSweep = thisHoloSweep - mean(thisHoloSweep(1:nPulseCoords(1)-100)); % Baseline the stimulation-timed window
   
        if max(thisHoloSweep(nPulseCoords:end)) > 0.00% IMPORTANT: comment out this if condition (or lower the threshold to 0) if we want to use all trials instead of just those that actually show some depolarization
            holoSortedDataAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, end), 1} = [holoSortedDataAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, end), 1}, thisHoloSweep];
        else
            holoSortedDataAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, end), 1} = [holoSortedDataAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, end), 1}, NaN(length(thisHoloSweep), 1)];
        end
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

% % Show average traces for each hologram
% holoSortedDataMean = cell(nConds, 1);
% for cc = 1:nConds
%     for hh = 1:nHolos(cc)
%         holoSortedDataMean{cc}(:, hh) = nanmean(holoSortedDataAllTrials{cc}{hh}, 2);
% 
%         % Baseline the mean holo traces
% %         holoSortedDataMean{cc}(:, hh) = holoSortedDataMean{cc}(:, hh) - mean(holoSortedDataMean{cc}((ipi+preStimWindow)/1000*Fs-100:(ipi+preStimWindow)/1000*Fs, hh));
% 
%         figure(cc*100+hh);
%         clf
%         for nn = 1:length(nPulseCoords)
%             xline(nPulseCoords(nn)/Fs, '--', 'LineWidth', 1.5, 'color', [.8 .8 .8]);
%         end
%         hold on
%         plot(linspace(0, length(holoSortedDataMean{cc}(:, hh))/Fs, length(holoSortedDataMean{cc}(:, hh))), holoSortedDataMean{cc}(:, hh), 'LineWidth', 1.5);
%     %     axis off
%         hold off
%         ylim([-0.2, 2]);
%         xlim([0, 0.15]);
%         xlabel('Time(s)')
%         ylabel('mV')
%         % pause
%     end
% end

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

% Combine all hologram traces, calculate one grand mean
holoCombinedDataAllTrials = [];
for cc = 1
    for hh = 1:nHolos(1)
        holoCombinedDataAllTrials = [holoCombinedDataAllTrials, holoSortedDataAllTrials{cc}{hh}(:, :)];
    end
end
holoComboDataGrandMean = nanmean(holoCombinedDataAllTrials, 2);

voltMapping.ephys.holoSeqIndex = holoSeqIndex; 
voltMapping.ephys.holoSortedDataAllTrials = holoSortedDataAllTrials;
voltMapping.ephys.nPulseCoords = nPulseCoords;
voltMapping.ephys.CIephysAllConds = CIephysAllConds;
voltMapping.ephys.holoSortedDataMean = holoSortedDataMean;
voltMapping.ephys.preStimWindow = preStimWindow;
voltMapping.ephys.postStimWindow = postStimWindow;
voltMapping.ephys.condSortedInputs = condSortedInputs;

%% Calculate ROI mask
% Step 1: Calculate mask for ROI 
if length(imagesIndex) < 100
    randTrialsForMask = randperm(length(imagesIndex), length(imagesIndex));
else
    randTrialsForMask = randperm(length(imagesIndex), 100); % Select random n trials/sweeps to be used to generate an average z-stacked image;
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
    currImgPath = [ImgfolderContents(imagesIndex(tt)).folder, '/', ImgfolderContents(imagesIndex(tt)).name];
    info = imfinfo(currImgPath);
    numFrames = numel(info);    

    % Preallocate the image stack
    imageStack = zeros(info(1).Height, info(1).Width, numFrames);
    imageStackMean = zeros(info(1).Height, info(1).Width);
    % Read each frame and store in the stack
    for frameIndex = 1:numFrames
        imageStack(: ,:, frameIndex) = imread(currImgPath, 'Index', frameIndex, 'Info', info);
    end
    imageStackMean = mean(imageStack(:, :, :), 3);
    %     for ff = 1:size(imageStack, 3); figure(101010); imagesc(imageStack(:, :, ff)); caxis([min(min(min(imageStack(:,:,:)))), max(max(max(imageStack(:,:,:))))]); axis equal; axis image; pause; end
    maxDvStack(:, :, tt) = imageStackMean;
end

meanMaxDvStack = mean(maxDvStack(:, :, :), 3); % The average of all the per trial average z-stacks
meanFluorMaxDvStack = meanMaxDvStack; % Grand average z-stack image to use for hand selecting ROI

figure(10); clf; colormap('winter'); imagesc(meanFluorMaxDvStack); axis equal; axis image; colorbar; set(gca, 'fontsize', 12);% caxis([-7 -4]);

% Hand select cell or area of interest, by freehand drawing
roiX = []; roiY = [];
roiHandSelect = drawfreehand;
roiHandSelectMask = createMask(roiHandSelect);
[roiX, roiY] = find(roiHandSelectMask);

% Calculate mean fluorescence of the area of interest and then form into final ROI mask
roiMeanMaxDvStack = zeros(size(maxDvStack, 1), size(maxDvStack, 2));
for rr = 1:length(roiX)
    roiMeanMaxDvStack(roiX(rr), roiY(rr)) = mean(maxDvStack(roiX(rr), roiY(rr),:), 3);
end

figure(11); clf;
imagesc(roiMeanMaxDvStack); axis equal; axis image; colorbar;

stdFluor = std(nonzeros(roiMeanMaxDvStack));
meanFluor = mean(nonzeros(roiMeanMaxDvStack));

% Designate cutoff fluorescence for pixels to be selected for ROI
cutOffFluor = meanFluor; %stdFluor*1 + meanFluor; % currently cutoff is 1 standard devs from mean fluorescence

roiMeanMaxDvStack(roiMeanMaxDvStack <= cutOffFluor) = 0;
roiMeanMaxDvStack(roiMeanMaxDvStack > 0) = 1;
roiStack = roiMeanMaxDvStack;

figure(12); clf;
imagesc(roiStack); axis equal; axis image;
[roiX, roiY] = find(roiStack); % XY coordinates of the pixels of interest

voltMapping.maxDvStack = maxDvStack;
voltMapping.meanFluorMaxDvStack = meanFluorMaxDvStack;
voltMapping.roiMeanMaxDvStack = roiMeanMaxDvStack;
voltMapping.cutOffFluor = cutOffFluor;
voltMapping.roiX = roiX;
voltMapping.roiY = roiY;

%% Calculate df and mean df for all holos sorted by trial conditions
% Use these lines if I need to further filter the df or dff trace
cutOffFreq = 50;   % Cutoff frequency
[b, a] = butter(4, cutOffFreq/(imagingFreq/2));  % 4th order Butterworth filter

% Preallocate for data sorted by holograms, all trials, across conditions
holoSortedImagingAllTrials = cell(nConds, 1);
filtHoloSortedImagingAllTrials = cell(nConds, 1);
for cc = 1:nConds
    holoSortedImagingAllTrials{cc} = cell(nHolos(cc), 1);
    filtHoloSortedImagingAllTrials{cc} = cell(nHolos(cc), 1);
end

counter = 0;
for tt = 1:nTrials %size(vsTest_inputs, 2)
    counter = counter+1;
    disp(['Trial number: ', num2str(counter)]);
%     if ismember(tt, excludeTrials)
%         continue
%     end
    
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
    
    % Baseline calculation
    % Fluorescence traces of pre-stim and post-stim periods
    startTimeImaging = floor(startTime*imagingFreq);
    postStimSeqTime = voltMapping.holoStimParams.postStimSeqTime;
    postStimSeqTimeImaging = ceil(postStimSeqTime/1000*imagingFreq);
    
    baselinePreImageStack = imageStack(roiX, roiY, 1:startTimeImaging);
    baselinePostImageStack = imageStack(roiX, roiY, (size(imageStack, 3)-postStimSeqTimeImaging):size(imageStack, 3));
    
    roiBaselinePreImageStack = [];
    roiBaselinePostImageStack = [];
    for ff = 1:size(baselinePreImageStack, 3)
        roiBaselinePreImageStack(ff, 1) = mean(mean(baselinePreImageStack(:, :, ff)));
    end
    for ff = 1:size(baselinePostImageStack, 3)
        roiBaselinePostImageStack(ff, 1) = mean(mean(baselinePostImageStack(:, :, ff)));
    end
    
    % Moving window variance during pre-stim period
    windowLimits = [1/imagingFreq, size(roiBaselinePreImageStack, 1)/imagingFreq]; % in sec The range within each sweep where variance is sampled by 10 sample windows, in ms. Will leave it as 500 ms total width, starting from 110ms after start (after the step pulse).
    firstLimit = windowLimits(1)*imagingFreq; % Start of the first sampling window (in samples).
    windowTime = ceil(size(roiBaselinePreImageStack, 1)/imagingFreq*1000)/1000;
    segmentTime = 0.050; % time length of each window (s)
    numSegments = windowTime/segmentTime; % number of sample windows within the limit
    windowWidth = floor((windowLimits(2) - windowLimits(1))/numSegments*imagingFreq); % Span of each sample window (in samples).
    
    varBaselinePre = movvar(roiBaselinePreImageStack, windowWidth); % moving window variance calculated across prestimulus 
    fanoBaselinePre = varBaselinePre/(mean(roiBaselinePreImageStack));
    % varBaselinePost = movvar(roiBaselinePostImageStack, windowWidth); % moving window variance calculated across poststim
    % fanoBaselinePost = varBaselinePost/(mean(roiBaselinePostImageStack));

    % Find moments of var/fano below threshold
    % Method 1: Threshold set at standard deviation from mean of entire period before stim begins
    % varThresholdPre = mean(varBaselinePre) - std(varBaselinePre); % variance threshold (mean-std) 
    % fanoThresholdPre = mean(fanoBaselinePre)-std(fanoBaselinePre); % fanofactor (variance/mean)
    % varThresholdPost = mean(varBaselinePost) - std(varBaselinePost);
    % fanoThresholdPost = mean(fanoBaselinePost)-std(fanoBaselinePost);

    % Method 2: Take quantile bottom 10% of fanofactor
    % Step 1: Define quantile threshold (e.g., bottom 10%)
    q = 0.10;  % Change to 0.05 for bottom 5%, etc.
    quantileCutoff = quantile(x, q);
    % Step 2: Select low points
    fanoLowestPre = find(x < quantileCutoff);
    roiBaselinePre = mean(roiBaselinePreImageStack(fanoLowestPre));

    % Calculate mean of fluorescence during moments of fanofactor below threshold, to be designated as baseline
    % [vLowestPre, ~] = find(varBaselinePre < varThresholdPre); % moments of variance below threshold
    % [fanoLowestPre, ~] = find(fanoBaselinePre < fanoThresholdPre); % moments of fanofactor below threshold
    % roiBaselinePre = mean(roiBaselinePreImageStack(fanoLowestPre)); % fluoresence where fanofactor is below thershold
    roiBaselineMean = roiBaselinePre; % mean([roiBaselinePre, roiBaselinePost]); % grand mean fluoresence of all the moments where fanofactor is below threshold

    % Break apart the imageStack for this trial into stim windows and rearrange according to hologram sequence 
    holoSeqThisTrial = (unique(voltMapping.outParams.sequenceThisTrial{tt}, 'stable') - min(unique(voltMapping.outParams.sequenceThisTrial{tt}, 'stable')) + 1)'; % Hologram sequence for this trial    
    for hh = 1:nHolos(voltMapping.trialCond(tt, 1))
        % Extract the frames associated with the current hologram, including pre and post stimulation windows
%         framesThisHolo = imageStack(:, :, ceil((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*imagingFreq):(ceil(voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)*imagingFreq)+ceil((ipi*nPulses+postStimWindow)/1000*imagingFreq)));
        roiFramesThisHolo = imageStack(roiX, roiY, ceil((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*imagingFreq):(ceil((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*imagingFreq)+ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)));

        % Set baseline(f0), this is based on a single mean across the preStimWindow, so it may be a flawed approach!
        % roiBaselineAllPixels = framesThisHolo(roiX, roiY, 1:ceil(preStimWindow/1000*imagingFreq));
        % roiBaselineMean = nanmean(roiBaselineAllPixels, 'all');
        
        dfThisHolo = [];
        dffThisHolo = [];
        for ff = 1:size(roiFramesThisHolo, 3)
            % ROI pixels for this frame
            currFrameRoi = roiFramesThisHolo(:, :, ff);
            currFrameRoiMean = nanmean(currFrameRoi, 'all'); % The mean across all pixels in the ROI for the select frame
            
            % Calculate df
            intensityChange = currFrameRoiMean - roiBaselineMean; % essentially, intensityChange = df, and dff is intensityChnage/roiBaselineMean
            dfThisHolo = [dfThisHolo; intensityChange];
            dffThisHolo = [dffThisHolo; intensityChange/roiBaselineMean];
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
            holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, NaN(ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)+1, 1)];
            filtHoloSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [filtHoloSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, NaN(ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)+1, 1)];
            
        else
%           holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, tt), 1}(:, tt) = dfThisHolo;
            holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, dffThisHolo];
            filtHoloSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [filtHoloSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, filtdffThisHolo];

        end
    end
end

%% Calculate mean response (and CI) for each hologram across trials and per condition

holoSortedImagingMean = cell(nConds, 1); % The mean response for each hologram across conditions
filtHoloSortedImagingMean = cell(nConds, 1);
for cc = 1:nConds
    holoSortedImagingMean{cc} = cell(nHolos(cc), 1);
    filtHoloSortedImagingMean{cc} = cell(nHolos(cc), 1);
end
for cc = 1:nConds
    for hh = 1:nHolos(cc)
        holoSortedImagingMean{cc}{hh} = nanmean(holoSortedImagingAllTrials{cc}{hh}, 2);
        filtHoloSortedImagingMean{cc}{hh} = nanmean(filtHoloSortedImagingAllTrials{cc}{hh}, 2);        
        % Baseline the mean holo traces
%         holoSortedImagingMean{cc}{hh} = holoSortedImagingMean{cc}{hh} - mean(holoSortedImagingMean{cc}{hh}(1:preStimWindow/1000*imagingFreq));
%         filtHoloSortedImagingMean{cc}{hh} = filtHoloSortedImagingMean{cc}{hh} - mean(filtHoloSortedImagingMean{cc}{hh}(1:preStimWindow/1000*imagingFreq));
      
    %     figure(30);
    %     clf
    %     for nn = 1:length(nPulseCoordsImaging)
    %         xline(nPulseCoordsImaging(nn)/imagingFreq, '--', 'LineWidth', 1.5, 'color', [.8 .8 .8]);
    %     end
    %     hold on
    %     plot(linspace(0, length(holoSortedImagingMean{cc}{hh})/imagingFreq, length(holoSortedImagingMean{cc}{hh})), holoSortedImagingMean{cc}{hh}, 'LineWidth', 1.5);
    %     plot(linspace(0, length(filtHoloSortedImagingMean{cc}{hh})/imagingFreq, length(filtHoloSortedImagingMean{cc}{hh})), filtHoloSortedImagingMean{cc}{hh}, 'LineWidth', 1.5);
    % %     axis off
    %     hold off
    %     xlabel('Time(s)')
    %     ylabel('df/f')
        % pause
    end
end

% Below is an alternative step I am experimenting with:
% Reject all Imaging Trials where cooresponding ephys traces were rejected by thresholding.
% Reminder: There's no point in rejecting imaging trials with thresholding
% because what the electrode reads (no matter the access resistance) isn't
% necessarily the same as what the GEVI reads out. The GEVI in theory should read
% out depolarization more consistently.
holoSortedImagingAllTrials_ALT = holoSortedImagingAllTrials;
filtHoloSortedImagingAllTrials_ALT = filtHoloSortedImagingAllTrials;
for cc = 1:nConds
    for hh = 1:nHolos(cc)
        for tt = 1:size(holoSortedDataAllTrials{cc}{hh}, 2)
            if isnan(holoSortedDataAllTrials{cc}{hh}(:, tt))
                holoSortedImagingAllTrials_ALT{cc}{hh}(:, tt) = NaN;
                filtHoloSortedImagingAllTrials_ALT{cc}{hh}(:, tt) = NaN;
            end
        end
        nonNanColumns = ~any(isnan(holoSortedImagingAllTrials_ALT{cc}{hh}), 1);
        holoSortedImagingAllTrials_ALT{cc}{hh} = holoSortedImagingAllTrials_ALT{cc}{hh}(:, nonNanColumns);
        filtHoloSortedImagingAllTrials_ALT{cc}{hh} = filtHoloSortedImagingAllTrials_ALT{cc}{hh}(:, nonNanColumns);        
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

CIDffAllConds = cell(nConds, 1);
filtCIDffAllConds = cell(nConds, 1);
for cc = 1:nConds 
    for hh = 1:nHolos(cc)
        confidence_level = 0.95;
        means = nanmean(holoSortedImagingAllTrials{cc}{hh, 1}, 2);
        filtMeans = nanmean(filtHoloSortedImagingAllTrials{cc}{hh, 1}, 2);       
        std_errors = std(holoSortedImagingAllTrials{cc}{hh, 1}, 0, 2, "omitnan") / sqrt(size(holoSortedImagingAllTrials{cc}{hh, 1}, 2));
        filtStd_errors = std(filtHoloSortedImagingAllTrials{cc}{hh, 1}, 0, 2, "omitnan") / sqrt(size(filtHoloSortedImagingAllTrials{cc}{hh, 1}, 2));
  
        t_score = tinv((1 + confidence_level) / 2, size(holoSortedImagingAllTrials{cc}{hh, 1}, 2) - 1);
        filtT_score = tinv((1 + confidence_level) / 2, size(filtHoloSortedImagingAllTrials{cc}{hh, 1}, 2) - 1);
        margin_of_error = t_score * std_errors;
        filtMargin_of_error = filtT_score * filtStd_errors;
        lower_bounds = means - margin_of_error;
        filtLower_bounds = filtMeans - filtMargin_of_error;
        upper_bounds = means + margin_of_error;
        filtUpper_bounds = filtMeans + filtMargin_of_error;
        if UpOrDown == '2'
            CIDffAllConds{cc}{hh, 1} = [lower_bounds, upper_bounds];
            filtCIDffAllConds{cc}{hh, 1} = [filtLower_bounds, filtUpper_bounds];
        elseif UpOrDown =='1'
            CIDffAllConds{cc}{hh, 1} = [-lower_bounds, -upper_bounds];
            filtCIDffAllConds{cc}{hh, 1} = [-filtLower_bounds, -filtUpper_bounds];

        end
    end
end

voltMapping.nPulseCoordsImaging                = nPulseCoordsImaging;
voltMapping.holoSortedImagingAllTrials         = holoSortedImagingAllTrials;
voltMapping.filtHoloSortedImagingAllTrials     = filtHoloSortedImagingAllTrials;
voltMapping.holoSortedImagingAllTrials_ALT     = holoSortedImagingAllTrials_ALT;
voltMapping.filtHoloSortedImagingAllTrials_ALT = filtHoloSortedImagingAllTrials_ALT;
voltMapping.holoSeqIndex                       = holoSeqIndex;
voltMapping.holoSortedImagingMean              = holoSortedImagingMean;
voltMapping.filtHoloSortedImagingMean          = filtHoloSortedImagingMean;
% voltMapping.holoComboImagingAllTrials        = holoComboImagingAllTrials;
% voltMapping.holoComboImagingGrandMean        = holoComboImagingGrandMean;
voltMapping.CIDffAllConds                      = CIDffAllConds;
voltMapping.filtCIDffAllConds                  = filtCIDffAllConds;
voltMapping.preStimWindow                      = preStimWindow;
voltMapping.postStimWindow                     = postStimWindow;

%% Align imaging with ephys traces
for cc = 1:nConds
    for hh = 1:nHolos(cc)
        figure(cc*1000+hh);
%         set(gcf, 'Position',  [100, 100, 600, 400])
        clf
        hold on
%         fill([linspace(0, size(CIDfAllConds{cc}{hh}, 1)/imagingFreq, size(CIDfAllConds{cc}{hh}, 1)), fliplr(linspace(0, size(CIDfAllConds{cc}{hh}, 1)/imagingFreq, size(CIDfAllConds{cc}{hh}, 1)))],...
%         [CIDfAllConds{cc}{hh}(:, 1)', fliplr(CIDfAllConds{cc}{hh}(:, 2)')], [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
%         % plot CI lowerbound
%         plot(linspace(0, size(CIDfAllConds{cc}{hh}, 1)/imagingFreq, size(CIDfAllConds{cc}{hh}, 1)), CIDfAllConds{cc}{hh}(:, 1), '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]); 
%         % plot CI upperbound
%         plot(linspace(0, size(CIDfAllConds{cc}{hh}, 1)/imagingFreq, size(CIDfAllConds{cc}{hh}, 1)), CIDfAllConds{cc}{hh}(:, 2), '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);   

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
%         yyaxis right
%         plot(linspace(0, size(holoSortedImagingMean{cc}{hh}, 1)/imagingFreq, size(holoSortedImagingMean{cc}{hh}, 1)), holoSortedImagingMean{cc}{hh}, '-', 'linewidth', 2, 'color', 'g');
%         plot(linspace(0, size(holoSortedImagingMean{cc}{hh}, 1)/imagingFreq, size(holoSortedImagingMean{cc}{hh}, 1)), filter(b, a, holoSortedImagingMean{cc}{hh})*100, '-', 'linewidth', 2, 'color', 'g');
        plot(linspace(0, size(filtHoloSortedImagingMean{cc}{hh}, 1)/imagingFreq, size(filtHoloSortedImagingMean{cc}{hh}, 1)), filtHoloSortedImagingMean{cc}{hh}*100, '-', 'linewidth', 2, 'color', 'g');
        ylabel('dF/F (%)');
        xlabel('Time (s)');
%         ylim([-1 max(holoSortedImagingMean{cc}{hh})*100])
%         ylim([min(holoSortedImagingMean{cc}{hh})*100 max(holoSortedImagingMean{cc}{hh})*100])
        ylim([-1.5 3])
        ax.YColor = [0 1 0];
        % xticks([0:0.05:size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);

        % axis off
%         
        % plot ephys trace only
%         yyaxis left
%     %     axes('Position',[.70 .12 .2 .2]);
%     %     box on
%         plot(linspace(0, size(holoSortedDataMean{cc}(:, hh), 1)/Fs, size(holoSortedDataMean{cc}(:, hh), 1)), holoSortedDataMean{cc}(:, hh), 'linewidth', 1.5, 'color', [0 0 0]);
%         gca;
%         set(gca,'xtick',[], 'fontsize', 18);
% %         ylim([min(holoSortedDataMean{cc}(:, hh)), max(holoSortedDataMean{cc}(:, hh))]);
%         ylim([-0.5 2]);
%         ylabel('dV');
%         ax.YColor = [0, 0, 0];
        % axis off
%         xlim([0 size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);
%         xticks([0:0.05:size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);
%         
        % show line at dff = 0
%         yline(0, '-', 'LineWidth', 1.5, 'color', [0.9 0.9 0.9]);
        for nn = 1:length(nPulseCoords)
            xline(nPulseCoords(nn)/Fs, '-', 'LineWidth', 10, 'color', [1 0 0]);
        end
        
        % plot([0.01; 0.035], [0.5; 0.5], '-k', 'LineWidth', 2);
        % plot([0.01; 0.01], [0.5; 1], '-k', 'LineWidth', 2);
        hold off
     pause
    end
end

%% Save Analysis Results
expID = num2str(mouseID);
voltMapping.mouseID = ['voltMapping_Analysis_', expID, '_newBaselining'];
directory = '/Volumes/ExData2/Voltage Imaging/VoltMapping/Analysis Results';
fileName = [num2str(voltMapping.mouseID), '.mat'];
save(fullfile(directory, fileName), 'voltMapping', '-v7.3');

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

cutOffFreq = 50;   % Cutoff frequency
[b, a] = butter(4, cutOffFreq/(imagingFreq/2));

%% Recreate imgClick and stimtarget 
figure(1010)
stimTargetFOV = flipdim(imgClick{1, 1}, 2);
stimTargetFOV = imrotate(stimTargetFOV, 90);
imagesc(stimTargetFOV);
hold on;
scatter(ExpStruct.holoRequest.actualtargets(:, 1), ExpStruct.holoRequest.actualtargets(:, 2), 'w');
set(gca,'XTick',[], 'YTick', [])
set(gca,'xdir','reverse','ydir','reverse')
camroll(-270)
axis('square');