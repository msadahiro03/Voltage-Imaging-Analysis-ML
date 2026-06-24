%% Voltage Imaging Mapping (VIM) Analysis
% This analysis is used for 2P optogenetic mapping experiments where the readout is either ephys or imaging.

%%
clear all
close all

%% Load files and setup
% Step 1: Read the ephys file
% ephysFilePath = char(uigetdir('E:\Data\Voltage Imaging\voltMapping\Ephys Data')); % Select and set root folder where all experiments with cells you want to analyze are located
ephysFilePath = char(uigetdir('D:\Data\Voltage Imaging\Mapping\Data\Ephys Data')); % Select and set root folder where all experiments with cells you want to analyze are located
% ephysFilePath = char(uigetdir('/Volumes/Seagate Backup Plus Drive/New folder/240925')); % Select and set root folder where all experiments with cells you want to analyze are located

ephysFileDir = dir(ephysFilePath);
load([ephysFileDir(end).folder, '/', ephysFileDir(end).name]);
disp(ephysFileDir(end).name);
 
% Step 2: Identify the folder containing the imaging files correspdonding to the ephys file
% % ImgsFilePath = char(uigetdir('E:\Data\Voltage Imaging\voltMapping\Imaging Data')); % Select and set root folder where all experiments with cells you want to analyze are located
ImgsFilePath = char(uigetdir('D:\Data\Voltage Imaging\Mapping\Data\Imaging Data')); % Select and set root folder where all experiments with cells you want to analyze are located
% ImgsFilePath = char(uigetdir('/Volumes/Seagate Backup Plus Drive/New folder/092524')); % Select and set root folder where all experiments with cells you want to analyze are located

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

imagingFreq = ExpStruct.sampleFreq;

% Step 3: Setup struct
voltMapping = ExpStruct;

if exist('ExpStruct2', 'var') % If the experiment has a second patch electrode
    voltMapping.ExpStruct2 = ExpStruct2;
end

% Stimulation properties
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
directory = 'D:\Data\Voltage Imaging\voltMapping\Analysis Results\';
fileName = ['voltMapping ', num2str(ExpStruct.mouseID), '.mat'];

lpCut =2000; % filtering data params
[blp,alp] = butter(4, [lpCut/Fs],'low');

voltMapping.imagesIndex = imagesIndex;
voltMapping.imagingFreq = imagingFreq;
voltMapping.UpOrDown = UpOrDown;
voltMapping.ephysFilePath = ephysFilePath;
voltMapping.ImgsFilePath = ImgsFilePath;

%% Trial exclusion and baselining
% Excludes all trials where the baseline Vm suddenly jumps (losing cell etc.)
vThreshold = -50; %set threshold Vm here
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

voltMapping.excludeTrials = excludeTrials;
voltMapping.ephys.baseVoltAllTrials = baseVoltAllTrials;
voltMapping.ephys.mappingInputsBaselined = mappingInputsBaselined;

%% Break apart each continuous ephys trace into stim windows and rearrange according to hologram sequence
holoSeqIndex = cell(nConds, 1);
holoSortedDataAllTrials = cell(nConds, 1);
for cc = 1:nConds
    holoSortedDataAllTrials{cc} = cell(nHolos(cc), 1);
end

preStimWindow = 24;
postStimWindow = 40; % time(ms) to add to stim window after final pulse and ipi

for tt = 1:nTrials
    if ismember(tt, excludeTrials)
        holoSeqIndex{voltMapping.trialCond(tt, 1)} = [holoSeqIndex{voltMapping.trialCond(tt, 1)}, NaN(nHolos(voltMapping.trialCond(tt, 1)), 1)]; % Find hologram sequences for every trial
        for hh = 1:nHolos(voltMapping.trialCond(tt, 1))
            holoSortedDataAllTrials{voltMapping.trialCond(tt, 1)}{hh, 1}(:, tt) = NaN(ipi*nPulses/1000*Fs+1, 1);
        end
        continue
    end

    holoSeqIndex{voltMapping.trialCond(tt, 1)} = [holoSeqIndex{voltMapping.trialCond(tt, 1)}, unique(voltMapping.outParams.sequenceThisTrial{tt}, 'stable')']; % Find hologram sequences for every trial
    sortingIndex = 1:length(holoSeqIndex{voltMapping.trialCond(tt, 1)}(:, end));
    for hh = 1:nHolos(voltMapping.trialCond(tt, 1))
        thisHoloSweep = mappingInputsBaselined{tt}((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*Fs:((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)*Fs)+((ipi*nPulses/1000+postStimWindow/1000)*Fs))); % Break whole sweep into stimulation-timed windows
        
        thisHoloSweep = thisHoloSweep - mean(thisHoloSweep((ipi+preStimWindow)/1000*Fs-100:(ipi+preStimWindow)/1000*Fs)); % Baseline the stimulation-timed window
%         holoSortedDataAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, tt), 1}(:, tt) = thisHoloSweep;
%         holoSortedDataAllTrials{voltMapping.trialCond(tt, 1)}{sortingIndex(hh), 1}(:, tt) = thisHoloSweep;
        if max(thisHoloSweep) > 1.25 % IMPORTANT: comment out this if condition (or lower the threshold to 0) if we want to use all trials instead of just those that actually show some depolarization
            holoSortedDataAllTrials{voltMapping.trialCond(tt, 1)}{sortingIndex(hh), 1} = [holoSortedDataAllTrials{voltMapping.trialCond(tt, 1)}{sortingIndex(hh), 1}, thisHoloSweep];
        else
            holoSortedDataAllTrials{voltMapping.trialCond(tt, 1)}{sortingIndex(hh), 1} = [holoSortedDataAllTrials{voltMapping.trialCond(tt, 1)}{sortingIndex(hh), 1}, NaN(length(thisHoloSweep), 1)];
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

% % NEEDS FIXING: Remove trials where dv goes below zero (not sure why some trials do that  but those end up messing up the mean)
% for pp = 1:nConds
%     for hh = 1:nHolos
%         for tt = 1:nTrials
% %             if mean(holoSortedDataAllTrials{pp}{hh}(:, tt)) < 0.1
% %                 holoSortedDataAllTrials{pp}{hh}(:, tt) = NaN;
% %             end
% 
%             percentNegative = sum(holoSortedDataAllTrials{pp}{hh}(:, tt) < 0) / numel(holoSortedDataAllTrials{pp}{hh}(:, tt)) * 100;
%             percentThresh = 50;
%             if percentNegative >= percentThresh
%                 holoSortedDataAllTrials{pp}{hh}(:, tt) = NaN;
%             end
%         end
%     end
% end

nPulseCoords = [];
for ii = 1:nPulses
    nPulseCoords = [nPulseCoords, (ii * ipi/1000*Fs)+preStimWindow/1000*Fs];
end

nPulseCoordsImaging = [];
for ii = 1:nPulses
    nPulseCoordsImaging = [nPulseCoordsImaging, (ii * ipi/1000*imagingFreq) + preStimWindow/1000*imagingFreq];
end

holoSortedDataMean = cell(nConds, 1);
for cc = 1:nConds
    for hh = 1:nHolos(cc)
        holoSortedDataMean{cc}(:, hh) = nanmean(holoSortedDataAllTrials{cc}{hh}, 2);
        
        % Baseline the mean holo traces
        holoSortedDataMean{cc}(:, hh) = holoSortedDataMean{cc}(:, hh) - mean(holoSortedDataMean{cc}((ipi+preStimWindow)/1000*Fs-100:(ipi+preStimWindow)/1000*Fs, hh));
        
        figure(1);
        clf
        for nn = 1:length(nPulseCoords)
            xline(nPulseCoords(nn)/Fs, '--', 'LineWidth', 1.5, 'color', [.8 .8 .8]);
        end
        hold on
        plot(linspace(0, length(holoSortedDataMean{cc}(:, hh))/Fs, length(holoSortedDataMean{cc}(:, hh))), filtfilt(blp,alp,holoSortedDataMean{cc}(:, hh)), 'LineWidth', 1.5);
    %     axis off
        hold off
        xlabel('Time(s)')
        ylabel('mV')
        pause
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
voltMapping.ephys.holoSortedDataMean = holoSortedDataMean;
voltMapping.ephys.preStimWindow = preStimWindow;
voltMapping.ephys.postStimWindow = postStimWindow;

%% Calculate ROI mask
% Step 1: Calculate mask for ROI 
maxDvStack = [];
for tt = 1:150 %length(imagesIndex)
    if ismember(imagesIndex(tt), excludeTrials)
        continue
    end
    currImgPath = [ImgfolderContents(imagesIndex(tt)).folder, '/', ImgfolderContents(imagesIndex(tt)).name];
    info = imfinfo(currImgPath);
    numFrames = numel(info);    

    % Preallocate the image stack
    imageStack = zeros(info(1).Height, info(1).Width, numFrames);

    % Read each frame and store in the stack
    for frameIndex = 1:numFrames
        imageStack(:,:,frameIndex) = imread(currImgPath, 'Index', frameIndex, 'Info', info);
    end
%     for ff = 1:size(imageStack, 3); figure(101010); imagesc(imageStack(:, :, ff)); caxis([min(min(min(imageStack(:,:,:)))), max(max(max(imageStack(:,:,:))))]); axis equal; axis image; pause; end
    maxDvStack = cat(3, maxDvStack, imageStack);
end

meanMaxDvStack = mean(maxDvStack(:, :, :), 3);
meanFluorMaxDvStack = meanMaxDvStack;

figure(10); clf; imagesc(meanFluorMaxDvStack); axis equal; axis image; colorbar; set(gca, 'fontsize', 12);% caxis([-7 -4]);

% Old Method: Automatically set cutoff from whole mean
% meanMaxDvStack(meanMaxDvStack <= cutOffFluor) = 0;
% meanMaxDvStack(meanMaxDvStack > 0) = 1;
% roiStack = meanMaxDvStack;
% figure(20); imagesc(roiStack); axis equal; axis image;
% [roiX, roiY] = find(roiStack);
% imagesc(meanMaxDvStack);
% axis equal
% axis image

%%%%%
% New Method: Hand select cell or area of interest
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

% Old version: Calculate mean fluorescence of the area of interest and then form into final ROI mask
% roiMeanMaxDvStack = mean(maxDvStack(min(roiX):max(roiX), roiY, :), 3);
% imagesc(roiMeanMaxDvStack)

stdFluor = std(nonzeros(roiMeanMaxDvStack));
meanFluor = mean(nonzeros(roiMeanMaxDvStack));

% Designate cutoff fluorescence for pixels to be selected for ROI
cutOffFluor = meanFluor; %stdFluor*1+ meanFluor; % currently cutoff is 1 standard devs from mean fluorescence

roiMeanMaxDvStack(roiMeanMaxDvStack <= cutOffFluor) = 0;
roiMeanMaxDvStack(roiMeanMaxDvStack > 0) = 1;
roiStack = roiMeanMaxDvStack;

figure(12); clf;
imagesc(roiStack); axis equal; axis image;
[roiX, roiY] = find(roiStack);

voltMapping.meanFluorMaxDvStack = meanFluorMaxDvStack;
voltMapping.roiMeanMaxDvStack = roiMeanMaxDvStack;
voltMapping.cutOffFluor = cutOffFluor;
voltMapping.roiX = roiX;
voltMapping.roiY = roiY;

%% Calculate df and mean df for all conditions
% Step 1: Calculate df across all trials across all conditions
dfAllTrialsCondSorted = cell(nConds, 1);
dfAllTrialsNonSorted = cell(nTrials, 1);
traceAllConds = cell(nConds, 1);
f0Start =  0.01; % in seconds
f0End = startTime;

for tt = 1:nTrials %size(vsTest_inputs, 2)
    if ismember(tt, excludeTrials)
        continue
    end
    
    % Read the multi-frame image
    currImgPath = [ImgfolderContents(imagesIndex(tt)).folder, '/', ImgfolderContents(imagesIndex(tt)).name];
    info = imfinfo(currImgPath);
    numFrames = numel(info);

    % Preallocate the image stack
    imageStack = zeros(info(1).Height, info(1).Width, numFrames);
    
    % Read each frame and store in the stack
    for frameIndex = 1:numFrames
        imageStack(:,:,frameIndex) = imread(currImgPath, 'Index', frameIndex, 'Info', info);
    end
     
%     for tt = 1:size(imageStack, 3); figure(101010); imagesc(imageStack(:, :, tt)); caxis([min(min(min(imageStack(:,:,:)))), max(max(max(imageStack(:,:,:))))]); axis equal; axis image; pause; end
   
%    Extract the pixel intensity values (new option 2)    
    df = [];
    roiPixels = [];
    previousRoiPixels = [];
    for ff = 2:size(imageStack, 3)
        % Read current frame
        currentFrame = imageStack(:, :, ff); 
        % Read previous frame
        previousFrame = imageStack(:, :, ff-1);

            for rr = 1:size(roiX, 1)
                roiPixels(rr) = currentFrame(roiX(rr), roiY(rr));
                previousRoiPixels(rr) = previousFrame(roiX(rr), roiY(rr));
            end

        intensityChange = mean(abs(roiPixels(:) - previousRoiPixels(:)));
        df = [df; intensityChange];
    end

    df = (df - mean(df(ceil(imagingFreq*f0Start):floor(imagingFreq*f0End))))/mean(df(ceil(imagingFreq*f0Start):floor(imagingFreq*f0End)))*100;

    if UpOrDown == '2'
        dfAllTrialsCondSorted{voltMapping.trialCond(tt, 1), 1} = [dfAllTrialsCondSorted{voltMapping.trialCond(tt, 1), 1}, -df];
        dfAllTrialsNonSorted{tt} = -df;
    elseif UpOrDown =='1'
        dfAllTrialsCondSorted{voltMapping.trialCond(tt, 1), 1} = [dfAllTrialsCondSorted{voltMapping.trialCond(tt, 1), 1}, df];
        dfAllTrialsNonSorted{tt} = df;
    end

end

voltMapping.dfAllTrialsCondSorted = dfAllTrialsCondSorted;
voltMapping.dfAllTrialsNonSorted = dfAllTrialsNonSorted;

%% Break apart trial-by-trial imaging data into stim windows and rearrange according to hologram sequence
holoSortedImagingAllTrials = cell(nConds, 1);
for cc = 1:nConds
    holoSortedImagingAllTrials{cc} = cell(nHolos(cc), 1);
end

for tt = 1:nTrials
    sortingIndex = [];
    sortingIndex = 1 : length(holoSeqIndex{voltMapping.trialCond(tt, 1)}(:, end));   
    for hh = 1:nHolos(voltMapping.trialCond(tt, 1))
        if isempty(dfAllTrialsNonSorted{tt})
            dfThisHolo = NaN(ceil(ipi*nPulses/1000*imagingFreq)+1, 1);
        else
            dfThisHolo = dfAllTrialsNonSorted{tt}(ceil((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*imagingFreq):(ceil(voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)*imagingFreq)+ceil((ipi*nPulses+postStimWindow)/1000*imagingFreq)));
            dfThisHolo = dfThisHolo - mean(dfThisHolo((ipi+preStimWindow)/1000*imagingFreq-2:(ipi+preStimWindow)/1000*imagingFreq));
        end
        
        if ismember(tt, excludeTrials) % isnan(holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, tt))
%             holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{hh, 1}(:, tt) = NaN(ceil(ipi*nPulses/1000*imagingFreq)+1, 1);
        holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{sortingIndex(hh), 1} = [holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{sortingIndex(hh), 1}, NaN(ceil(ipi*nPulses/1000*imagingFreq)+1, 1)];
        else
%         holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, tt), 1}(:, tt) = dfThisHolo;
        
        holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{sortingIndex(hh), 1} = [holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{sortingIndex(hh), 1}, dfThisHolo];

        end
    end
end

% NEEDS FIXING: If strictly NaN-ing all holo stim trials where ephys goes negative
% holoSortedImagingAllTrials_ALT = holoSortedImagingAllTrials;
% for cc = 1:nConds
%     for hh = 1:nHolos
%         for tt = 1:nTrials
%             if isnan(holoSortedDataAllTrials{cc}{hh}(:, tt))
%                 holoSortedImagingAllTrials_ALT{cc}{hh}(:, tt) = NaN;
%             end
%         end
%     end
% end

% Reject all Imaging Trials where cooresponding ephys traces were rejected by thresholding.
% Reminder: There's no point in rejecting imaging trials with thresholding
% because what the electrode reads (no matter the access resistance) isn't
% necessarily the same as what the GEVI reads out. The GEVI in theory should read
% out depolarization more consistently.
holoSortedImagingAllTrials_ALT = holoSortedImagingAllTrials;
for cc = 1:nConds
    for hh = 1:nHolos(cc)
        for tt = 1:size(holoSortedDataAllTrials{cc}{hh}, 2)
            if isnan(holoSortedDataAllTrials{cc}{hh}(:, tt))
                holoSortedImagingAllTrials_ALT{cc}{hh}(:, tt) = NaN;
            end
        end
        nonNanColumns = ~any(isnan(holoSortedImagingAllTrials_ALT{cc}{hh}), 1);
        holoSortedImagingAllTrials_ALT{cc}{hh} = holoSortedImagingAllTrials_ALT{cc}{hh}(:, nonNanColumns);
    end
end  


holoSortedImagingMean = cell(nConds, 1);
for cc = 1:nConds
    for hh = 1:nHolos(cc)
        holoSortedImagingMean{cc}(:, hh) = nanmean(holoSortedImagingAllTrials_ALT{cc}{hh}, 2);
        
        % Baseline the mean holo traces
        holoSortedImagingMean{cc}(:, hh) = holoSortedImagingMean{cc}(:, hh) - mean(holoSortedImagingMean{cc}((ipi+preStimWindow)/1000*imagingFreq-2:(ipi+preStimWindow)/1000*imagingFreq, hh));
        
        figure(30);
        clf
        for nn = 1:length(nPulseCoordsImaging)
            xline(nPulseCoordsImaging(nn)/imagingFreq, '--', 'LineWidth', 1.5, 'color', [.8 .8 .8]);
        end
        hold on
        plot(linspace(0, length(holoSortedImagingMean{cc}(:, hh))/imagingFreq, length(holoSortedImagingMean{cc}(:, hh))), holoSortedImagingMean{cc}(:, hh), 'LineWidth', 1.5);
    %     axis off
        hold off
        xlabel('Time(s)')
        ylabel('df/f')
        pause
    end
end

% % Combine all hologram traces, calculate one grand mean
holoComboImagingAllTrials = [];
for cc = 1
    for hh = 1:nHolos(1)
        holoComboImagingAllTrials = [holoComboImagingAllTrials, holoSortedImagingAllTrials_ALT{cc}{hh}(:, :)];
    end
end
holoComboImagingGrandMean = nanmean(holoComboImagingAllTrials, 2);

CIDfAllConds = cell(nConds, 1);
for cc = 1:nConds 
    for hh = 1:nHolos(cc)
        confidence_level = 0.95;
        means = nanmean(holoSortedImagingAllTrials_ALT{cc}{hh, 1}, 2);
        std_errors = std(holoSortedImagingAllTrials_ALT{cc}{hh, 1}, 0, 2, "omitnan") / sqrt(size(holoSortedImagingAllTrials_ALT{cc}{hh, 1}, 2));
    
        t_score = tinv((1 + confidence_level) / 2, size(holoSortedImagingAllTrials_ALT{cc}{hh, 1}, 2) - 1);
        margin_of_error = t_score * std_errors;
        lower_bounds = means - margin_of_error;
        upper_bounds = means + margin_of_error;
        if UpOrDown == '2'
            CIDfAllConds{cc}{hh, 1} = [lower_bounds, upper_bounds];
        elseif UpOrDown =='1'
            CIDfAllConds{cc}{hh, 1} = [-lower_bounds, -upper_bounds];
        end
    end
end

voltMapping.holoSortedImagingAllTrials     = holoSortedImagingAllTrials;
voltMapping.holoSortedImagingAllTrials_ALT = holoSortedImagingAllTrials_ALT;
voltMapping.holoSeqIndex                   = holoSeqIndex;
voltMapping.holoSortedImagingMean          = holoSortedImagingMean;
voltMapping.holoComboImagingAllTrials      = holoComboImagingAllTrials;
voltMapping.holoComboImagingGrandMean      = holoComboImagingGrandMean;
voltMapping.CIDfAllConds                   = CIDfAllConds;
voltMapping.preStimWindow = preStimWindow;
voltMapping.postStimWindow = postStimWindow;

%% Align imaging with ephys traces
for cc = 1:nConds
    for hh = 1:nHolos(cc)
        figure(100);
%         set(gcf, 'Position',  [100, 100, 600, 400])
        clf
        hold on
        fill([linspace(0, size(CIDfAllConds{cc}{hh}, 1)/imagingFreq, size(CIDfAllConds{cc}{hh}, 1)), fliplr(linspace(0, size(CIDfAllConds{cc}{hh}, 1)/imagingFreq, size(CIDfAllConds{cc}{hh}, 1)))],...
        [CIDfAllConds{cc}{hh}(:, 1)', fliplr(CIDfAllConds{cc}{hh}(:, 2)')], [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
        % plot CI lowerbound
        plot(linspace(0, size(CIDfAllConds{cc}{hh}, 1)/imagingFreq, size(CIDfAllConds{cc}{hh}, 1)), CIDfAllConds{cc}{hh}(:, 1), '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]); 
        % plot CI upperbound
        plot(linspace(0, size(CIDfAllConds{cc}{hh}, 1)/imagingFreq, size(CIDfAllConds{cc}{hh}, 1)), CIDfAllConds{cc}{hh}(:, 2), '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);   
        
        % show line at dff = 0
        yline(0, '-', 'LineWidth', 1.5, 'color', [0.9 0.9 0.9]);
        for nn = 1:length(nPulseCoords)
            xline(nPulseCoords(nn)/Fs, '--', 'LineWidth', 1.5, 'color', [.8 .8 .8]);
        end

        % plot ephys and voltage traces
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

        % plot voltage trace only
        plot(linspace(0, size(holoSortedImagingMean{cc}(:, hh), 1)/imagingFreq, size(holoSortedImagingMean{cc}(:, hh), 1)), holoSortedImagingMean{cc}(:, hh), '-', 'linewidth', 2, 'color', 'g');
        ylabel('dF/F (%)');
        xlabel('Time (s)');
%         ylim([-1 max(holoSortedImagingMean{pp}(:, hh))])
%         ylim([min(holoSortedImagingMean{cc}(:, hh)) max(holoSortedImagingMean{cc}(:, hh))])
%         
%         % plot ephys trace only
%         yyaxis right
%     %     axes('Position',[.70 .12 .2 .2]);
%     %     box on
%     %     plot(linspace(0, size(meanTraceAllConds, 1)/Fs, size(meanTraceAllConds, 1)), meanTraceAllConds(:, cc), 'linewidth', 1, 'color', [0.3010 0.7450 0.9330]);
%         plot(linspace(0, size(holoSortedDataMean{cc}(:, hh), 1)/Fs, size(holoSortedDataMean{cc}(:, hh), 1)), holoSortedDataMean{cc}(:, hh), 'linewidth', 1, 'color', [0.3010 0.7450 0.9330]);
%         gca;
%         set(gca,'xtick',[], 'fontsize', 18);
    %     ylim([min(meanTraceAllConds(:, cc)), max(meanTraceAllConds(:, cc))]);
%         ylabel('dV');
        
        hold off
     pause
    end
end

%% Save Analysis Results
directory = 'D:\Data\Voltage Imaging\Mapping\Mapping Analysis Results';
expID = num2str(ExpStruct.mouseID);
voltMapping.mouseID = ['voltMapping_Analysis_', expID];
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