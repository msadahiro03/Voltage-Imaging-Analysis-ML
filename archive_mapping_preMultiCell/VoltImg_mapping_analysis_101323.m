%% Voltage Imaging Mapping (VIM) Analysis
% This analysis is used for 2P optogenetic mapping experiments where the readout is either ephys or imaging.

%%
clear all

%% Load files and setup
% Step 1: Read the ephys file
% ephysFilePath = char(uigetdir('/Volumes/Elements/Data/Voltage Imaging/Mapping/231016/SSTCre_tdTomato_Slice_101623_Cell1_HoloInputTest')); % Select and set root folder where all experiments with cells you want to analyze are located
% ephysFilePath = char(uigetdir('E:\Data\Voltage Imaging\Mapping\Ephys data')); % Select and set root folder where all experiments with cells you want to analyze are located
ephysFilePath = char(uigetdir('/Volumes/Elements/Data/Voltage Imaging/Mapping/Ephys data')); % Select and set root folder where all experiments with cells you want to analyze are located

ephysFileDir = dir(ephysFilePath);
load([ephysFileDir(3).folder, '/', ephysFileDir(3).name]);
disp(ephysFileDir(3).name);
 
% Step 2: Identify the folder containing the imaging files correspdonding to the ephys file
% ImgsFilePath = char(uigetdir('E:\Data\Voltage Imaging\Mapping\Imaging data')); % Select and set root folder where all experiments with cells you want to analyze are located
ImgsFilePath = char(uigetdir('/Volumes/Elements/Data/Voltage Imaging/Mapping/Imaging data')); % Select and set root folder where all experiments with cells you want to analyze are located

ImgfolderContents = dir(ImgsFilePath);

% Step 2a: Avoid hidden files and non image files in the image folder
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

imagingFreq = 330.31;

% Step 3: Setup struct
voltMapping = ExpStruct;

if exist('ExpStruct2', 'var')
    voltMapping.ExpStruct2 = ExpStruct2;
end

% Singlespot grid properties
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

% Step 4: Identify GEVI
UpOrDown = input('1 for upward GEVI, 2 for downward GEVI ', 's');

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
voltMapping.UpOrDown = UpOrDown;
voltMapping.ephysFilePath = ephysFilePath;
voltMapping.ImgsFilePath = ImgsFilePath;

%% Trial exclusion and baselining
mappingInputs = ExpStruct.inputs;

for tt = 1:length(ExpStruct.inputs)
    ExpStruct.inputs{tt} = ExpStruct.inputs{tt}*2;
end

baseVoltAllTrials = [];
excludeTrials = [];
mappingInputsBaselined =[];
for tt = 1:size(mappingInputs, 1)
    baseVolt = mean(mappingInputs{tt}(1:startTime*Fs));
    baseVoltAllTrials = [baseVoltAllTrials; baseVolt];
    mappingInputsBaselined{tt, 1} = mappingInputs{tt} - baseVolt;

    if baseVoltAllTrials(tt) > -45
        excludeTrials = [excludeTrials, tt];
    end
end

voltMapping.excludeTrials = excludeTrials;
voltMapping.ephys.baseVoltAllTrials = baseVoltAllTrials;
voltMapping.ephys.mappingInputsBaselined = mappingInputsBaselined;

%% Break apart each continuous ephys trace into stim windows and rearrange according to hologram sequence
holoSeqIndex = cell(nConds, 1);
holoSortedDataAllTrials = cell(nConds, 1);
    for tt = 1:nTrials
        holoSeqIndex{voltMapping.trialCond(tt, 1)} = [holoSeqIndex{voltMapping.trialCond(tt, 1)}, unique(voltMapping.outParams.sequenceThisTrial{tt}, 'stable')']; % Find hologram sequences for every trial

        for hh = 1:nHolos
            holoSortedDataAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, tt), 1}(:, tt) = mappingInputsBaselined{tt}(voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)*Fs:(voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)*Fs)+(ipi*nPulses/1000*Fs));
        end
    end

% Remove trials where dv goes below zero (not sure why some trials do that  but those end up messing up the mean)
for pp = 1:nConds
    for hh = 1:nHolos
        for tt = 1:nTrials
            if mean(holoSortedDataAllTrials{pp}{hh}(:, tt)) < 0.1
                holoSortedDataAllTrials{pp}{hh}(:, tt) = NaN;
            end
        end
    end
end

nPulseCoords = [];
for ii = 1:nPulses
    nPulseCoords = [nPulseCoords, ii * ipi/1000*Fs];
end

nPulseCoordsImaging = [];
for ii = 1:nPulses
    nPulseCoordsImaging = [nPulseCoordsImaging, ii * ipi/1000*imagingFreq];
end

holoSortedDataMean = cell(nConds, 1);
for pp = 1:nConds
    for hh = 1:nHolos
        holoSortedDataMean{pp}(:, hh) = nanmean(holoSortedDataAllTrials{pp}{hh}, 2);
        
        % Baseline the mean holo traces
        holoSortedDataMean{pp}(:, hh) = holoSortedDataMean{pp}(:, hh) - mean(holoSortedDataMean{pp}(1:50, hh));
        
        figure(1);
        clf
        for nn = 1:length(nPulseCoords)
            xline(nPulseCoords(nn)/Fs, '--', 'LineWidth', 1.5, 'color', [.8 .8 .8]);
        end
        hold on
        plot(linspace(0, length(holoSortedDataMean{pp}(:, hh))/Fs, length(holoSortedDataMean{pp}(:, hh))), filtfilt(blp,alp,holoSortedDataMean{pp}(:, hh)), 'LineWidth', 1.5);
    %     axis off
        hold off
        xlabel('Time(s)')
        ylabel('mV')
        pause
    end
end

voltMapping.ephys.holoSeqIndex = holoSeqIndex; 
voltMapping.ephys.holoSortedDataAllTrials = holoSortedDataAllTrials;
voltMapping.ephys.nPulseCoords = nPulseCoords;
voltMapping.ephys.holoSortedDataMean = holoSortedDataMean;

%% Calculate ROI mask
% Step 1: Calculate mask for ROI 
maxDvStack = [];
for tt = 1:20%length(imagesIndex)
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

figure(20); clf;
imagesc(roiStack); axis equal; axis image;
[roiX, roiY] = find(roiStack);

voltMapping.meanFluorMaxDvStack = meanFluorMaxDvStack;
voltMapping.roiMeanMaxDvStack = roiMeanMaxDvStack;
voltMapping.roiX = roiX;
voltMapping.roiY = roiY;

%% Calculate df and mean df for all conditions
% Step 1: Calculate df across all trials across all conditions

dfAllTrialsCondSorted = cell(nConds, 1);
dfAllTrialsNonSorted = cell(nTrials, 1);
traceAllConds = cell(nConds, 1);
f0Start =  0.01; % in seconds
f0End = startTime;

for tt = 1:length(imagesIndex) %size(vsTest_inputs, 2)
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

    dfAllTrialsCondSorted{voltMapping.trialCond(tt, 1), 1} = [dfAllTrialsCondSorted{voltMapping.trialCond(tt, 1), 1}, df];
    dfAllTrialsNonSorted{tt} = df;
end

% Break apart trial-by-trial imaging data into stim windows and rearrange according to hologram sequence
holoSortedImagingAllTrials = cell(nConds, 1);
for tt = 1:nTrials
    for hh = 1:nHolos
        dfThisHolo = dfAllTrialsNonSorted{tt}(voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)*imagingFreq:(voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)*imagingFreq)+(ipi*nPulses/1000*imagingFreq));
        dfThisHolo = dfThisHolo - mean(dfThisHolo(1:3));

        holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, tt), 1}(:, tt) = dfThisHolo;
    end
end

holoSortedImagingMean = cell(nConds, 1);
for pp = 1:nConds
    for hh = 1:nHolos
        holoSortedImagingMean{pp}(:, hh) = nanmean(holoSortedImagingAllTrials{pp}{hh}, 2);
        
        figure(30);
        clf
        for nn = 1:length(nPulseCoordsImaging)
            xline(nPulseCoordsImaging(nn)/imagingFreq, '--', 'LineWidth', 1.5, 'color', [.8 .8 .8]);
        end
        hold on
        plot(linspace(0, length(holoSortedImagingMean{pp}(:, hh))/imagingFreq, length(holoSortedImagingMean{pp}(:, hh))), holoSortedImagingMean{pp}(:, hh), 'LineWidth', 1.5);
    %     axis off
        hold off
        xlabel('Time(s)')
        ylabel('mV')
%         pause
    end
end

CIDfAllConds = cell(nConds, 1);
for pp = 1:nConds 
    for hh = 1:nHolos
        confidence_level = 0.95;
        means = mean(holoSortedImagingAllTrials{pp}{hh, 1}, 2);
        std_errors = std(holoSortedImagingAllTrials{pp}{hh, 1}, 0, 2) / sqrt(size(holoSortedImagingAllTrials{pp}{hh, 1}, 2));
    
        t_score = tinv((1 + confidence_level) / 2, size(dfAllConds{cc}, 2) - 1);
        margin_of_error = t_score * std_errors;
        lower_bounds = means - margin_of_error;
        upper_bounds = means + margin_of_error;
        if UpOrDown == '2'
            CIDfAllConds{pp}{hh, 1} = [-lower_bounds, -upper_bounds];
        elseif UpOrDown =='1'
            CIDfAllConds{pp}{hh, 1} = [lower_bounds, upper_bounds];
        end
    end
end


%% Align imaging with ephys traces
for pp = 1:nConds
    for hh = 1:nHolos
        figure(40+hh)
        clf
        % plot voltage trace
        plot(linspace(0, size(holoSortedImagingMean{pp}(:, hh), 1)/imagingFreq, size(holoSortedImagingMean{pp}(:, hh), 1)), holoSortedImagingMean{pp}(:, hh), '-', 'linewidth', 2, 'color', 'k');
        ylim([-1 2.5])
        ylabel('dF/F (%)');
        xlabel('Time (s)');
        
        % show line at dff = 0
        yline(0, '-', 'LineWidth', 1.5, 'color', [0.9 0.9 0.9]);
        hold off
            for nn = 1:length(nPulseCoords)
            xline(nPulseCoords(nn)/Fs, '--', 'LineWidth', 1.5, 'color', [.8 .8 .8]);
        end
        % plot ephys trace
        yyaxis right
    %     axes('Position',[.70 .12 .2 .2]);
    %     box on
    %     plot(linspace(0, size(meanTraceAllConds, 1)/Fs, size(meanTraceAllConds, 1)), meanTraceAllConds(:, cc), 'linewidth', 1, 'color', [0.3010 0.7450 0.9330]);
        plot(linspace(0, size(holoSortedDataMean{pp}(:, hh), 1)/Fs, size(holoSortedDataMean{pp}(:, hh), 1)), holoSortedDataMean{pp}(:, hh), 'linewidth', 1, 'color', [0.3010 0.7450 0.9330]);
        gca;
        set(gca,'xtick',[], 'fontsize', 18);
    %     ylim([min(meanTraceAllConds(:, cc)), max(meanTraceAllConds(:, cc))]);
        ylabel('dV');
     
    end
end

%% Save Analysis Results
directory = 'D:\Data\Voltage Imaging\Mapping\Mapping Analysis Results';
voltMapping.mouseID = ['voltMapping_Analysis_', num2str(ExpStruct.mouseID)];
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