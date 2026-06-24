%% Voltage Imaging Poisson Test (on slice) Analysis Code
% This spinoff of the slice dvTest is for when the injected current is done
% in poisson fashion. All injected current should be depolarizing enough
% for spiking...
%%
clear all
%% Load files and setup
% Step 1: Read the ephys file
ephysFilePath = char(uigetdir('H:\Data\Voltage Imaging\VoltImg_slice_test\Ephys Data')); % Select and set root folder where all experiments with cells you want to analyze are located
% ephysFilePath = char(uigetdir('/Volumes/Elements/Data/Voltage Imaging/VoltImg_slice_test/Ephys Data')); % Select and set root folder where all experiments with cells you want to analyze are located
ephysFileDir = dir(ephysFilePath);
load([ephysFileDir(3).folder, '/', ephysFileDir(3).name]);

% Get dv parameters and condition sequence
dvCondSequence = ExpStruct.dvStepParams.dvCondSequence; % sequence of randomized dv trials
dvToTest = ExpStruct.dvStepParams.dvToTest; % dv steps to be simulated by current injection
nConds = length(unique(dvCondSequence)); % number of conditions (dv steps)

pulseStart = ExpStruct.dvStepParams.pulseStart;
sweepDur = ExpStruct.dvStepParams.sweepDur;
pulseFreq = ExpStruct.dvStepParams.pulseFreq;
imagingFreq = 330.31;
Fs = ExpStruct.Fs;

% Step 2: Identify the folder containing the imaging files correspdonding to the electrophysiology recording
ImgsFilePath = char(uigetdir('H:\Data\Voltage Imaging\VoltImg_slice_test\Imaging Data')); % Select and set root folder where all experiments with cells you want to analyze are located
% ImgsFilePath = char(uigetdir('/Volumes/Elements/Data/Voltage Imaging/VoltImg_slice_test/Imaging Data')); % Select and set root folder where all experiments with cells you want to analyze are located

ImgfolderContents = dir(ImgsFilePath);

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

% Step 2b: revise cond sequence to match actual trials recorded
dvCondSequence = dvCondSequence(1:(size(imagesIndex)));
% dvCondSequence = ones(1, 2000)

% Step 3: Identify bad trials, to skip in analyses below
vsTest_inputs = ExpStruct.dvStepParams.vsTest_inputs;
baselineAllTrials = [];
excludeTrials = [];
for tt = 1:size(vsTest_inputs, 2)
    baseline = mean(vsTest_inputs(1:0.001*pulseStart*Fs, tt));
    baselineAllTrials = [baselineAllTrials, baseline];
    
    if baselineAllTrials(tt) > -50
        excludeTrials = [excludeTrials, tt];
    end
end

% Step 3: Identify GEVI
UpOrDown = input('1 for upward GEVI, 2 for downward GEVI ', 's');

% Step 4: Calculate frameclock offsets
frameClock_inputs = [];
risingEdgeIndices = [];
risingEdgeTimes = [];
if isfield(ExpStruct.dvStepParams, 'frameClock_inputs')
    frameClock_inputs = ExpStruct.dvStepParams.frameClock_inputs;
    for tt = 1:size(frameClock_inputs, 2)
        trialTime = (0:length(frameClock_inputs(:, tt))-1)*(1/Fs);
        riseThreshold = 0; % rising edge threshold, any number well above baseline will do
        if any(frameClock_inputs(:, tt) == 1)
            risingEdgeIndices(tt) = find(frameClock_inputs(:, tt) > riseThreshold, 1); % find the index for first rising edge
            risingEdgeTimes(tt) = trialTime(risingEdgeIndices(tt));
        else
            risingEdgeIndices(tt) = 0;
            risingEdgeTimes(tt) = 0;
        end
    end
else
    risingEdgeTimes = repmat(0.0011, 1, size(vsTest_inputs, 2) - size(excludeTrials, 2));
end

voltImgTest_Poisson_Analysis.ephysData.Fs = Fs;
voltImgTest_Poisson_Analysis.pulseParams.pulseStart = pulseStart;
voltImgTest_Poisson_Analysis.pulseParams.sweepDur = sweepDur;
voltImgTest_Poisson_Analysis.pulseParams.pulseFreq = pulseFreq;
voltImgTest_Poisson_Analysis.imagingFreq = imagingFreq;
voltImgTest_Poisson_Analysis.dvCondSequence = dvCondSequence;
voltImgTest_Poisson_Analysis.dvToTest = dvToTest;
voltImgTest_Poisson_Analysis.Rinput = ExpStruct.dvStepParams.Rinput;
voltImgTest_Poisson_Analysis.fRinput = ExpStruct.dvStepParams.fRinput;
voltImgTest_Poisson_Analysis.ephysData.vsTest_inputs = vsTest_inputs;
voltImgTest_Poisson_Analysis.ephysData.baselineAllTrials = baselineAllTrials;
voltImgTest_Poisson_Analysis.ephysData.excludeTrials = excludeTrials;

%% Calculate ROI mask
% Step 1: Calculate mask for ROI using spiking resp
maxDvTrials = [];
maxDvTrials = find(dvCondSequence == max(unique(dvCondSequence)));

maxDvStack = [];
for tt = 1:length(maxDvTrials)
    if ismember(maxDvTrials(tt), excludeTrials)
        continue
    end
    currImgPath = [ImgfolderContents(imagesIndex(maxDvTrials(tt))).folder, '/', ImgfolderContents(imagesIndex(maxDvTrials(tt))).name];
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

figure(10); clf; imagesc(meanFluorMaxDvStack); axis equal; axis image; colorbar; % caxis([-7 -4]);

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

voltImgTest_Poisson_Analysis.maxDvTrials = maxDvTrials;
voltImgTest_Poisson_Analysis.meanFluorMaxDvStack = meanFluorMaxDvStack;
voltImgTest_Poisson_Analysis.roiMeanMaxDvStack = roiMeanMaxDvStack;
voltImgTest_Poisson_Analysis.roiX = roiX;
voltImgTest_Poisson_Analysis.roiY = roiY;

%%
dfAllConds = cell(nConds, 1);
traceAllConds = cell(nConds, 1);
f0Start =  0.01; % in seconds
f0End = ExpStruct.dvStepParams.pulseStart/1000;% 0.03; in seconds

% SpikingIndices = [];
% SpikingTimes = [];
%     for tt = 1:size(vsTest_inputs, 2)
%         trialTime = (0:length(vsTest_inputs(:, tt))-1)*(1/Fs);
%         riseThreshold = 0; % rising edge threshold, any number well above baseline will do
%         if any(frameClock_inputs(:, tt) == 1)
%             SpikingIndices{tt} = find(vsTest_inputs(:, tt) > riseThreshold, 1); % find the index for first rising edge
%             SpikingTimes(tt) = trialTime(SpikingIndices(tt));
%         else
%             SpikingIndices(tt) = 0;
%             SpikingTimes(tt) = 0;
%         end
%     end


for tt = 1:size(vsTest_inputs, 2)
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
        df = -df;
    end
    
    dfAllConds{dvCondSequence(tt), 1} = [dfAllConds{dvCondSequence(tt), 1}, df];
    traceAllConds{dvCondSequence(tt), 1} = [traceAllConds{dvCondSequence(tt), 1}, vsTest_inputs(:, tt)];
end


for tt = 1:10 %size(vsTest_inputs, 2)
    figure(30+tt)
    clf
    hold on;
    plot(linspace(0, size(traceAllConds{1}, 1)/Fs, size(traceAllConds{1}, 1)), traceAllConds{1}(:, tt), 'linewidth', 1, 'color', [0.8 0.8 0.8])
    plot(linspace(0, size(dfAllConds{1}, 1)/imagingFreq, size(dfAllConds{1}, 1)), dfAllConds{1}(:, tt), 'linewidth', 1, 'color', 'k')
    hold off
    ylim([-40 50])
    xlabel(['time(s)'])
    ylabel(['df/f'])
end

voltImgTest_Analysis.maxDvStack = maxDvStack;
voltImgTest_Analysis.meanFluorMaxDvStack = meanFluorMaxDvStack;
voltImgTest_Analysis.stdFluor = stdFluor;
voltImgTest_Analysis.meanFluor = meanFluor;
voltImgTest_Analysis.cutOffFluor = cutOffFluor; 
voltImgTest_Analysis.roiStack = roiStack;
voltImgTest_Analysis.dfAllConds = dfAllConds;
voltImgTest_Analysis.traceAllConds = traceAllConds;

%% Save Analysis Results
directory = 'D:\Data\Voltage Imaging\voltImg_slice_test\Analysis Results';
voltImgTest_Analysis.mouseID = ['voltImgTest_Analysis_', num2str(ExpStruct.mouseID)];
fileName = [num2str(voltImgTest_Analysis.mouseID), '.mat'];
save(fullfile(directory, fileName), 'voltImgTest_Analysis', '-v7.3');

