%% Voltage Imaging Mapping (VIM) Analysis
% This analysis is used for 2P optogenetic mapping experiments where the readout is either ephys or imaging.

%%
clear all
close all

%% Load files and setup
% Step 1: Read the ephys file
% ephysFilePath = char(uigetdir('/Volumes/ExData2/Voltage Imaging/VisStim_GroundTruth/Ephys Data'));
ephysFilePath = char(uigetdir('E:\Voltage Imaging\VisStim_GroundTruth\Ephys Data')); % Select and set root folder where all experiments with cells you want to analyze are located
% ephysFilePath = char(uigetdir('E:\Data\Voltage Imaging\dvTest\Ephys Data')); 
% ephysFilePath = char(uigetdir('/Volumes/Untitled/Voltage Imaging/dvTest/Ephys')); % Select and set root folder where all experiments with cells you want to analyze are located
ephysFileDir = dir(ephysFilePath);
load([ephysFileDir(3).folder, '/', ephysFileDir(3).name]);

% Step 2: Identify the folder containing the imaging files
% ImgsFilePath = char(uigetdir('E:\Voltage Imaging\voltMapping\Imaging')); % Select and set root folder where all experiments with cells you want to analyze are located
ImgsFilePath = char(uigetdir('E:\Voltage Imaging\VisStim_GroundTruth\Imaging Data')); % Select and set root folder where all experiments with cells you want to analyze are located
% ImgsFilePath = char(uigetdir('/Volumes/ExData2/Voltage Imaging/VisStim_GroundTruth/Imaging Data'));  % Select and set root folder where all experiments with cells you want to analyze are located
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

imagingFreq = 330.22;
% 330.33, 512x32
% 165.11, 512x64
% 149.56, 512x90
% 299.12, 512x37

% Step 3: Load visual stimulus sequence file
visStimSeqFilePath = char(uigetdir('E:\Voltage Imaging\VisStim_GroundTruth\visStimSequences'));
% visStimSeqFilePath = char(uigetdir('/Volumes/ExData2/Voltage Imaging/VisStim_GroundTruth/visStimSequences'));
stimSeqFileName = uigetfile(visStimSeqFilePath);
load([visStimSeqFilePath,'/', stimSeqFileName]);

% Step 4: Setup
nStims = length(unique(stimulus_sequence));
nTrials = length(imagesIndex); % ALTERNATIVELY "length(find(cellfun(@isempty, ExpStruct.inputs)==0))". Instead of "length(ExpStruct.inputs)" this puts out true number of trials successfully recorded
preStimWindow = 250; % Pre-stim time = 0.1 or 0.5s

% Step 5: Identify GEVI type
UpOrDown = input('1 for upward GEVI, 2 for downward GEVI: ', 's');

% Step 6: Save directory
directory = 'D:\Data\Voltage Imaging\8visStim\Analysis Results\';
fileName = ['voltMapping ', '.mat'];

voltImg_8VisStim.imagesIndex = imagesIndex;
voltImg_8VisStim.imagingFreq = imagingFreq;
voltImg_8VisStim.UpOrDown = UpOrDown;
voltImg_8VisStim.ImgsFilePath = ImgsFilePath;

%% Removal of trials where no spikes occurred 
positiveTrials = [];
for tt = 1:nTrials
    [~, locs] = findpeaks(-ExpStruct.inputs(:, tt), ExpStruct.Fs,'MinPeakHeight', 5);
    if ~isempty(locs)
        positiveTrials = [positiveTrials, tt];
    end
end

%% Calculate ROI mask
% Step 1: Calculate mask for ROI 
randTrialsForMask = 1:nTrials; % randperm(length(imagesIndex), 400);

% Preallocate maxDvStack
currImgPath = [ImgfolderContents(imagesIndex(1)).folder, '/', ImgfolderContents(imagesIndex(1)).name];
info = imfinfo(currImgPath);
maxDvStack = zeros(info(1).Height, info(1).Width, 400);

counter = 0;
for tt = randTrialsForMask
    counter = counter+1;
    disp(['Trial number: ', num2str(counter), ' of ', num2str(nTrials) ]);

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

meanMaxDvStack = mean(maxDvStack(:, :, :), 3);
meanFluorMaxDvStack = meanMaxDvStack;

greenScale = [zeros(256,1), (0:255)'/255, zeros(256,1)];
figure(10); clf; imagesc(meanFluorMaxDvStack); colormap(greenScale); axis equal; axis image; axis off; colorbar off; %colorbar; set(gca, 'fontsize', 12); 

% figure(10); clf; colormap('winter'); imagesc(meanFluorMaxDvStack); axis equal; axis image; colorbar; set(gca, 'fontsize', 12);% caxis([-7 -4]);

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

voltImg_8VisStim.meanFluorMaxDvStack = meanFluorMaxDvStack;
voltImg_8VisStim.roiMeanMaxDvStack = roiMeanMaxDvStack;
voltImg_8VisStim.cutOffFluor = cutOffFluor;
voltImg_8VisStim.roiX = roiX;
voltImg_8VisStim.roiY = roiY;

%% Calculate df and mean df for trials
EightVisStimImagingAllTrials = [];
filtEightVisStimImagingAllTrials = [];

counter = 0;
for tt = 1:nTrials
    counter = counter+1;
    disp(['Trial number: ', num2str(counter), ' of ', num2str(nTrials) ]);
    % Read the multi-frame image
    currImgPath = [ImgfolderContents(imagesIndex(tt)).folder, '/', ImgfolderContents(imagesIndex(tt)).name];
    info = imfinfo(currImgPath);
    numFrames = numel(info);

    % Preallocate the image stack
    imageStack = zeros(info(1).Height, info(1).Width, numFrames);
    
    % Read each frame and store in the stack
    for frameIndex = 1:numFrames
        imageStack(:,:,frameIndex) = imread(currImgPath, 'Index', frameIndex, 'Info', info);

        % Run these lines if artifact removal is necessary, comment out if not needed
        % [cleanFrame] = VoltImg_mapping_removeArtifact(imageStack(:,:,frameIndex));
        % imageStack(:,:,frameIndex) = cleanFrame; %replace the raw frame with cleaned up frame where artifact-corrupt lines are NaN'd
    end
    
    % Set baseline(f0), this is based on a single mean across the preStimWindow, so it may be a flawed approach!
    roiBaselineAllPixels = imageStack(roiX, roiY, floor(imagingFreq*0.1):ceil(preStimWindow/1000*imagingFreq)); % holistically set to sample from frame #50 instead of 1
    roiBaselineMean = nanmean(roiBaselineAllPixels, 'all');
    
    dfThisTrial = [];
    dffThisTrial = [];
    for ff = 1:size(imageStack, 3)
        % ROI pixels for this frame
        currFrameRoi = imageStack(roiX, roiY, ff);
        currFrameRoiMean = nanmean(currFrameRoi, 'all');
        
        % Calculate df
        intensityChange = currFrameRoiMean - roiBaselineMean;
        dfThisTrial = [dfThisTrial; intensityChange];
        dffThisTrial = [dffThisTrial; intensityChange/roiBaselineMean];
    end
    
    if UpOrDown == '2'
        dfThisTrial = -dfThisTrial;
        dffThisTrial = -dffThisTrial;
    elseif UpOrDown =='1'
        dfThisTrial = dfThisTrial;
        dffThisTrial = dffThisTrial;
    end
        
    % Use these lines if I need to further filter the df or dff trace
    % set filter for imaging
    cutOffFreq = 60;   % Cutoff frequency
    [b, a] = butter(4, cutOffFreq/(imagingFreq/2));  % 4th order Butterworth filter
    filtdffThisTrial = filter(b, a, dffThisTrial);
    
    % Or this filter
    %     f0 = 105;     % Frequency to remove (Hz)
    %     Q = 20;      % Quality factor (higher Q results in a narrower notch)
    %     % Design a notch filter (band-stop)
    %     w0 = f0/(imagingFreq/2);  % Normalize the frequency to the Nyquist frequency
    %     [b, a] = iirnotch(w0, w0/Q);
    %     filtdffThisTrial = filter(b, a, dffThisTrial);

    % Re-normalize the dff data
    dffThisTrial = dffThisTrial - nanmean(dffThisTrial(floor(imagingFreq*0.1):ceil(preStimWindow/1000*imagingFreq)));
    filtdffThisTrial = filtdffThisTrial - nanmean(filtdffThisTrial(floor(imagingFreq*0.1):ceil(preStimWindow/1000*imagingFreq)));

    % Use this condition if non-vis-responsive trials need to be counted out
    if max(dffThisTrial(75:end)*100)<8
        dffThisTrial(:) = NaN;
    end
  
    if max(filtdffThisTrial(75:end)*100)<8
        filtdffThisTrial(:) = NaN;
    end
    
    EightVisStimImagingAllTrials = [EightVisStimImagingAllTrials, dffThisTrial];
    filtEightVisStimImagingAllTrials = [filtEightVisStimImagingAllTrials, filtdffThisTrial];
  
end

[numRows, numCols] = size(filtEightVisStimImagingAllTrials);
filtVisStimImagingSorted_8Vis = cell(nStims, 1);
VisStimImagingSorted_8Vis = cell(nStims, 1);
VisStimEphysSorted_8Vis = cell(nStims, 1);
% Loop through the columns and place them in the correct group
for tt = 1:numCols
    filtVisStimImagingSorted_8Vis{stimulus_sequence(tt)} = [filtVisStimImagingSorted_8Vis{stimulus_sequence(tt)} filtEightVisStimImagingAllTrials(:, tt)];
    VisStimImagingSorted_8Vis{stimulus_sequence(tt)} = [VisStimImagingSorted_8Vis{stimulus_sequence(tt)} EightVisStimImagingAllTrials(:, tt)];
    VisStimEphysSorted_8Vis{stimulus_sequence(tt)} = [VisStimEphysSorted_8Vis{stimulus_sequence(tt)} ExpStruct.inputs(:, tt)];
end
    
% Now merge the 8 groups into 4 groups (1+5, 2+6, 3+7, 4+8), combining 2 directions into 1
numFinalGroups = 4;
VisStimImagingSorted_4Vis = cell(numFinalGroups, 1);
filtVisStimImagingSorted_4Vis = cell(numFinalGroups, 1);
VisStimEphysSorted_4Vis = cell(numFinalGroups, 1);
for vv = 1:numFinalGroups
    % Combine groups i and i+4
    VisStimImagingSorted_4Vis{vv} = [VisStimImagingSorted_4Vis{vv}, VisStimImagingSorted_8Vis{vv}, VisStimImagingSorted_8Vis{vv+4}];
    filtVisStimImagingSorted_4Vis{vv} = [filtVisStimImagingSorted_4Vis{vv}, filtVisStimImagingSorted_8Vis{vv}, filtVisStimImagingSorted_8Vis{vv+4}];
    VisStimEphysSorted_4Vis{vv} = [VisStimEphysSorted_4Vis{vv}, VisStimEphysSorted_8Vis{vv}, VisStimEphysSorted_8Vis{vv+4}];
end
    
meanfiltVisStim_8Vis = [];
meanVisStim_8Vis = [];
meanVisStim_Ephys_8Vis = [];
for vv = 1:length(filtVisStimImagingSorted_8Vis)
    meanfiltVisStim_8Vis(:, vv) = nanmean(filtVisStimImagingSorted_8Vis{vv}, 2);
    meanVisStim_8Vis(:, vv) = nanmean(VisStimImagingSorted_8Vis{vv}, 2);
    meanVisStim_Ephys_8Vis(:, vv) = nanmean(VisStimEphysSorted_8Vis{vv}, 2);
end

meanfiltVisStim_4Vis = [];
meanVisStim_4Vis = [];
meanVisStim_Ephys_4Vis = [];
for vv = 1:length(filtVisStimImagingSorted_4Vis)
    meanfiltVisStim_4Vis(:, vv) = nanmean(filtVisStimImagingSorted_4Vis{vv}, 2);
    meanVisStim_4Vis(:, vv) = nanmean(VisStimImagingSorted_4Vis{vv}, 2);
    meanVisStim_Ephys_4Vis(:, vv) = nanmean(VisStimEphysSorted_4Vis{vv}, 2);
end
    
figure(99);
clf
plot(linspace(0, length(meanfiltVisStim_4Vis(:, 1))/imagingFreq, length(meanfiltVisStim_4Vis(:, 1))), meanfiltVisStim_4Vis(:, 1), 'LineWidth', 1.5);
hold on;
plot(linspace(0, length(meanfiltVisStim_4Vis(:, 2))/imagingFreq, length(meanfiltVisStim_4Vis(:, 2))), meanfiltVisStim_4Vis(:, 2), 'LineWidth', 1.5);
plot(linspace(0, length(meanfiltVisStim_4Vis(:, 3))/imagingFreq, length(meanfiltVisStim_4Vis(:, 3))), meanfiltVisStim_4Vis(:, 3), 'LineWidth', 1.5);
plot(linspace(0, length(meanfiltVisStim_4Vis(:, 4))/imagingFreq, length(meanfiltVisStim_4Vis(:, 4))), meanfiltVisStim_4Vis(:, 4), 'LineWidth' , 1.5);
xlabel('Time (s)');
ylabel('dF/F');
legend('0', '45', '90', '135');
xlim([0.25, 1]);
xline(0.3, ':', 'Linewidth', 3, 'color', [0.5 0.5 0.5]);
hold off

VisStimSpikeNum_4Vis = cell(numFinalGroups, 1);
meanVisStimSpikeNum_4Vis = [];
for vv = 1:length(VisStimEphysSorted_4Vis)
    for tt = 1:size(VisStimEphysSorted_4Vis{vv}, 2)
    [~, locs] = findpeaks(-VisStimEphysSorted_4Vis{vv}(:, tt), ExpStruct.Fs,'MinPeakHeight',6);
    if isempty(locs)
        VisStimSpikeNum_4Vis{vv}(tt) = 0;
    else
        VisStimSpikeNum_4Vis{vv}(tt) = length(locs);
    end
    end
    meanVisStimSpikeNum_4Vis(vv, 1) = mean(VisStimSpikeNum_4Vis{vv});
end


%%
close all;
for vv = 1:length(VisStimImagingSorted_4Vis)
    randTrials = randi([1, size(VisStimImagingSorted_4Vis{vv}, 2)], 1, 10);
    for tt = randTrials %1:size(VisStimImagingSorted_4Vis{vv}, 2) 
    figure(100*vv+tt);
    set(gcf,'Position',[100 100 750 750])
    clf

    yyaxis right
    plot(linspace(0, size(VisStimEphysSorted_4Vis{vv}(:, tt), 1)/ExpStruct.Fs, size(VisStimEphysSorted_4Vis{vv}(:, tt), 1)), -VisStimEphysSorted_4Vis{vv}(:, tt), '-', 'linewidth', 1.5, 'color', [0 0 0]);
    ylim([min(-VisStimEphysSorted_4Vis{vv}(:, tt))+1, max(-VisStimEphysSorted_4Vis{vv}(:, tt))]);
    ylabel('pA')
    hold on
%     yyaxis right
    plot(linspace(0, size(VisStimImagingSorted_4Vis{vv}(:, tt), 1)/imagingFreq, size(VisStimImagingSorted_4Vis{vv}(:, tt), 1)), VisStimImagingSorted_4Vis{vv}(:, tt)*100, '-', 'linewidth', 1.5, 'color', [0 0.75 0]);
    ylabel('dF/F (%)');
    xlabel('Time (s)');
    ylim([min(VisStimImagingSorted_4Vis{vv}(:, tt))*100 - 2 , max(VisStimImagingSorted_4Vis{vv}(:, tt))*100 + 0.1*(max(VisStimImagingSorted_4Vis{vv}(:, tt))*100)]);
%     hold off;

    set(gca, 'fontsize', 20);
%     plot([0; 0.2], [-10; -10], '-k', 'LineWidth', 2);
%     plot([0; 0], [-10; -5], '-k', 'LineWidth', 2);
    % axis off
    xline(0.3, ':', 'Linewidth', 3, 'color', [0.5 0.5 0.5]);
    end
end

%%
CIDfAllConds = cell(nConds, 1);
filtCIDfAllConds = cell(nConds, 1);
for cc = 1:nConds 
    for hh = 1:nHolos(cc)
        confidence_level = 0.95;
        means = nanmean(holoSortedImagingAllTrials_ALT{cc}{hh, 1}, 2);
        filtMeans = nanmean(filtHoloSortedImagingAllTrials_ALT{cc}{hh, 1}, 2);       
        std_errors = std(holoSortedImagingAllTrials_ALT{cc}{hh, 1}, 0, 2, "omitnan") / sqrt(size(filtHoloSortedImagingAllTrials_ALT{cc}{hh, 1}, 2));
        filtStd_errors = std(filtHoloSortedImagingAllTrials_ALT{cc}{hh, 1}, 0, 2, "omitnan") / sqrt(size(filtHoloSortedImagingAllTrials_ALT{cc}{hh, 1}, 2));
  
        t_score = tinv((1 + confidence_level) / 2, size(holoSortedImagingAllTrials_ALT{cc}{hh, 1}, 2) - 1);
        filtT_score = tinv((1 + confidence_level) / 2, size(filtHoloSortedImagingAllTrials_ALT{cc}{hh, 1}, 2) - 1);
        margin_of_error = t_score * std_errors;
        filtMargin_of_error = filtT_score * filtStd_errors;
        lower_bounds = means - margin_of_error;
        filtLower_bounds = filtMeans - filtMargin_of_error;
        upper_bounds = means + margin_of_error;
        filtUpper_bounds = filtMeans + filtMargin_of_error;
        if UpOrDown == '2'
            CIDfAllConds{cc}{hh, 1} = [lower_bounds, upper_bounds];
            filtCIDfAllConds{cc}{hh, 1} = [filtLower_bounds, filtUpper_bounds];
        elseif UpOrDown =='1'
            CIDfAllConds{cc}{hh, 1} = [-lower_bounds, -upper_bounds];
            filtCIDfAllConds{cc}{hh, 1} = [-filtLower_bounds, -filtUpper_bounds];

        end
    end
end

voltImg_8VisStim.holoSortedImagingAllTrials     = holoSortedImagingAllTrials;
voltImg_8VisStim.filtHoloSortedImagingAllTrials = filtHoloSortedImagingAllTrials;
voltImg_8VisStim.holoSortedImagingAllTrials_ALT = holoSortedImagingAllTrials_ALT;
voltImg_8VisStim.filtHoloSortedImagingAllTrials_ALT = filtHoloSortedImagingAllTrials_ALT;
voltImg_8VisStim.holoSeqIndex                   = holoSeqIndex;
voltImg_8VisStim.holoSortedImagingMean          = holoSortedImagingMean;
voltImg_8VisStim.filtHoloSortedImagingMean = filtHoloSortedImagingMean;
% voltMapping.holoComboImagingAllTrials      = holoComboImagingAllTrials;
% voltMapping.holoComboImagingGrandMean      = holoComboImagingGrandMean;
voltImg_8VisStim.CIDfAllConds                   = CIDfAllConds;
voltImg_8VisStim.filtCIDfAllConds              = filtCIDfAllConds ;
voltImg_8VisStim.preStimWindow = preStimWindow;
voltImg_8VisStim.postStimWindow = postStimWindow;

%% Save Analysis Results
directory = '/Volumes/Untitled/Voltage Imaging/voltMapping/Analysis Results';
expID = num2str(ExpStruct.mouseID);
voltImg_8VisStim.mouseID = ['voltMapping_Analysis_', expID];
fileName = [num2str(voltImg_8VisStim.mouseID), '.mat'];
save(fullfile(directory, fileName), 'voltMapping', '-v7.3');

TimeNow = datetime('now', 'Format', 'dd-MMM-yyyy HH:mm:ss');
disp(['finished saving at: ' char(TimeNow)])
