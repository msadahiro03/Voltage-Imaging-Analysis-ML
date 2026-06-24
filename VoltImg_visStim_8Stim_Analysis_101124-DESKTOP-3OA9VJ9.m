%% Voltage Imaging Mapping (VIM) Analysis
% This analysis is used for 2P optogenetic mapping experiments where the readout is either ephys or imaging.

%%
clear all
close all

%% Load files and setup
% Step 1: Identify the folder containing the imaging files
% ImgsFilePath = char(uigetdir('E:\Voltage Imaging\voltMapping\Imaging')); % Select and set root folder where all experiments with cells you want to analyze are located
% ImgsFilePath = char(uigetdir('F:\Voltage Imaging\In Vivo Visual Stim')); % Select and set root folder where all experiments with cells you want to analyze are located
ImgsFilePath = char(uigetdir('/Volumes/ExData2/Voltage Imaging/VisStim_GroundTruth/Imaging Data'));  % Select and set root folder where all experiments with cells you want to analyze are located
ImgfolderContents = dir(ImgsFilePath);
disp(ImgfolderContents(end).name);

% Step 1a: Avoid hidden files and non image files
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

imagingFreq = 330.33;
% 330.33, 512x32
% 165.11, 512x64
% 149.56, 512x90
% 299.12, 512x37

% Step1b: Load visual stimulus sequence file
% visStimSeqFilePath = char(uigetdir('F:\Voltage Imaging\visStimSequences'));
visStimSeqFilePath = char(uigetdir('/Volumes/ExData2/Voltage Imaging/VisStim_GroundTruth/visStimSequences'));
stimSeqFileName = uigetfile(visStimSeqFilePath);
load([visStimSeqFilePath,'/', stimSeqFileName]);

% Step 2: Setup
nStims = length(unique(stimulus_sequence));
nTrials = length(imagesIndex); % ALTERNATIVELY "length(find(cellfun(@isempty, ExpStruct.inputs)==0))". Instead of "length(ExpStruct.inputs)" this puts out true number of trials successfully recorded
preStimWindow = 500; % Pre-stim time = 0.1 or 0.5s

% Step 3: Identify GEVI type
UpOrDown = input('1 for upward GEVI, 2 for downward GEVI: ', 's');

% Step 4: Save directory
directory = 'D:\Data\Voltage Imaging\8visStim\Analysis Results\';
fileName = ['voltMapping ', '.mat'];

voltImg_8VisStim.imagesIndex = imagesIndex;
voltImg_8VisStim.imagingFreq = imagingFreq;
voltImg_8VisStim.UpOrDown = UpOrDown;
voltImg_8VisStim.ImgsFilePath = ImgsFilePath;

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
    roiBaselineAllPixels = imageStack(roiX, roiY, 50:ceil(preStimWindow/1000*imagingFreq)); % holistically set to sample from frame #50 instead of 1
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
    dffThisTrial = dffThisTrial - nanmean(dffThisTrial(75:150));
    filtdffThisTrial = filtdffThisTrial - nanmean(filtdffThisTrial(75:150));

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
% Loop through the columns and place them in the correct group
for tt = 1:numCols
    filtVisStimImagingSorted_8Vis{stimulus_sequence(tt)} = [filtVisStimImagingSorted_8Vis{stimulus_sequence(tt)} filtEightVisStimImagingAllTrials(:, tt)];
    VisStimImagingSorted_8Vis{stimulus_sequence(tt)} = [VisStimImagingSorted_8Vis{stimulus_sequence(tt)} EightVisStimImagingAllTrials(:, tt)];
end
    
% Now merge the 8 groups into 4 groups (1+5, 2+6, 3+7, 4+8), combining 2 directions into 1
numFinalGroups = 4;
VisStimImagingSorted_4Vis = cell(numFinalGroups, 1);
filtVisStimImagingSorted_4Vis = cell(numFinalGroups, 1);
for i = 1:numFinalGroups
    % Combine groups i and i+4
    VisStimImagingSorted_4Vis{i} = [VisStimImagingSorted_8Vis{i}, VisStimImagingSorted_8Vis{i+4}];
    filtVisStimImagingSorted_4Vis{i} = [filtVisStimImagingSorted_8Vis{i}, filtVisStimImagingSorted_8Vis{i+4}];
end
    
for vv = 1:length(filtVisStimImagingSorted_8Vis)
    meanfiltVisStim_8Vis(:, vv) = nanmean(filtVisStimImagingSorted_8Vis{vv}, 2);
    meanVisStim_8Vis(:, vv) = nanmean(VisStimImagingSorted_8Vis{vv}, 2);
end

for vv = 1:length(filtVisStimImagingSorted_4Vis)
    meanfiltVisStim_4Vis(:, vv) = nanmean(filtVisStimImagingSorted_4Vis{vv}, 2);
    meanVisStim_4Vis(:, vv) = nanmean(VisStimImagingSorted_4Vis{vv}, 2);
end
    
    % Combine the 2 directions per orientation
%     for vv = 1:
    
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
hold off

%%
% for vv = 1:length(VisStimImagingSorted_4Vis)
%     randTrials = randi([1, size(VisStimImagingSorted_4Vis{vv}, 2)], 1, 10);
%     for tt = randTrials %1:size(VisStimImagingSorted_4Vis{vv}, 2) 
%     figure(100*vv+tt);
%     set(gcf,'Position',[100 100 1000 1000])
%     clf
%     plot(linspace(0, size(VisStimImagingSorted_4Vis{vv}(:, tt), 1)/imagingFreq, size(VisStimImagingSorted_4Vis{vv}(:, tt), 1)), VisStimImagingSorted_4Vis{vv}(:, tt)*100, '-', 'linewidth', 1.5, 'color', [0 0.75 0]);
%     ylabel('dF/F (%)');
%     xlabel('Time (s)');
%     ylim([-15, 50])
%     hold on;
% 
%     plot([0; 0.2], [-10; -10], '-k', 'LineWidth', 2);
%     plot([0; 0], [-10; -5], '-k', 'LineWidth', 2);
% 
%     % axis off
% %     xline(1.5, 'Linewidth', 3, 'color', [0.5 0.5 0.5]);
%     end
% end
% hold off;

    figure(3333);
    set(gcf,'Position',[100 100 1000 500])
    clf
    plot(linspace(0, size(VisStimImagingSorted_4Vis{1, 1}(:, 38)  , 1)/imagingFreq, size(VisStimImagingSorted_4Vis{1, 1}(:, 38)   , 1)), VisStimImagingSorted_4Vis{1, 1}(:, 38)*100, '-', 'linewidth', 1.5, 'color', [0 0.75 0]);
    ylabel('dF/F (%)');
    xlabel('Time (s)');
    ylim([-25, 30])
    hold on;
    
    figure(3334);
    set(gcf,'Position',[100 100 1000 500])
    clf
    plot(linspace(0, size(VisStimImagingSorted_4Vis{1, 1}(:, 34)  , 1)/imagingFreq, size(VisStimImagingSorted_4Vis{1, 1}(:, 34)   , 1)), VisStimImagingSorted_4Vis{1, 1}(:, 34)*100, '-', 'linewidth', 1.5, 'color', [0 0.75 0]);
    ylabel('dF/F (%)');
    xlabel('Time (s)');
    ylim([-25, 30])
    hold on;
    
    figure(3335);
    set(gcf,'Position',[100 100 1000 500])
    clf
    plot(linspace(0, size(VisStimImagingSorted_4Vis{2, 1}(:, 33)  , 1)/imagingFreq, size(VisStimImagingSorted_4Vis{2, 1}(:, 33)   , 1)), VisStimImagingSorted_4Vis{2, 1}(:, 33)*100, '-', 'linewidth', 1.5, 'color', [0 0.75 0]);
    ylabel('dF/F (%)');
    xlabel('Time (s)');
    ylim([-25, 30])
    hold on;

    figure(3336);
    set(gcf,'Position',[100 100 1000 500])
    clf
    plot(linspace(0, size(VisStimImagingSorted_4Vis{2, 1}(:, 19)  , 1)/imagingFreq, size(VisStimImagingSorted_4Vis{2, 1}(:, 19)   , 1)), VisStimImagingSorted_4Vis{2, 1}(:, 19)*100, '-', 'linewidth', 1.5, 'color', [0 0.75 0]);
    ylabel('dF/F (%)');
    xlabel('Time (s)');
    ylim([-25, 30])
    hold on;

    plot([0; 0.2], [-10; -10], '-k', 'LineWidth', 2);
    plot([0; 0], [-10; -5], '-k', 'LineWidth', 2);
hold off;

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

%% Align imaging with ephys traces
for cc = 1:nConds
    for hh = 1:nHolos(cc)
        figure(100+10*cc+hh);
%         set(gcf, 'Position',  [100, 100, 600, 400])
        clf
        hold on
%         fill([linspace(0, size(CIDfAllConds{cc}{hh}, 1)/imagingFreq, size(CIDfAllConds{cc}{hh}, 1)), fliplr(linspace(0, size(CIDfAllConds{cc}{hh}, 1)/imagingFreq, size(CIDfAllConds{cc}{hh}, 1)))],...
%         [CIDfAllConds{cc}{hh}(:, 1)', fliplr(CIDfAllConds{cc}{hh}(:, 2)')], [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
%         % plot CI lowerbound
%         plot(linspace(0, size(CIDfAllConds{cc}{hh}, 1)/imagingFreq, size(CIDfAllConds{cc}{hh}, 1)), CIDfAllConds{cc}{hh}(:, 1), '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]); 
%         % plot CI upperbound
%         plot(linspace(0, size(CIDfAllConds{cc}{hh}, 1)/imagingFreq, size(CIDfAllConds{cc}{hh}, 1)), CIDfAllConds{cc}{hh}(:, 2), '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);   

        fill([linspace(0, size(filtCIDfAllConds{cc}{hh}, 1)/imagingFreq, size(filtCIDfAllConds{cc}{hh}, 1)), fliplr(linspace(0, size(filtCIDfAllConds{cc}{hh}, 1)/imagingFreq, size(filtCIDfAllConds{cc}{hh}, 1)))],...
        [filtCIDfAllConds{cc}{hh}(:, 1)'*100, fliplr(filtCIDfAllConds{cc}{hh}(:, 2)'*100)], [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
        % plot CI lowerbound
        plot(linspace(0, size(filtCIDfAllConds{cc}{hh}, 1)/imagingFreq, size(filtCIDfAllConds{cc}{hh}, 1)), filtCIDfAllConds{cc}{hh}(:, 1)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]); 
        % plot CI upperbound
        plot(linspace(0, size(filtCIDfAllConds{cc}{hh}, 1)/imagingFreq, size(filtCIDfAllConds{cc}{hh}, 1)), filtCIDfAllConds{cc}{hh}(:, 2)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);   
        
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
%         plot(linspace(0, size(holoSortedImagingMean{cc}{hh}, 1)/imagingFreq, size(holoSortedImagingMean{cc}{hh}, 1)), holoSortedImagingMean{cc}{hh}, '-', 'linewidth', 2, 'color', 'g');
        plot(linspace(0, size(holoSortedImagingMean{cc}{hh}, 1)/imagingFreq, size(holoSortedImagingMean{cc}{hh}, 1)), filter(b, a, holoSortedImagingMean{cc}{hh})*100, '-', 'linewidth', 2, 'color', 'g');
        ylabel('dF/F (%)');
        xlabel('Time (s)');
%         ylim([-1 max(holoSortedImagingMean{cc}{hh})*100])
%         ylim([min(holoSortedImagingMean{cc}{hh})*100 max(holoSortedImagingMean{cc}{hh})*100])
        ylim([-0.3 3.25])
        ax.YColor = [0 1 0];
        % axis off
%         
        % plot ephys trace only
        yyaxis right
    %     axes('Position',[.70 .12 .2 .2]);
    %     box on
        plot(linspace(0, size(holoSortedDataMean{cc}(:, hh), 1)/Fs, size(holoSortedDataMean{cc}(:, hh), 1)), holoSortedDataMean{cc}(:, hh), 'linewidth', 1, 'color', [0 0 0]);
        gca;
        set(gca,'xtick',[], 'fontsize', 18);
%         ylim([min(holoSortedDataMean{cc}(:, hh)), max(holoSortedDataMean{cc}(:, hh))]);
        ylim([-0.3 3.25]);
        ylabel('dV');
        ax.YColor = [0, 0, 0];
        % axis off
        
        % show line at dff = 0
%         yline(0, '-', 'LineWidth', 1.5, 'color', [0.9 0.9 0.9]);
        for nn = 1:length(nPulseCoords)
            xline(nPulseCoords(nn)/Fs, '-', 'LineWidth', 10, 'color', [1 0 0]);
        end
        
        plot([0.01; 0.035], [0.5; 0.5], '-k', 'LineWidth', 2);
        plot([0.01; 0.01], [0.5; 1], '-k', 'LineWidth', 2);
        hold off
     pause
    end
end

%% Save Analysis Results
directory = '/Volumes/ExData2/Voltage Imaging/VisStim_GroundTruth/Analysis Results';
expID = num2str(ExpStruct.mouseID);
voltImg_8VisStim.mouseID = ['voltVis_Analysis_', expID];
fileName = [num2str(voltImg_8VisStim.mouseID), '.mat'];
save(fullfile(directory, fileName), 'voltMapping', '-v7.3');

TimeNow = datetime('now', 'Format', 'dd-MMM-yyyy HH:mm:ss');
disp(['finished saving at: ' char(TimeNow)])

%% Load Analysis Results
% Run this section after loading the specific cell analysis file, then
% re-run the above sections to regenerate figures

names = fieldnames(voltImg_8VisStim);
for i = 1:numel(names)
    assignin('caller', names{i}, voltImg_8VisStim.(names{i}));
end
nConds = length(outParams.power);

names = fieldnames(ephys);
for i = 1:numel(names)
    assignin('caller', names{i}, ephys.(names{i}));
end

Fs = voltImg_8VisStim.daqParams.Fs;
trialTime = voltImg_8VisStim.daqParams.maxSweepLengthSec;
nTrials = length(voltImg_8VisStim.trialCond); % ALTERNATIVELY "length(find(cellfun(@isempty, ExpStruct.inputs)==0))". Instead of "length(ExpStruct.inputs)" this puts out true number of trials successfully recorded
powers = voltImg_8VisStim.outParams.power; % ALTERNATIVELY "unique(ExpStruct.trialCond)" what powers were used
nConds = length(voltImg_8VisStim.outParams.power); % ALTERNATIVELY "length(unique(ExpStruct.trialCond))" total number of powers used
nHolos = voltImg_8VisStim.holoStimParams.nHolos; % number of holograms in grid
pulseDurs = unique(voltImg_8VisStim.outParams.pulseDur);
nPulses = unique(voltImg_8VisStim.outParams.nPulses);
ipi = voltImg_8VisStim.outParams.ipi;
totalPulses = nHolos*nPulses;
% SpotCoordinates = SortedData.holoRequest.targets;
startTime = (voltImg_8VisStim.holoStimParams.startTime)/1000;
imagesIndex = voltImg_8VisStim.imagesIndex;
UpOrDown = voltImg_8VisStim.UpOrDown;
ephysFilePath = voltImg_8VisStim.ephysFilePath;
ImgsFilePath = voltImg_8VisStim.ImgsFilePath;

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