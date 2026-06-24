%% Voltage Imaging Slice Test Analysis Code
% Voltage imaging plots now show 95% CI
% Changed the indexing method for image folder so that the code
% automatically avoids hidden files and non-image files.
% Remember: dF/F =( F(t) - F0)/F0 
%%
clear all
% close all

%% Load files and setup
% Step 1: Read the ephys file
% ephysFilePath = char(uigetdir('D:\Data\Voltage Imaging\voltImg_slice_test\Ephys'));
% ephysFilePath = char(uigetdir('E:\Data\Voltage Imaging\dvTest\Ephys Data')); % Select and set root folder where all experiments with cells you want to analyze are located
% ephysFilePath = char(uigetdir('E:\Voltage Imaging\VoltImg_slice_test\Ephys Data')); 
ephysFilePath = char(uigetdir('/Volumes/ExData2/Voltage Imaging/VoltImg_slice_test/Ephys Data')); % Select and set root folder where all experiments with cells you want to analyze are located
ephysFileDir = dir(ephysFilePath);
load([ephysFileDir(3).folder, '/', ephysFileDir(3).name]);

% Get dv parameters and condition sequence
dvCondSequence = ExpStruct.dvStepParams.dvCondSequence; % sequence of randomized dv trials
dvToTest = ExpStruct.dvStepParams.dvToTest; % dv steps to be simulated by current injection
nConds = length(unique(dvCondSequence)); % number of conditions (dv steps)

pulseStart = ExpStruct.dvStepParams.pulseStart;
sweepDur = ExpStruct.dvStepParams.sweepDur;
nPulses = ExpStruct.dvStepParams.nPulses;
pulseFreq = ExpStruct.dvStepParams.pulseFreq;
imagingFreq = 330.22;
Fs = ExpStruct.Fs;

% Step 2: Identify the folder containing the imaging files correspdonding to the electrophysiology recording
% ImgsFilePath = char(uigetdir('D:\Data\Voltage Imaging\voltImg_slice_test\Imaging'));
% ImgsFilePath = char(uigetdir('E:\Data\Voltage Imaging\dvTest\Imaging Data')); % Select and set root folder where all experiments with cells you want to analyze are located
% ImgsFilePath = char(uigetdir('E:\Voltage Imaging\VoltImg_slice_test\Imaging Data')); 
ImgsFilePath = char(uigetdir('/Volumes/ExData2/Voltage Imaging/VoltImg_slice_test/Imaging Data')); % Select and set root folder where all experiments with cells you want to analyze are located
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
dvCondSequence = dvCondSequence(1:(length(imagesIndex)));
% dvCondSequence = ones(1, 2000)

% Step 3: Identify bad trials, to skip in analyses below
if isfield(ExpStruct.dvStepParams,'vsTest_inputs')
    vsTest_inputs = ExpStruct.dvStepParams.vsTest_inputs;
else
    vsTest_inputs = ExpStruct.dvStepParams.vsTest_inputs_Ch1;
end

baselineAllTrials = [];
excludeTrials = [];
for tt = 1:size(vsTest_inputs, 2)
    baseline = mean(vsTest_inputs(1:0.001*pulseStart*Fs, tt));
    baselineAllTrials = [baselineAllTrials, baseline];
    
    if baselineAllTrials(tt) > -55
        excludeTrials = [excludeTrials, tt];
    end
end

% Step 3: Identify GEVI
UpOrDown = input('1 for upward GEVI (ASAP6, FORCE, etc.), 2 for downward GEVI (Jedi2Pd, ASAP5, 7y etc.) ', 's');

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

voltImgTest_Analysis.ephysData.Fs = Fs;
voltImgTest_Analysis.pulseParams.pulseStart = pulseStart;
voltImgTest_Analysis.pulseParams.sweepDur = sweepDur;
voltImgTest_Analysis.pulseParams.nPulses = nPulses;
voltImgTest_Analysis.pulseParams.pulseFreq = pulseFreq;
voltImgTest_Analysis.imagingFreq = imagingFreq;
voltImgTest_Analysis.dvCondSequence = dvCondSequence;
voltImgTest_Analysis.dvToTest = dvToTest;
voltImgTest_Analysis.Rinput = ExpStruct.dvStepParams.Rinput;
voltImgTest_Analysis.fRinput = ExpStruct.dvStepParams.fRinput;
voltImgTest_Analysis.ephysData.vsTest_inputs = vsTest_inputs;
voltImgTest_Analysis.ephysData.baselineAllTrials = baselineAllTrials;
voltImgTest_Analysis.ephysData.excludeTrials = excludeTrials;

%% Calculate ROI mask
% Step 1: Calculate mask for ROI using spiking resp
% The logic here is that downward GEVI dim with depolarization, so
% efficient ROI selection with downward GEVI should be done on all trials
% that are NOT spiking (spiking = dimmmer cell).
% randperm(length(imagesIndex), 400);

maxDvTrials = [];
if UpOrDown == '1'
    maxDvTrials = find(dvCondSequence == max(unique(dvCondSequence)));
elseif UpOrDown == '2'
    maxDvTrials = find(dvCondSequence == min(unique(dvCondSequence)));
end

if length(maxDvTrials) > 100 
    maxDvTrials = randsample(maxDvTrials, 100);
end

maxDvStack = [];
counter = 0;
for tt = 1:length(maxDvTrials) 
    counter = counter+1;
    disp(['Trial number: ', num2str(counter), ' of ', num2str(length(maxDvTrials))]);
    if ismember(maxDvTrials(tt), excludeTrials)
        continue
    end
    currImgPath = [ImgfolderContents(imagesIndex(maxDvTrials(tt))).folder, '/', ImgfolderContents(imagesIndex(maxDvTrials(tt))).name];
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
    
%     % Preallocate the image stack
%     imageStack = zeros(info(1).Height, info(1).Width, numFrames);
% 
%     % Read each frame and store in the stack
%     for frameIndex = 1:numFrames
%         imageStack(:,:,frameIndex) = imread(currImgPath, 'Index', frameIndex, 'Info', info);
%     end
% %     for ff = 1:size(imageStack, 3); figure(101010); imagesc(imageStack(:, :, ff)); caxis([min(min(min(imageStack(:,:,:)))), max(max(max(imageStack(:,:,:))))]); axis equal; axis image; pause; end
%     maxDvStack = cat(3, maxDvStack, imageStack);
end

meanMaxDvStack = mean(maxDvStack(:, :, :), 3);
meanFluorMaxDvStack = meanMaxDvStack;

greenScale = [zeros(256,1), (0:255)'/255, zeros(256,1)];
figure(10); clf; imagesc(meanFluorMaxDvStack); colormap(greenScale); axis equal; axis image; axis off; colorbar off; %colorbar; set(gca, 'fontsize', 12); 

% figure(10); clf; imagesc(meanFluorMaxDvStack); axis equal; axis image; colorbar; set(gca, 'fontsize', 12); %axis off; colorbar off;% caxis([-7 -4]);

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
imagesc(roiMeanMaxDvStack); axis equal; axis image; axis off; colorbar off; %colorbar;

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
imagesc(roiStack); axis equal; axis image; axis off; colorbar off;
[roiX, roiY] = find(roiStack);

voltImgTest_Analysis.maxDvTrials = maxDvTrials;
voltImgTest_Analysis.meanFluorMaxDvStack = meanFluorMaxDvStack;
voltImgTest_Analysis.roiMeanMaxDvStack = roiMeanMaxDvStack;
voltImgTest_Analysis.roiX = roiX;
voltImgTest_Analysis.roiY = roiY;

%% Calculate df and mean df for all trials across all conditions
% Step 1: Calculate mean df across all trials across all conditions
cutOffFreq = 50;   % Cutoff frequency
[b, a] = butter(4, cutOffFreq/(imagingFreq/2)); % 4th order Butterworth filter

dfAllConds = cell(nConds, 1);
dffAllConds = cell(nConds, 1);
filtdffAllConds = cell(nConds, 1);
traceAllConds = cell(nConds, 1);
    
counter = 0;
for tt = 1:size(vsTest_inputs, 2)
    counter = counter+1;
    disp(['Trial number: ', num2str(counter), ' of ', num2str(size(vsTest_inputs, 2))]);
%     if ismember(tt, excludeTrials)
%         continue
%     end
    
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
%         [cleanFrame] = VoltImg_mapping_removeArtifact(imageStack(:,:,frameIndex));
%         imageStack(:,:,frameIndex) = cleanFrame; %replace the raw frame with cleaned up frame where artifact-corrupt lines are NaN'd    
    end
    
%     for tt = 1:size(imageStack, 3); figure(101010); imagesc(imageStack(:, :, tt)); caxis([min(min(min(imageStack(:,:,:)))), max(max(max(imageStack(:,:,:))))]); axis equal; axis image; pause; end
   
	% Set baseline(f0), this is calculated off ROI pixels and during baseline before the pulses 
if dvCondSequence(tt) == max(unique(dvCondSequence))
    roiBaselineAllPixels = imageStack(roiX, roiY, floor((pulseStart-pulseStart/2)/1000*imagingFreq:pulseStart/1000*imagingFreq)); %1:ceil(pulseStart/1000*imagingFreq));
    roiBaselineMean = nanmean(roiBaselineAllPixels, 'all');    %#ok<NANMEAN>
else
    roiBaselineAllPixels = imageStack(roiX, roiY, floor((pulseStart-pulseStart/2)/1000*imagingFreq:pulseStart/1000*imagingFreq)); %1:ceil(pulseStart/1000*imagingFreq));
    roiBaselineMean = nanmean(roiBaselineAllPixels, 'all'); %#ok<NANMEAN>
end

    df = [];
    dff = [];    
    for ff = 1:size(imageStack, 3)
        % ROI pixels for this frame, and calculate their mean
        currFrameRoi = imageStack(roiX, roiY, ff);
        currFrameRoiMean = nanmean(currFrameRoi, 'all');
        % Calculate df
        intensityChange = currFrameRoiMean - roiBaselineMean;
        df = [df; intensityChange];
        dff = [dff; intensityChange/roiBaselineMean];
    end
    
    if UpOrDown == '2'
        df = -df;
        dff = -dff;
    elseif UpOrDown =='1'
        df = df;
        dff = dff;
    end
    
    % Use these lines if I need to further filter the df or dff trace
    filtdff = filter(b, a, dff);

    % dfAllConds{dvCondSequence(tt), 1} = [dfAllConds{dvCondSequence(tt), 1}, df];
    % dffAllConds{dvCondSequence(tt), 1} = [dffAllConds{dvCondSequence(tt), 1}, dff];
    % filtdffAllConds{dvCondSequence(tt), 1} = [filtdffAllConds{dvCondSequence(tt), 1}, filtdff];
    % traceAllConds{dvCondSequence(tt), 1} = [traceAllConds{dvCondSequence(tt), 1}, vsTest_inputs(:, tt)];

    dfAllConds{dvCondSequence(tt), 1}(:, end+1) = df;
    dffAllConds{dvCondSequence(tt), 1}(:, end+1)  = dff;
    filtdffAllConds{dvCondSequence(tt), 1}(:, end+1)  = filtdff;
    traceAllConds{dvCondSequence(tt), 1}(:, end+1)  = vsTest_inputs(:, tt);

end

meanDfAllConds = [];
meanDffAllConds = [];
meanFiltDffAllConds = [];
meanTraceAllConds = [];
for cc = 1:nConds
    meanDfAllConds = [meanDfAllConds, mean(dfAllConds{cc}, 2)];
    meanDffAllConds = [meanDffAllConds, mean(dffAllConds{cc}, 2)];
    meanFiltDffAllConds = [meanFiltDffAllConds, mean(filtdffAllConds{cc}, 2)];
    meanTraceAllConds = [meanTraceAllConds, mean(traceAllConds{cc}, 2)];
end

% Calculate peak poststim df/f
peakDff = [];
peakFiltDff = [];
for cc = 1:nConds 
    for pp = 1:nPulses
        peakDff{cc}(:, pp) = max(meanDffAllConds((floor(pulseStart+((1000/pulseFreq)*(pp-1)))/1000*imagingFreq):floor((pulseStart+((1000/pulseFreq)*pp))/1000*imagingFreq), cc));
        peakFiltDff{cc}(:, pp) = max(meanFiltDffAllConds((floor(pulseStart+((1000/pulseFreq)*(pp-1)))/1000*imagingFreq):floor((pulseStart+((1000/pulseFreq)*pp))/1000*imagingFreq), cc));
    end
end

CIDfAllConds = [];
CIDffAllConds = [];
CIFiltDffAllConds = [];
for cc = 1:nConds 
    confidence_level = 0.95;
    meansDff = mean(dffAllConds{cc}, 2);
    std_errorsDff = std(dffAllConds{cc}, 0, 2) / sqrt(size(dffAllConds{cc}, 2));
    t_scoreDff = tinv((1 + confidence_level) / 2, size(dffAllConds{cc}, 2) - 1);
    margin_of_errorDff = t_scoreDff * std_errorsDff;
    lower_boundsDff = meansDff - margin_of_errorDff;
    upper_boundsDff = meansDff + margin_of_errorDff;

    meansDf = mean(dfAllConds{cc}, 2);
    std_errorsDf = std(dfAllConds{cc}, 0, 2) / sqrt(size(dfAllConds{cc}, 2));
    t_scoreDf = tinv((1 + confidence_level) / 2, size(dfAllConds{cc}, 2) - 1);
    margin_of_errorDf = t_scoreDf * std_errorsDf;
    lower_boundsDf = meansDf - margin_of_errorDf;
    upper_boundsDf = meansDf + margin_of_errorDf;

    meansFiltDff = mean(filtdffAllConds{cc}, 2);
    std_errorsFiltDff = std(filtdffAllConds{cc}, 0, 2) / sqrt(size(filtdffAllConds{cc}, 2));
    t_scoreFiltDff = tinv((1 + confidence_level) / 2, size(filtdffAllConds{cc}, 2) - 1);
    margin_of_errorFiltDff = t_scoreFiltDff * std_errorsFiltDff;
    lower_boundsFiltDff = meansFiltDff - margin_of_errorFiltDff;
    upper_boundsFiltDff = meansFiltDff + margin_of_errorFiltDff;

    if UpOrDown == '2'
        CIDffAllConds{cc, 1} = [lower_boundsDff, upper_boundsDff];
        CIDfAllConds{cc, 1} = [lower_boundsDf, upper_boundsDf];
        CIFiltDffAllConds{cc, 1} = [lower_boundsFiltDff, upper_boundsFiltDff];
    elseif UpOrDown =='1'
        CIDffAllConds{cc, 1} = [lower_boundsDff, upper_boundsDff];
        CIDfAllConds{cc, 1} = [lower_boundsDf, upper_boundsDf];
        CIFiltDffAllConds{cc, 1} = [lower_boundsFiltDff, upper_boundsFiltDff];
    end
end

voltImgTest_Analysis.maxDvStack = maxDvStack;
voltImgTest_Analysis.meanFluorMaxDvStack = meanFluorMaxDvStack;
voltImgTest_Analysis.stdFluor = stdFluor;
voltImgTest_Analysis.meanFluor = meanFluor;
voltImgTest_Analysis.cutOffFluor = cutOffFluor; 
voltImgTest_Analysis.roiStack = roiStack;
voltImgTest_Analysis.dfAllConds = dfAllConds;
voltImgTest_Analysis.dffAllConds = dffAllConds;
voltImgTest_Analysis.filtdffAllConds = filtdffAllConds;
voltImgTest_Analysis.traceAllConds = traceAllConds;
voltImgTest_Analysis.meanDffAllConds = meanDffAllConds;
voltImgTest_Analysis.meanDfAllConds = meanDfAllConds;
voltImgTest_Analysis.meanFiltDffAllConds = meanFiltDffAllConds;
voltImgTest_Analysis.meanTraceAllConds = meanTraceAllConds;
voltImgTest_Analysis.peakDff = peakDff;
voltImgTest_Analysis.peakFiltDff = peakFiltDff;
voltImgTest_Analysis.CIDffAllConds = CIDffAllConds;
voltImgTest_Analysis.CIDfAllConds = CIDfAllConds;
voltImgTest_Analysis.CIFiltDffAllConds = CIFiltDffAllConds;

%% Align imaging with ephys traces (one figure, subplots per condition)
nImgFrames = size(meanFiltDffAllConds, 1);
nEphysSamps = size(meanTraceAllConds, 1);
tImg = linspace(0, nImgFrames / imagingFreq, nImgFrames);
tEphys = linspace(0, nEphysSamps / Fs, nEphysSamps);
ephysBlStart = 1;
ephysBlEnd = pulseStart/1000*20000;

nColsSub = ceil(sqrt(nConds));
nRowsSub = ceil(nConds / nColsSub);
figure(31);
clf

for cc = 1:nConds
    subplot(nRowsSub, nColsSub, cc)
    hold on

    vmMean = mean(meanTraceAllConds(ephysBlStart:ephysBlEnd, cc), 'omitnan');
    vmTrace = meanTraceAllConds(:, cc) - vmMean;

    yyaxis left
    ciMat = CIDffAllConds{cc, 1};
    nCi = size(ciMat, 1);
    tCi = linspace(0, nCi / imagingFreq, nCi);
    ciLo = ciMat(:, 1) * 100;
    ciHi = ciMat(:, 2) * 100;
    fill([tCi, fliplr(tCi)], [ciLo', fliplr(ciHi')], [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
    plot(tCi, ciLo, '-', 'LineWidth', 1, 'Color', [0.9, 0.9, 0.9]);
    plot(tCi, ciHi, '-', 'LineWidth', 1, 'Color', [0.9, 0.9, 0.9]);
    plot(tImg, meanFiltDffAllConds(:, cc) * 100, '-', 'LineWidth', 2, 'Color', 'g');
    yline(0, '-', 'LineWidth', 1.5, 'Color', [0.9, 0.9, 0.9]);
    ylabel('dF/F (%)');

    yyaxis right
    plot(tEphys, vmTrace, '-', 'LineWidth', 2, 'Color', [0.8, 0.8, 0.8]);
    ylabel('Vm');

    for pp = 1:nPulses
        xline((pulseStart / 1000) + (pp - 1) / pulseFreq, '--', 'LineWidth', 1.5, 'Color', [0.5, 0.5, 0.5]);
    end

    % Same vertical position for y=0 on left (dF/F) and right (Vm): share fraction t0 of axis span
    yPad = 0.05;
    Lmin = min([0; ciLo(:); ciHi(:); meanFiltDffAllConds(:, cc) * 100]);
    Lmax = max([0; ciLo(:); ciHi(:); meanFiltDffAllConds(:, cc) * 100]);
    Rmin = min([0; vmTrace(:)]);
    Rmax = max([0; vmTrace(:)]);
    yyaxis left
    yLauto = ylim;
    if yLauto(1) < 0 && yLauto(2) > 0
        t0 = (0 - yLauto(1)) / (yLauto(2) - yLauto(1));
    else
        t0 = 0.5;
    end
    t0 = min(max(t0, 0.02), 0.98);
    spanL = max(-Lmin / t0, Lmax / (1 - t0)) * (1 + yPad);
    spanR = max(-Rmin / t0, Rmax / (1 - t0)) * (1 + yPad);
    ylim([-t0 * spanL, (1 - t0) * spanL]);
    yyaxis right
    ylim([-t0 * spanR, (1 - t0) * spanR]);

    titleStr = ['dV = ', num2str(abs(mean(meanTraceAllConds(1:(pulseStart / 1000) * Fs, cc))) + (max(meanTraceAllConds(((pulseStart / 1000) * Fs):((pulseStart / 1000) * Fs + ((1000 / pulseFreq) / 1000 * Fs)), cc)))), 'mV,', ' (1st pulse) ', ...
        num2str(abs(mean(meanTraceAllConds(1:(pulseStart / 1000) * Fs, cc))) + (max(meanTraceAllConds(((pulseStart / 1000) * Fs + ((1000 / pulseFreq) / 1000 * Fs)):end, cc)))), 'mV,', ' (2nd pulse)'];
    title(titleStr, 'FontSize', 10);
    xlabel('Time (s)');
    set(gca, 'FontSize', 16);
    hold off
end

%% Show some individual trials
for cc = 5 %1:nConds
    for tt = 1:size(dffAllConds{cc}, 2)
        figure(300+cc*10+tt)
        % set(gcf,'Position',[100 100 1000 1000])
        clf
        hold on
    %     axes('Position',[.70 .12 .2 .2]);
    %     box on
%         plot(linspace(0, size(traceAllConds{:, cc}(:, tt), 1)/Fs, size(traceAllConds{:, cc}(:, tt), 1)), traceAllConds{:, cc}(:, tt), 'linewidth', 2, 'color', [0.7 0.7 0.7]);
        plot(linspace(0, size(traceAllConds{cc}(:, tt), 1)/Fs, size(traceAllConds{cc}(:, tt), 1)), traceAllConds{cc}(:, tt)+abs(mean(traceAllConds{cc}(1:(pulseStart/1000)*Fs, cc)))-70, 'linewidth', 2, 'color', [0.8 0.8 0.8]);
        gca;
        set(gca,'xtick',[], 'fontsize', 18);
        % ylim([-1, max(meanTraceAllConds(:, cc)+abs(mean(meanTraceAllConds(1:(pulseStart/1000)*Fs, cc))))]);
    %     ylim([-0.5, 1])
        % ylim([min(meanTraceAllConds(:, cc)), max(meanTraceAllConds(:, cc))]);
        % ylim([min(min(meanDffAllConds)) max(max(meanDffAllConds))]);
        % ylim([-80 120]);
        ylabel('Voltage (mV)');
        xticks([0:0.25:size(traceAllConds{cc}(:, tt), 1)/Fs]);
        % ylabel('mV');
        hold off
        % ylim([-20 100]);
    
    %     plot ephys trace
        yyaxis right
        for pp = 1:1
            xline((pulseStart/1000)+(pp-1)/pulseFreq, '--', 'LineWidth', 1.5, 'color', [.5 .5 .5]);
        end  
        
        % [negpks1, peakLociX] = findpeaks(dffAllConds{:, cc}(:, tt), 'MinPeakHeight', 70);
        
        % plot voltage trace
%         yyaxis left
        plot(linspace(0, size(dffAllConds{cc}(:, tt), 1)/imagingFreq, size(dffAllConds{cc}(:, tt), 1)), dffAllConds{cc}(:, tt), '-', 'linewidth', 2, 'color', [0 0.75 0]);
        % scatter(peakLociX/imagingFreq , negpks1+2, 'k', '*');
        ylabel('dF/F (%)');
        xlim([0.1 sweepDur]);
    %     ylim([-0.5 0.5]);
        set(gca, 'fontsize', 16);
        xlabel('Time (s)');
        
        % show line at dff = 0
        yline(0, '-', 'LineWidth', 1.5, 'color', [0.9 0.9 0.9]);
        hold on
    
        % plot([0.15; 0.25], [0.05; 0.05], '-k', 'LineWidth', 2); % x axis scale
        % plot([0.15; 0.15], [0.15; 0.05], '-k', 'LineWidth', 2); % y axis scale
        
        % axis off
        hold off
    end
end

%% Save Analysis Results
% directory = 'E:\Voltage Imaging\VoltImg_slice_test\Analysis Results';
directory = '/Volumes/ExData2/Voltage Imaging/VoltImg_slice_test/Analysis Results/dvTest_Analysis_AH';
voltImgTest_Analysis.mouseID = ['voltImgTest_Analysis_', num2str(ExpStruct.mouseID)];
fileName = [num2str(voltImgTest_Analysis.mouseID), '.mat'];
save(fullfile(directory, fileName), 'voltImgTest_Analysis', '-v7.3');

TimeNow = datetime('now', 'Format', 'dd-MMM-yyyy HH:mm:ss');
disp(['finished saving at: ' char(TimeNow)])

%% Load Analysis Results
% Run this section after loading the specific cell analysis file, then
% re-run the above sections to regenerate figures

names = fieldnames(voltImgTest_Analysis);
for i = 1:numel(names)
    assignin('caller', names{i}, voltImgTest_Analysis.(names{i}));
end
nConds = length(dfAllConds);

names = fieldnames(pulseParams);
for i = 1:numel(names)
    assignin('caller', names{i}, pulseParams.(names{i}));
end

names = fieldnames(ephysData);
for i = 1:numel(names)
    assignin('caller', names{i}, ephysData.(names{i}));
end

UpOrDown = input('1 for upward GEVI (ASAP6, FORCE, etc.), 2 for downward GEVI (Jedi2Pd, ASAP5, 7y etc.) ', 's');

%% %%%%%%%%%%%%% %%
%% %%%%%%%%%%%%% %%
%% %%%%%%%%%%%%% %%
%% dF over voltage
figure(40);
clf
hold on
for cc = 1:nConds
    scatter(abs(mean(meanTraceAllConds(1:(ExpStruct.dvStepParams.pulseStart/1000)*Fs, cc)))+(max(meanTraceAllConds(((ExpStruct.dvStepParams.pulseStart/1000)*Fs):((ExpStruct.dvStepParams.pulseStart/1000)*Fs+((1000/ExpStruct.dvStepParams.pulseFreq)/1000*Fs)), cc))),...
        max(-meanDfAllConds(ceil((ExpStruct.dvStepParams.pulseStart/1000)*imagingFreq):floor(((ExpStruct.dvStepParams.pulseStart/1000)*imagingFreq+((1000/ExpStruct.dvStepParams.pulseFreq)/1000*imagingFreq))), cc)), 50, 'k', 'filled');
    scatter(abs(mean(meanTraceAllConds(1:(ExpStruct.dvStepParams.pulseStart/1000)*Fs, cc)))+(max(meanTraceAllConds(((ExpStruct.dvStepParams.pulseStart/1000)*Fs+((1000/ExpStruct.dvStepParams.pulseFreq)/1000*Fs)):end, cc))),...
        max(-meanDfAllConds(ceil(((ExpStruct.dvStepParams.pulseStart/1000)*imagingFreq+((1000/ExpStruct.dvStepParams.pulseFreq)/1000*imagingFreq))):end, cc)), 50, 'k', 'filled');
end
xlabel('dV (mV)')
ylabel('dF')
axis square;
set(gca,'linewidth',1.5, 'fontsize', 12)
hold off

%% In case dvCondSequence is screwed up
temp = [];
temp = vsTest_inputs/2;
dvCondSequence = [];
for tt = 1:size(temp, 2)
    if max(temp(:, tt))> -20
        dvCondSequence(tt) = 2;
    else
        dvCondSequence(tt) = 1;
    end
end

%%
clear all
%% Load files and setup
% Step 1: Read the ephys file
analyzedFilePath = char(uigetdir('D:\Data\Voltage Imaging\voltImg_slice_test\Analysis Results')); % Select and set root folder where all experiments with cells you want to analyze are located
% ephysFilePath = char(uigetdir('/Volumes/Elements/Data/Voltage Imaging/VoltImg_slice_test/Ephys Data')); % Select and set root folder where all experiments with cells you want to analyze are located
analyzedFileDir = dir(analyzedFilePath);
% load([analyzedFileDir(3).folder, '/', analyzedFileDir(3).name]);

% Step 2: Avoid hidden files and non image files
fileNames = [];
fileType = '.mat';
for ii = 1:length(analyzedFileDir)
    % Check if the entry is a regular file (not a directory) and its name doesn't start with a period
    if ~analyzedFileDir(ii).isdir && ~startsWith(analyzedFileDir(ii).name, '.') && endsWith(analyzedFileDir(ii).name, fileType)
        % Add the file name to the cell array
        fileNames{ii, 1} = analyzedFileDir(ii).name;
    end
end
filesIndex = find(~cellfun(@isempty, fileNames));

peakDffAllCells = [];
peakdVAllCells = [];
for ff = 1:length(filesIndex)
    load([analyzedFileDir(filesIndex(ff)).folder, '/', analyzedFileDir(filesIndex(ff)).name]);
        Fs = voltImgTest_Analysis.ephysData.Fs;
        pulseStart = voltImgTest_Analysis.pulseParams.pulseStart;
        sweepDur = voltImgTest_Analysis.pulseParams.sweepDur;
        nPulses = voltImgTest_Analysis.pulseParams.nPulses;
        pulseFreq = voltImgTest_Analysis.pulseParams.pulseFreq;
        imagingFreq = voltImgTest_Analysis.imagingFreq;
        dvCondSequence = voltImgTest_Analysis.dvCondSequence;
        dvToTest = voltImgTest_Analysis.dvToTest;
        Rinput = voltImgTest_Analysis.Rinput;
        fRinput = voltImgTest_Analysis.fRinput;
        vsTest_inputs = voltImgTest_Analysis.ephysData.vsTest_inputs;
        baselineAllTrials = voltImgTest_Analysis.ephysData.baselineAllTrials;
        excludeTrials = voltImgTest_Analysis.ephysData.excludeTrials;

%         maxDvTrials = voltImgTest_Analysis.maxDvTrials;
        meanFluorMaxDvStack = voltImgTest_Analysis.meanFluorMaxDvStack;
%         roiMeanMaxDvStack = voltImgTest_Analysis.roiMeanMaxDvStack;
%         roiX = voltImgTest_Analysis.roiX;
%         roiY = voltImgTest_Analysis.roiY;

        maxDvStack = voltImgTest_Analysis.maxDvStack;
        meanFluorMaxDvStack = voltImgTest_Analysis.meanFluorMaxDvStack;
        stdFluor = voltImgTest_Analysis.stdFluor;
        meanFluor = voltImgTest_Analysis.meanFluor;
        cutOffFluor = voltImgTest_Analysis.cutOffFluor;
        roiStack = voltImgTest_Analysis.roiStack;
        dfAllConds = voltImgTest_Analysis.dfAllConds;
        traceAllConds = voltImgTest_Analysis.traceAllConds;
        meanDfAllConds = voltImgTest_Analysis.meanDfAllConds;
        meanTraceAllConds = voltImgTest_Analysis.meanTraceAllConds;
%         peakDff = voltImgTest_Analysis.peakDff;
        CIDfAllConds = voltImgTest_Analysis.CIDfAllConds;
        
        nConds = length(unique(dvCondSequence));
        
        % Calculate peak poststim df/f
        peakDff = [];
        peakdV = [];
        for cc = 1:nConds 
            for pp = 1:nPulses
                peakDff{cc}(:, pp) = max(meanDfAllConds((floor(pulseStart+((1000/pulseFreq)*(pp-1)))/1000*imagingFreq):floor((pulseStart+((1000/pulseFreq)*pp))/1000*imagingFreq), cc));
                peakdV{cc}(:, pp) = abs(mean(meanTraceAllConds(1:(pulseStart/1000)*Fs, cc))) + max(meanTraceAllConds((floor(pulseStart+((1000/pulseFreq)*(pp-1)))/1000*Fs):floor((pulseStart+((1000/pulseFreq)*pp))/1000*Fs), cc));
            end
        end
        
        peakDffAllCells{ff} = peakDff;
        peakdVAllCells{ff} = peakdV;
end

%%
% jediPeakDffs = peakDffAllCells;
% jediPeakdVs = peakdVAllCells;
% asap6PeakDffs = peakDffAllCells;
% asap6PeakdVs = peakdVAllCells;
% asap5PeakDffs = peakDffAllCells;
% asap5PeakdVs = peakdVAllCells;

peakDffComparison.jediPeakDffs = jediPeakDffs;
peakDffComparison.jediPeakdVs = jediPeakdVs;
peakDffComparison.asap6PeakDffs = asap6PeakDffs;
peakDffComparison.asap6PeakdVs = asap6PeakdVs;
peakDffComparison.asap5PeakDffs = asap5PeakDffs;
peakDffComparison.asap5PeakdVs = asap5PeakdVs;

%%
jediDfftoDV = cell(1, 2);
for nn = 1:size(jediPeakDffs, 2) %cell
    for cc = 1:size(jediPeakDffs{nn}, 2) %condition
        for pp = 1:size(jediPeakDffs{nn}{cc}, 2)
            jediDfftoDV{1, 1} = [jediDfftoDV{1, 1}; jediPeakDffs{nn}{cc}(pp)];
            jediDfftoDV{1, 2} = [jediDfftoDV{1, 2}; jediPeakdVs{nn}{cc}(pp)];
        end
    end
end
jediDfftoDVMerged = [jediDfftoDV{1, 1}, jediDfftoDV{1, 2}];

asap6DfftoDV = cell(1, 2);
for nn = 1:size(asap6PeakDffs, 2) %cell
    for cc = 1:size(asap6PeakDffs{nn}, 2) %condition
        for pp = 1:size(asap6PeakDffs{nn}{cc}, 2)
            asap6DfftoDV{1, 1} = [asap6DfftoDV{1, 1}; asap6PeakDffs{nn}{cc}(pp)];
            asap6DfftoDV{1, 2} = [asap6DfftoDV{1, 2}; asap6PeakdVs{nn}{cc}(pp)];
        end
    end
end
asap6DfftoDVMerged = [asap6DfftoDV{1, 1}, asap6DfftoDV{1, 2}];

asap5DfftoDV = cell(1, 2);
for nn = 1:size(asap5PeakDffs, 2) %cell
    for cc = 1:size(asap5PeakDffs{nn}, 2) %condition
        for pp = 1:size(asap5PeakDffs{nn}{cc}, 2)
            asap5DfftoDV{1, 1} = [asap5DfftoDV{1, 1}; asap5PeakDffs{nn}{cc}(pp)];
            asap5DfftoDV{1, 2} = [asap5DfftoDV{1, 2}; asap5PeakdVs{nn}{cc}(pp)];
        end
    end
end
asap5DfftoDVMerged = [asap5DfftoDV{1, 1}, asap5DfftoDV{1, 2}];

figure(100);
clf
hold on
scatter(jediDfftoDVMerged(:,2),jediDfftoDVMerged(:,1), 'b');
    JediFit = fitlm(jediDfftoDVMerged(:,2), jediDfftoDVMerged(:,1), 'Intercept', false);
%     JediFitLine = plot(JediFit, 'Marker', 'none', 'Color', 'b');
scatter(asap6DfftoDVMerged(:,2),asap6DfftoDVMerged(:,1), 'g');
    asap6Fit = fitlm(asap6DfftoDVMerged(:,2), asap6DfftoDVMerged(:,1), 'Intercept', false);
%     asap6FitLine = plot(asap6Fit, 'Marker', 'none', 'Color', 'g');
scatter(asap5DfftoDVMerged(:,2),asap5DfftoDVMerged(:,1), 'r');
    asap5Fit = fitlm(asap5DfftoDVMerged(:,2), asap5DfftoDVMerged(:,1), 'Intercept', false);
%     asap5FitLine = plot(asap5Fit , 'Marker', 'none', 'Color', 'r');
set(gca, 'LineWidth', 1)
axis square;
ylabel('dF/F(%)');
xlabel('dV');
title(['Comparison of peak dF/F']);
legend('Jedi2P', 'ASAP6', 'ASAP5');
hold off

figure(101);
clf
hold on
scatter(jediDfftoDVMerged(:,2),jediDfftoDVMerged(:,1), 'b');
    JediFit = fitlm(jediDfftoDVMerged(:,2), jediDfftoDVMerged(:,1), 'Intercept', false);
%     JediFitLine = plot(JediFit, 'Marker', 'none', 'Color', 'b');
scatter(asap6DfftoDVMerged(:,2),asap6DfftoDVMerged(:,1), 'g');
    asap6Fit = fitlm(asap6DfftoDVMerged(:,2), asap6DfftoDVMerged(:,1), 'Intercept', false);
%     asap6FitLine = plot(asap6Fit, 'Marker', 'none', 'Color', 'g');
scatter(asap5DfftoDVMerged(:,2),asap5DfftoDVMerged(:,1), 'r');
    asap5Fit = fitlm(asap5DfftoDVMerged(:,2), asap5DfftoDVMerged(:,1), 'Intercept', false);
%     asap5FitLine = plot(asap5Fit , 'Marker', 'none', 'Color', 'r');
set(gca, 'LineWidth', 1);
axis square;
ylabel('dF/F(%)');
xlabel('dV');
xlim([0 10]);
ylim([0 1.5]);
title(['Comparison of peak dF/F at low dV']);
legend('Jedi2P', 'ASAP6', 'ASAP5');
hold off




%% Save Analysis Results
directory = 'E:\Voltage Imaging\dvTest\Analysis Results';
fileName = ['peakDffComparison', '.mat'];
save(fullfile(directory, fileName), 'peakDffComparison', '-v7.3');

TimeNow = datetime('now', 'Format', 'dd-MMM-yyyy HH:mm:ss');
disp(['finished saving at: ' char(TimeNow)])

