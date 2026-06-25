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
    if numFrames > imagingFreq*4 % Just take frames from first 4 seconds (if actual sweeps are longer than 4s), because not entire sweeps are needed
        numFrames = floor(imagingFreq*4);
    end

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
maxDvStackBrighter = meanFluorMaxDvStack + 150; % Brighter z-stack image if needed
maxDvStackBrighter(maxDvStackBrighter > 255) = 255; % Clamping for brightness effect

figure(10); clf; colormap('winter'); imagesc(maxDvStackBrighter); axis equal; axis image; colorbar; set(gca, 'fontsize', 12); % caxis([-7 -4]);
nCells = input('How many neurons to analyze?: ', 's');
% Hand select cell or area of interest, by freehand drawing
roiXAllCells = cell(nCells, 1);
roiYAllCells = cell(nCells, 1);

for cc = 1:nCells
    figure(11);
    set(gcf, 'Position',  [100, 100, 1800, 600])
    clf
    subplot(1,2,1)
    colormap('winter'); imagesc(maxDvStackBrighter); axis equal; axis image; colorbar; set(gca, 'fontsize', 12);
    roiX = []; roiY = [];
    roiHandSelect = drawfreehand;
    roiHandSelectMask = createMask(roiHandSelect);
    [roiX, roiY] = find(roiHandSelectMask);
    
    % Calculate mean fluorescence of the area of interest and then form into final ROI mask
    roiMeanMaxDvStack = zeros(size(maxDvStack, 1), size(maxDvStack, 2));
    for rr = 1:length(roiX)
        roiMeanMaxDvStack(roiX(rr), roiY(rr)) = mean(maxDvStack(roiX(rr), roiY(rr),:), 3);
    end
    hold on;
    
    subplot(1,2,2);
    imagesc(roiMeanMaxDvStack); axis equal; axis image; colorbar;
    
    stdFluor = std(nonzeros(roiMeanMaxDvStack));
    meanFluor = mean(nonzeros(roiMeanMaxDvStack));
    
    % Designate cutoff fluorescence for pixels to be selected for ROI
    cutOffFluor = meanFluor; %stdFluor*1 + meanFluor; % currently cutoff is 1 standard devs from mean fluorescence
    
    roiMeanMaxDvStack(roiMeanMaxDvStack <= cutOffFluor) = 0;
    roiMeanMaxDvStack(roiMeanMaxDvStack > 0) = 1;
    roiStack = roiMeanMaxDvStack;
    
    subplot(1,2,3);
    imagesc(roiStack); axis equal; axis image;
    [roiX, roiY] = find(roiStack); % XY coordinates of the pixels of interest
    hold off;

    roiXAllCells{cc} = roiX;
    roiYAllCells{cc} = roiY;
end

voltMapping.maxDvStack = maxDvStack;
voltMapping.meanFluorMaxDvStack = meanFluorMaxDvStack;
voltMapping.maxDvStackBrighter = maxDvStackBrighter;
voltMapping.roiMeanMaxDvStack = roiMeanMaxDvStack;
voltMapping.cutOffFluor = cutOffFluor;
voltMapping.roiXAllCells = roiXAllCells;
voltMapping.roiYAllCells = roiYAllCells;

%% Calculate Background region
% figure(10); clf; colormap('winter'); imagesc(meanFluorMaxDvStack); axis equal; axis image; colorbar; set(gca, 'fontsize', 12);% caxis([-7 -4]);
% 
% % Hand select cell or area of interest, by freehand drawing
% BkgrndX = []; BkgrndY = [];
% roiHandSelect = drawfreehand;
% roiHandSelectMask = createMask(roiHandSelect);
% [BkgrndX, BkgrndY] = find(roiHandSelectMask);
% 
% % Calculate mean fluorescence of the area of interest and then form into final ROI mask
% BkgrndMeanMaxDvStack = zeros(size(maxDvStack, 1), size(maxDvStack, 2));
% for rr = 1:length(BkgrndX)
%     BkgrndMeanMaxDvStack(BkgrndX(rr), BkgrndY(rr)) = mean(maxDvStack(BkgrndX(rr), BkgrndY(rr),:), 3);
% end
% 
% figure(11); clf;
% imagesc(BkgrndMeanMaxDvStack); axis equal; axis image; colorbar;
% 
% stdFluor = std(nonzeros(BkgrndMeanMaxDvStack));
% meanFluor = mean(nonzeros(BkgrndMeanMaxDvStack));
% 
% % Designate cutoff fluorescence for pixels to be selected for ROI
% bkgrndCutOffFluor = meanFluor; %stdFluor*1 + meanFluor; % currently cutoff is 1 standard devs from mean fluorescence
% 
% BkgrndMeanMaxDvStack(BkgrndMeanMaxDvStack <= bkgrndCutOffFluor) = 0;
% BkgrndMeanMaxDvStack(BkgrndMeanMaxDvStack > 0) = 1;
% roiStack = BkgrndMeanMaxDvStack;
% 
% figure(12); clf;
% imagesc(roiStack); axis equal; axis image;
% [BkgrndX, BkgrndY] = find(roiStack); % XY coordinates of the pixels of interest
% 
% voltMapping.BkgrndMeanMaxDvStack = BkgrndMeanMaxDvStack;
% voltMapping.bkgrndCutOffFluor = bkgrndCutOffFluor;
% voltMapping.BkgrndX = BkgrndX;
% voltMapping.BkgrndY = BkgrndY;

%% Calculate df and mean df for all holos sorted by trial conditions
% Use these lines if I need to further filter the df or dff trace
cutOffFreq = 50;   % Cutoff frequency
[b, a] = butter(4, cutOffFreq/(imagingFreq/2));  % 4th order Butterworth filter

% Preallocate for data for every cell, sorted by holograms, all trials, across conditions
holoSortedImagingCellNames = cell(nCells, 1);
filtHoloSortedImagingCellNames = cell(nCells, 1);
for nn = 1:nCells
    holoSortedImagingCellNames{nn} = ['holoSortedImagingAllTrials_', 'cell', num2str(nn)];
    filtHoloSortedImagingCellNames{nn} = ['filtHoloSortedImagingAllTrials_', 'cell', num2str(nn)];
    voltMapping.(holoSortedImagingCellNames{nn}) = cell(nConds, 1);
    voltMapping.(filtHoloSortedImagingCellNames{nn}) = cell(nConds, 1);
    
    for cc = 1:nConds
        voltMapping.(holoSortedImagingCellNames{nn}){cc} = cell(nHolos(cc), 1);
        voltMapping.(filtHoloSortedImagingCellNames{nn}){cc} = cell(nHolos(cc), 1);
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

    for cc = 1:nCells
        % Baselining to prestimulus at beginning of sweep: basic parameters
        startTimeImaging = floor(startTime*imagingFreq);  
        baselinePreImageStack = imageStack(roiXAllCells{cc}, roiYAllCells{cc}, 1:startTimeImaging); % image stack in ROI before start of stimulation
        roiBaselinePreImageStack = [];
        for ff = 1:size(baselinePreImageStack, 3)
            roiBaselinePreImageStack(ff, 1) = mean(mean(baselinePreImageStack(:, :, ff))); % Generate mean trace of ROI fluorescence during period before stimulation
        end
    
        % Baselining to prestimulus at beginning of sweep: parameters for moving window variance during pre-stim period
        windowLimits = [1/imagingFreq, size(roiBaselinePreImageStack, 1)/imagingFreq]; % in sec The range within each sweep where variance is sampled by 10 sample windows, in ms. Will leave it as 500 ms total width, starting from 110ms after start (after the step pulse).
        firstLimit = windowLimits(1)*imagingFreq; % Start of the first sampling window (in samples).
        windowTime = ceil(size(roiBaselinePreImageStack, 1)/imagingFreq*1000)/1000;
        segmentTime = 0.010; % time length of each window (s)
        numSegments = windowTime/segmentTime; % number of sample windows within the limit
        windowWidth = floor((windowLimits(2) - windowLimits(1))/numSegments*imagingFreq); % Span of each sample window (in samples).
    
        % Calculate variances/fanofactors in the prestimulus period using rolling window
        varBaselinePre = movvar(roiBaselinePreImageStack, windowWidth); % moving window variance calculated across prestimulus 
        fanoBaselinePre = varBaselinePre/(mean(roiBaselinePreImageStack));
        
        % % Whole sweep baseline Method1: Baseline entire trace to period before first stimulation by thresholding the fanofactor
        % % Fluorescence traces of pre-stim and post-stim periods
        % % Choose variance/fanofactor threshold
        % varThresholdPre = mean(varBaselinePre) - std(varBaselinePre); % variance threshold (mean-std) 
        % fanoThresholdPre = mean(fanoBaselinePre) - std(fanoBaselinePre); % fanofactor (variance/mean)
        % % Find variance/fanofactor beneath threshold
        % [vLowestPre, ~] = find(varBaseline < varThreshold);
        % [fanoLowestPre, ~] = find(fanoBaseline < fanoThreshold);
        % % Establish the baseline fluorescence based on all the points corresponding to the lowest moments of fano
        % roiBaselineMean = mean(roiBaselinePreImageStack(fanoLowestPre));
    
        % Whole sweep baseline Method2: Calculate baseline value before first stimulation based on quantile bottom 10% of fanofactor 
        % Step 1: Define quantile threshold (e.g., bottom 10%)
        q = 0.10;  % Change to 0.05 for bottom 5%, etc.
        quantileCutoff = quantile(fanoBaselinePre, q);
        % Step 2: Select low points
        [fanoLowestPre, ~] = find(fanoBaselinePre < quantileCutoff);
        roiBaselineThisHolo = mean(roiBaselinePreImageStack(fanoLowestPre)); % baseline value to subtract
    
        % Break apart the imageStack for this trial into stim windows and rearrange according to hologram sequence 
        holoSeqThisTrial = (unique(voltMapping.outParams.sequenceThisTrial{tt}, 'stable') - min(unique(voltMapping.outParams.sequenceThisTrial{tt}, 'stable')) + 1)'; % Hologram sequence for this trial    
        for hh = 1:nHolos(voltMapping.trialCond(tt, 1))
            windowLimits = []; firstLimit = []; windowTime = []; segmentTime = []; numSegments = []; windowWidth = [];
    
            % Extract the frames associated with the current hologram, including pre and post stimulation windows
            roiFramesThisHolo = imageStack(roiXAllCells{cc}, roiYAllCells{cc}, floor((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*imagingFreq):(ceil((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*imagingFreq)+ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)));
    
            % Stim window baselining method: Same as whole sweep baselining but baseline value is set to moment of least variability right before specific target/hologram
            % Image stack in ROI before stim of this target
            roiFramesThisHoloPreStim = roiFramesThisHolo(:, :, 1:preStimWindow/1000*imagingFreq);
        
            roiTraceThisHoloPreStim = [];
            for ff = 1:size(roiFramesThisHoloPreStim, 3)
                roiTraceThisHoloPreStim(ff, 1) = mean(mean(roiFramesThisHoloPreStim(:, :, ff))); % Generate mean trace of ROI fluorescence during period before stimulation
            end
        
            % Baselining to prestimulus at beginning of sweep: parameters for moving window variance during pre-stim period
            windowLimits = [1/imagingFreq, size(roiTraceThisHoloPreStim, 1)/imagingFreq]; % in sec The range within each sweep where variance is sampled by 10 sample windows, in ms. Will leave it as 500 ms total width, starting from 110ms after start (after the step pulse).
            firstLimit = windowLimits(1)*imagingFreq; % Start of the first sampling window (in samples).
            windowTime = ceil(size(roiTraceThisHoloPreStim, 1)/imagingFreq*1000)/1000;
            segmentTime = 0.010; % time length of each window (s)
            numSegments = windowTime/segmentTime; % number of sample windows within the limit
            windowWidth = floor((windowLimits(2) - windowLimits(1))/numSegments*imagingFreq); % Span of each sample window (in samples).
        
            % Calculate variances/fanofactors in the prestimulus period using rolling window
            varBaselineThisHolo = movvar(roiTraceThisHoloPreStim, windowWidth); % moving window variance calculated across prestimulus 
            fanoBaselineThisHolo = varBaselineThisHolo/(mean(roiTraceThisHoloPreStim));
    
            % Calculate baseline value from window before this hologram stimulation based on quantile bottom 10% of fanofactor 
            % Step 1: Define quantile threshold (e.g., bottom 10%)
            q = 0.10;  % Change to 0.05 for bottom 5%, etc.
            quantileCutoff = quantile(fanoBaselineThisHolo, q);
            % Step 2: Select low points
            [fanoLowestThisHolo, ~] = find(fanoBaselineThisHolo < quantileCutoff);
            roiBaselineThisHolo = mean(roiTraceThisHoloPreStim(fanoLowestThisHolo)); % baseline value to subtract
    
            dfThisHolo = [];
            dffThisHolo = [];
            for ff = 1:size(roiFramesThisHolo, 3)
                % ROI pixels for this frame
                currFrameRoi = roiFramesThisHolo(:, :, ff);
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
                voltMapping.(holoSortedImagingCellNames{cc}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, NaN(ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)+1, 1)];
                voltMapping.(filtHoloSortedImagingCellNames{cc}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [filtHoloSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, NaN(ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)+1, 1)];
            else
                %           holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, tt), 1}(:, tt) = dfThisHolo;
                voltMapping.(holoSortedImagingCellNames{cc}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, dffThisHolo];
                voltMapping.(filtHoloSortedImagingCellNames{cc}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [filtHoloSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, filtdffThisHolo];
            end        
        end % nHolos
    end % nCells
end % nTrials
