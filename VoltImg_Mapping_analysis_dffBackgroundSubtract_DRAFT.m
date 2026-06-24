%% Calculate df and mean df for all holos sorted by trial conditions
% Use these lines if I need to further filter the df or dff trace
input('Run df/f? (ctrl+c to stop!)');
cutOffFreq = 50;   % Cutoff frequency
[b, a] = butter(4, cutOffFreq/(imagingFreq/2));  % 4th order Butterworth filter

% Preallocate for data for every cell, sorted by holograms, all trials, across conditions
holoSortedImagingCellNames = cell(nCells, 1);
filtHoloSortedImagingCellNames = cell(nCells, 1);
holoSortedImagingCellNames_bkgrndsubtrct = cell(nCells, 1); 
filtHoloSortedImagingCellNames_bkgrndsubtrct = cell(nCells, 1);
for nn = 1:nCells
    holoSortedImagingCellNames{nn} = ['holoSortedImagingAllTrials_', 'cell', num2str(nn)];
    filtHoloSortedImagingCellNames{nn} = ['filtHoloSortedImagingAllTrials_', 'cell', num2str(nn)];
    holoSortedImagingCellNames_bkgrndsubtrct{nn} = ['holoSortedImagingAllTrials_bkgrndsubtrct_', 'cell', num2str(nn)];
    filtHoloSortedImagingCellNames_bkgrndsubtrct{nn} = ['filtHoloSortedImagingAllTrials_bkgrndsubtrct_', 'cell', num2str(nn)];   
    
    analysisStruct.(holoSortedImagingCellNames{nn}) = cell(nConds, 1);
    analysisStruct.(filtHoloSortedImagingCellNames{nn}) = cell(nConds, 1);
    analysisStruct.(holoSortedImagingCellNames_bkgrndsubtrct{nn}) = cell(nConds, 1);
    analysisStruct.(filtHoloSortedImagingCellNames_bkgrndsubtrct{nn}) = cell(nConds, 1);   
    
    for cc = 1:nConds
        analysisStruct.(holoSortedImagingCellNames{nn}){cc} = cell(nHolos(cc), 1);
        analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc} = cell(nHolos(cc), 1);
        analysisStruct.(holoSortedImagingCellNames_bkgrndsubtrct{nn}){cc} = cell(nHolos(cc), 1);
        analysisStruct.(filtHoloSortedImagingCellNames_bkgrndsubtrct{nn}){cc} = cell(nHolos(cc), 1);
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
        % Baselining to prestimulus at beginning of sweep: basic parameters
        startTimeImaging = floor(startTime*imagingFreq);  
        baselinePreImageStack = imageStack(roiXAllCells{nn}, roiYAllCells{nn}, 1:startTimeImaging); % image stack in ROI before start of stimulation
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
    
%         % Calculate variances/fanofactors in the prestimulus period using rolling window
%         varBaselinePre = movvar(roiBaselinePreImageStack, windowWidth); % moving window variance calculated across prestimulus 
%         fanoBaselinePre = varBaselinePre/(mean(roiBaselinePreImageStack));
        
%         % Whole sweep baseline Method1: Baseline entire trace to period before first stimulation by thresholding the fanofactor
%         % Fluorescence traces of pre-stim and post-stim periods
%         % Choose variance/fanofactor threshold
%         varThresholdPre = mean(varBaselinePre) - std(varBaselinePre); % variance threshold (mean-std) 
%         fanoThresholdPre = mean(fanoBaselinePre) - std(fanoBaselinePre); % fanofactor (variance/mean)
%         % Find variance/fanofactor beneath threshold
%         [vLowestPre, ~] = find(varBaseline < varThreshold);
%         [fanoLowestPre, ~] = find(fanoBaseline < fanoThreshold);
%         % Establish the baseline fluorescence based on all the points corresponding to the lowest moments of fano
%         roiBaselineMean = mean(roiBaselinePreImageStack(fanoLowestPre));
    
%         % Whole sweep baseline Method2: Calculate baseline value before first stimulation based on quantile bottom 10% of fanofactor 
%         % Step 1: Define quantile threshold (e.g., bottom 10%)
%         q = 0.10;  % Change to 0.05 for bottom 5%, etc.
%         quantileCutoff = quantile(fanoBaselinePre, q);
%         % Step 2: Select low points
%         [fanoLowestPre, ~] = find(fanoBaselinePre < quantileCutoff);
%         roiBaselineThisHolo = mean(roiBaselinePreImageStack(fanoLowestPre)); % baseline value to subtract
    
        % Break apart the imageStack for this trial into stim windows and rearrange according to hologram sequence 
        holoSeqThisTrial = (unique(voltMapping.outParams.sequenceThisTrial{tt}, 'stable') - min(unique(voltMapping.outParams.sequenceThisTrial{tt}, 'stable')) + 1)'; % Hologram sequence for this trial    
        for hh = 1:nHolos(voltMapping.trialCond(tt, 1))
            windowLimits = []; firstLimit = []; windowTime = []; segmentTime = []; numSegments = []; windowWidth = [];
    
            % Extract the frames associated with the current hologram, including pre and post stimulation windows
            roiFramesThisHolo = imageStack(roiXAllCells{nn}, roiYAllCells{nn}, floor((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*imagingFreq):(ceil((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*imagingFreq)+ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)));
            bkgrndFramesThisHolo = imageStack(BkgrndX, BkgrndY, floor((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*imagingFreq):(ceil((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*imagingFreq)+ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)));
            % Stim window baselining method: Same as whole sweep baselining but baseline value is set to moment of least variability right before specific target/hologram
            % Image stack in ROI before stim of this target
            roiFramesThisHoloPreStim = roiFramesThisHolo(:, :, 1:preStimWindow/1000*imagingFreq);
            bkgrndFramesThisHoloPreStim = bkgrndFramesThisHolo(:, :, 1:preStimWindow/1000*imagingFreq);
            
            roiTraceThisHoloPreStim = [];
            bkgrndTraceThisHoloPreStim = [];
            for ff = 1:size(roiFramesThisHoloPreStim, 3)
                roiTraceThisHoloPreStim(ff, 1) = mean(mean(roiFramesThisHoloPreStim(:, :, ff))); % Generate mean trace of ROI fluorescence during period before stimulation
                bkgrndTraceThisHoloPreStim(ff, 1) = mean(mean(bkgrndFramesThisHoloPreStim(:, :, ff))); % Generate mean trace of ROI fluorescence during period before stimulation
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
            
            % Background: Baselining to prestimulus at beginning of sweep: parameters for moving window variance during pre-stim period
            windowLimits = [1/imagingFreq, size(bkgrndTraceThisHoloPreStim, 1)/imagingFreq]; % in sec The range within each sweep where variance is sampled by 10 sample windows, in ms. Will leave it as 500 ms total width, starting from 110ms after start (after the step pulse).
            firstLimit = windowLimits(1)*imagingFreq; % Start of the first sampling window (in samples).
            windowTime = ceil(size(bkgrndTraceThisHoloPreStim, 1)/imagingFreq*1000)/1000;
            segmentTime = 0.010; % time length of each window (s)
            numSegments = windowTime/segmentTime; % number of sample windows within the limit
            windowWidth = floor((windowLimits(2) - windowLimits(1))/numSegments*imagingFreq); % Span of each sample window (in samples).
        
            % Background: Calculate variances/fanofactors in the prestimulus period using rolling window
            varBaselineThisHolo = movvar(bkgrndTraceThisHoloPreStim, windowWidth); % moving window variance calculated across prestimulus 
            fanoBaselineThisHolo = varBaselineThisHolo/(mean(bkgrndTraceThisHoloPreStim));
    
            % Background: Calculate baseline value from window before this hologram stimulation based on quantile bottom 10% of fanofactor 
            % Step 1: Define quantile threshold (e.g., bottom 10%)
            q = 0.10;  % Change to 0.05 for bottom 5%, etc.
            quantileCutoff = quantile(fanoBaselineThisHolo, q);
            % Step 2: Select low points
            [fanoLowestThisHolo, ~] = find(fanoBaselineThisHolo < quantileCutoff);
            bkgrndBaselineThisHolo = mean(bkgrndTraceThisHoloPreStim(fanoLowestThisHolo)); % baseline value to subtract
    
            dfThisHolo = [];
            dffThisHolo = [];
            bkgrnddfThisHolo = [];
            bkgrnddffThisHolo = [];
            dfThisHoloBkgdsub = [];
            dffThisHoloBkgdsub = [];
            for ff = 1:size(roiFramesThisHolo, 3)
                % ROI pixels for this frame
                currFrameRoi = roiFramesThisHolo(:, :, ff);
                currFrameRoiMean = nanmean(currFrameRoi, 'all'); % The mean across all pixels in the ROI for the select frame
                
                % Background pixels for this frame
                currFrameBkgrnd = bkgrndFramesThisHolo(:, :, ff);
                currFrameBkgrndMean = nanmean(currFrameBkgrnd, 'all');
                
                % Calculate df
                intensityChange = currFrameRoiMean - roiBaselineThisHolo; % essentially, intensityChange = df, and dff is intensityChnage/roiBaselineMean
                bkgrndIntensityChange = currFrameBkgrndMean - bkgrndBaselineThisHolo;
                
                dfThisHolo = [dfThisHolo; intensityChange];
                dffThisHolo = [dffThisHolo; intensityChange/roiBaselineThisHolo];
                
                bkgrnddfThisHolo = [bkgrnddfThisHolo; bkgrndIntensityChange];
                bkgrnddffThisHolo = [bkgrnddffThisHolo; bkgrndIntensityChange/bkgrndBaselineThisHolo];
                
                dfThisHoloBkgdsub = [dfThisHoloBkgdsub; intensityChange-bkgrndIntensityChange];
                dffThisHoloBkgdsub = [dffThisHoloBkgdsub; (intensityChange/roiBaselineThisHolo)-(bkgrndIntensityChange/bkgrndBaselineThisHolo)];
                
            end
            
            if UpOrDown == '2'
                dfThisHolo = -dfThisHolo;
                dffThisHolo = -dffThisHolo;
                dfThisHoloBkgdsub = -dfThisHoloBkgdsub;
                dffThisHoloBkgdsub = -dffThisHoloBkgdsub;
            end
            
            filtdffThisHolo = filter(b, a, dffThisHolo); % Run this if the voltage imaging trace need filtering
            filtdffThisHoloBkgdsub = filter(b, a, dffThisHoloBkgdsub);
            
            % Re-baseline the dff data
    %         dffThisHolo = dffThisHolo - mean(dffThisHolo(1:ceil(preStimWindow/1000*imagingFreq)));
            if ismember(tt, excludeTrials) % isnan(holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, tt))
                %           holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{hh, 1}(:, tt) = NaN(ceil(ipi*nPulses/1000*imagingFreq)+1, 1);
                analysisStruct.(holoSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(holoSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, NaN(ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)+1, 1)];
                analysisStruct.(filtHoloSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(filtHoloSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, NaN(ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)+1, 1)];
                analysisStruct.(holoSortedImagingCellNames_bkgrndsubtrct{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(holoSortedImagingCellNames_bkgrndsubtrct{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, NaN(ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)+1, 1)];
                analysisStruct.(filtHoloSortedImagingCellNames_bkgrndsubtrct{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(filtHoloSortedImagingCellNames_bkgrndsubtrct{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, NaN(ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)+1, 1)];
            
            else
                %           holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, tt), 1}(:, tt) = dfThisHolo;
                analysisStruct.(holoSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(holoSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, dffThisHolo];
                analysisStruct.(filtHoloSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(filtHoloSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, filtdffThisHolo];
                analysisStruct.(holoSortedImagingCellNames_bkgrndsubtrct{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(holoSortedImagingCellNames_bkgrndsubtrct{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, dffThisHoloBkgdsub];
                analysisStruct.(filtHoloSortedImagingCellNames_bkgrndsubtrct{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(filtHoloSortedImagingCellNames_bkgrndsubtrct{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, filtdffThisHoloBkgdsub];
            end        
        end % nHolos
    end % nCells
end % nTrials

%% Calculate mean response (and CI) for each hologram across trials and per condition
holoSortedMeanCellNames_bkgrndsubtrct = cell(nCells, 1);
filtHoloSortedMeanCellNames_bkgrndsubtrct = cell(nCells, 1);
for nn = 1:nCells
    holoSortedMeanCellNames_bkgrndsubtrct{nn} = ['holoSortedImagingMean_bkgrndsubtrct_', 'cell', num2str(nn)];
    filtHoloSortedMeanCellNames_bkgrndsubtrct{nn} = ['filtHoloSortedImagingMean_bkgrndsubtrct_', 'cell', num2str(nn)];
    analysisStruct.(holoSortedMeanCellNames_bkgrndsubtrct{nn}) = cell(nConds, 1);
    analysisStruct.(filtHoloSortedMeanCellNames_bkgrndsubtrct{nn}) = cell(nConds, 1);
    
    for cc = 1:nConds
        analysisStruct.(holoSortedMeanCellNames_bkgrndsubtrct{nn}){cc} = cell(nHolos(cc), 1);
        analysisStruct.(filtHoloSortedMeanCellNames_bkgrndsubtrct{nn}){cc} = cell(nHolos(cc), 1);
    end
end

for nn = 1:nCells
    for cc = 1:nConds
        for hh = 1:nHolos(cc)
            analysisStruct.(holoSortedMeanCellNames_bkgrndsubtrct{nn}){cc}{hh} = nanmean(analysisStruct.(holoSortedImagingCellNames_bkgrndsubtrct{nn}){cc}{hh}, 2);
            analysisStruct.(filtHoloSortedMeanCellNames_bkgrndsubtrct{nn}){cc}{hh} = nanmean(analysisStruct.(filtHoloSortedImagingCellNames_bkgrndsubtrct{nn}){cc}{hh}, 2);        
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

CIDffAllCondsCellNames_bkgrndsubtrct = cell(nCells, 1);
filtCIDffAllCondsCellNames_bkgrndsubtrct = cell(nCells, 1);
for nn = 1:nCells
    CIDffAllCondsCellNames_bkgrndsubtrct{nn} = ['CIDffAllConds_bkgrndsubtrct_', 'cell', num2str(nn)];
    filtCIDffAllCondsCellNames_bkgrndsubtrct{nn} = ['filtCIDffAllConds__bkgrndsubtrct_', 'cell', num2str(nn)];
    analysisStruct.(CIDffAllCondsCellNames_bkgrndsubtrct{nn}) = cell(nConds, 1);
    analysisStruct.(filtCIDffAllCondsCellNames_bkgrndsubtrct{nn}) = cell(nConds, 1);
end

for nn = 1:nCells
    for cc = 1:nConds 
        for hh = 1:nHolos(cc)
            confidence_level = 0.95;
            means = nanmean(analysisStruct.(holoSortedImagingCellNames_bkgrndsubtrct{nn}){cc}{hh, 1}, 2);
            filtMeans = nanmean(analysisStruct.(filtHoloSortedImagingCellNames_bkgrndsubtrct{nn}){cc}{hh, 1}, 2);       
            std_errors = std(analysisStruct.(holoSortedImagingCellNames_bkgrndsubtrct{nn}){cc}{hh, 1}, 0, 2, "omitnan") / sqrt(size(analysisStruct.(holoSortedImagingCellNames_bkgrndsubtrct{nn}){cc}{hh, 1}, 2));
            filtStd_errors = std(analysisStruct.(filtHoloSortedImagingCellNames_bkgrndsubtrct{nn}){cc}{hh, 1}, 0, 2, "omitnan") / sqrt(size(analysisStruct.(filtHoloSortedImagingCellNames_bkgrndsubtrct{nn}){cc}{hh, 1}, 2));
      
            t_score = tinv((1 + confidence_level) / 2, size(analysisStruct.(holoSortedImagingCellNames_bkgrndsubtrct{nn}){cc}{hh, 1}, 2) - 1);
            filtT_score = tinv((1 + confidence_level) / 2, size(analysisStruct.(filtHoloSortedImagingCellNames_bkgrndsubtrct{nn}){cc}{hh, 1}, 2) - 1);
            margin_of_error = t_score * std_errors;
            filtMargin_of_error = filtT_score * filtStd_errors;
            lower_bounds = means - margin_of_error;
            filtLower_bounds = filtMeans - filtMargin_of_error;
            upper_bounds = means + margin_of_error;
            filtUpper_bounds = filtMeans + filtMargin_of_error;
            if UpOrDown == '2'
                analysisStruct.(CIDffAllCondsCellNames_bkgrndsubtrct{nn}){cc}{hh, 1} = [lower_bounds, upper_bounds];
                analysisStruct.(filtCIDffAllCondsCellNames_bkgrndsubtrct{nn}){cc}{hh, 1} = [filtLower_bounds, filtUpper_bounds];
            elseif UpOrDown =='1'
                analysisStruct.(CIDffAllCondsCellNames_bkgrndsubtrct{nn}){cc}{hh, 1} = [-lower_bounds, -upper_bounds];
                analysisStruct.(filtCIDffAllCondsCellNames_bkgrndsubtrct{nn}){cc}{hh, 1} = [-filtLower_bounds, -filtUpper_bounds];
            end
        end
    end
end

voltMapping.nPulseCoordsImaging = nPulseCoordsImaging;
voltMapping.holoSeqIndex        = holoSeqIndex;
voltMapping.preStimWindow       = preStimWindow;
voltMapping.postStimWindow      = postStimWindow;

%% Align Exclusion-applied imaging with ephys traces
nn = double(input('which cell number? '));

exclFiltCIDffAllConds = voltMapping.(cellID{nn}).filtCIDffAllConds_bkgrndsubtrct;
exclFiltHoloSortedImagingMean = voltMapping.(cellID{nn}).filtHoloSortedImagingMean_bkgrndsubtrct;

for cc = 1:nConds
    for hh = 1:nHolos(cc)
        figure(cc*1000+hh);
        set(gcf, 'Position',  [100, 100, 800, 400])
        clf
 
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
        ylim([-1.5 3])
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
  