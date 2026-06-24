%% Calculate df and mean df for all holos sorted by trial conditions
% Use these lines if I need to further filter the df or dff trace
input('This step splits whole trace by holograms stimmed and calculates df/f in each stim window (ctrl+c to stop!)');
cutOffFreq = 50;   % Cutoff frequency
[b, a] = butter(4, cutOffFreq/(imagingFreq/2));  % 4th order Butterworth filter

%% Preallocate data for each cell
% A bit cumbersome but in the final step this extra work helps reorganize data for each cell analyzed
holoSortedImagingCellNames = cell(nCells, 1);
filtHoloSortedImagingCellNames = cell(nCells, 1);
for nn = 1:nCells
    F0CellNames{nn} = ['F0AllTrials_', 'cell', num2str(nn)];
    roiMeanFCellNames{nn} = ['roiMeanF_', 'cell', num2str(nn)];
    bkgrndMeanFCellNames{nn} = ['bkgrndMeanF_', 'cell', num2str(nn)];
    subScalarCellNames{nn} = ['subScalar_', 'cell', num2str(nn)];
    roiMeanFCorrectedCellNames{nn} = ['roiMeanFCorrected_', 'cell', num2str(nn)];
    globalF0CellNames{nn} = ['globalF0_', 'cell', num2str(nn)];
    dFCellNames{nn} = ['dF_', 'cell', num2str(nn)];
    dFFCellNames{nn} = ['dFF', 'cell', num2str(nn)];
    holoSortedImagingCellNames{nn} = ['holoSortedImagingAllTrials_', 'cell', num2str(nn)];
    filtHoloSortedImagingCellNames{nn} = ['filtHoloSortedImagingAllTrials_', 'cell', num2str(nn)];

    analysisStruct.(F0CellNames{nn}) = [];
    analysisStruct.(roiMeanFCellNames{nn}) = [];
    analysisStruct.(bkgrndMeanFCellNames{nn}) = [];
    analysisStruct.(subScalarCellNames{nn}) = [];
    analysisStruct.(roiMeanFCorrectedCellNames{nn}) = [];
    analysisStruct.(globalF0CellNames{nn}) = [];
    analysisStruct.(dFCellNames{nn}) = [];
    analysisStruct.(dFFCellNames{nn}) = [];
    analysisStruct.(holoSortedImagingCellNames{nn}) = cell(nConds, 1);
    analysisStruct.(filtHoloSortedImagingCellNames{nn}) = cell(nConds, 1);
    
    for cc = 1:nConds
        analysisStruct.(holoSortedImagingCellNames{nn}){cc} = cell(nHolos(cc), 1);
        analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc} = cell(nHolos(cc), 1);
    end
end

%% F, F0, dF, dF/F0 Calculation
counter = 0;
for tt = 1:nTrials %size(vsTest_inputs, 2)
    counter = counter+1;

    if ismember(tt, excludeTrials) % A skipping mechanism if trial is excluded (only relevant if ephys was done in this experiment)
        continue
    end
    
    disp(['Trial number: ', num2str(counter)]);    
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

    % Start time of first hologram stimulation in image sample
    startTimeImaging = floor(startTime*imagingFreq);    

    counter2 = 0;
    for nn = 1:nCells
        counter2 = counter2+1;
        disp(['Getting F traces and F0 for Cell: ', num2str(counter2)]);
        
        rawWholeRoiF = imageStack(roiXAllCells{nn}, roiYAllCells{nn}, :);
        roiMeanF = [];
        for ff = 1:size(rawWholeRoiF, 3)
            roiMeanF(ff, 1) = mean(mean(rawWholeRoiF(:, :, ff))); % Generate mean trace of ROI fluorescence
        end

        rawWholeBkgrndF = imageStack(bkgrndRoiXAllCells{nn}, bkgrndRoiYAllCells{nn}, :);
        bkgrndMeanF = [];
        for ff = 1:size(rawWholeBkgrndF, 3)
            bkgrndMeanF(ff, 1) = mean(mean(rawWholeBkgrndF(:, :, ff))); % Generate mean trace of np/background fluorescence 
        end        
        
        % Neuropil correction subscalar (statistically sound way of deciding)
        % This is strictly decided on prestim baseline at start of trial, and technically should involve all the 
        % prestim periods before each presynaptic candidate stimmed
        baselineIndices = 1:startTimeImaging;
        
        % Robust linear fit with intercept
        b = robustfit(bkgrndMeanF(baselineIndices), roiMeanF(baselineIndices));   % b(1)=intercept, b(2)=alpha
        subScalar = b(2);
        
        % Clamp to sane range (adjust if your optics suggest otherwise)
        subScalar = min(max(subScalar, 0), 1);
    
        % Forcibly clamp subScalar to 0.9 if above calculation decides 90%+ scalar
        % to prevent oversubtracting
        if subScalar > 0.8
            subScalar = 0.8;
        end

        subScalar = 0.7;

        roiMeanFCorrected = roiMeanF - subScalar*bkgrndMeanF;

        % F0 calculation
        % Restrict to prestim baseline period (500ms), calculate "baseline" based on fluorescence at time point of lowest fano
        windowLimits = []; firstLimit = []; windowTime = []; segmentTime = []; numSegments = []; windowWidth = [];
        baselineF = roiMeanFCorrected(1:startTimeImaging); % image stack in ROI before start of stimulation
        % Baselining to prestimulus at beginning of sweep: parameters for moving window variance during pre-stim period
        windowLimits = [1/imagingFreq, size(baselineF, 1)/imagingFreq]; % in sec The range within each sweep where variance is sampled by 10 sample windows, in ms. Will leave it as 500 ms total width, starting from 110ms after start (after the step pulse).
        firstLimit = windowLimits(1)*imagingFreq; % Start of the first sampling window (in samples).
        windowTime = ceil(size(baselineF, 1)/imagingFreq*1000)/1000; % Time of entire window (in seconds)
        segmentTime = 0.025; % time length of each window (s)
        numSegments = windowTime/segmentTime; % number of sample windows within the limit
        windowWidth = floor((windowLimits(2) - windowLimits(1))/numSegments*imagingFreq); % Span of each sample window (in samples).
        % Calculate variances/fanofactors in the prestimulus period using rolling window
        varBaseline = movvar(baselineF, windowWidth); % moving window variance calculated across prestimulus 
        fanoBaseline = varBaseline/(mean(baselineF));
        % Whole sweep baseline Method2: Calculate baseline value before first stimulation based on quantile bottom 10% of fanofactor 
        % Step 1: Define quantile threshold (e.g., bottom 10%)
        q = 0.10;  % Change to 0.05 for bottom 5%, etc.
        quantileCutoff = quantile(fanoBaseline, q);
        % Step 2: Select low points
        [fanoLowest, ~] = find(fanoBaseline < quantileCutoff);
        F0 = mean(baselineF(fanoLowest)); % baseline value to subtract

        if ismember(nn, excludeTrials) % isnan(holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, tt))
            analysisStruct.(F0CellNames{nn})(tt, 1) = NaN;
            analysisStruct.(roiMeanFCellNames{nn})(:, tt) = NaN(numFrames, 1);
            analysisStruct.(bkgrndMeanFCellNames{nn})(:, tt) = NaN(numFrames, 1);
            analysisStruct.(subScalarCellNames{nn})(tt, 1) = NaN;
            analysisStruct.(roiMeanFCorrectedCellNames{nn})(:, tt) = NaN(numFrames, 1);
        else
            analysisStruct.(F0CellNames{nn})(tt, 1) = F0;
            analysisStruct.(roiMeanFCellNames{nn})(:, tt) = roiMeanF;
            analysisStruct.(bkgrndMeanFCellNames{nn})(:, tt) = bkgrndMeanF;
            analysisStruct.(subScalarCellNames{nn})(tt, 1) = subScalar;
            analysisStruct.(roiMeanFCorrectedCellNames{nn})(:, tt) = roiMeanFCorrected;
        end
    end
end

% Global F0 calculation
counter3 = 0;
for nn = 1:nCells
    counter3 = counter3 + 1;
    disp(['Getting dF/F0 for Cell: ', num2str(counter2)]);

    globalF0ThisCell = mean(analysisStruct.(F0CellNames{nn})); % Mean but I could also use median (depending on skew)
    
    analysisStruct.(globalF0CellNames{nn}) = globalF0ThisCell;
end

% dF/F0 calculation
for nn = 1:nCells
    for tt = 1:nTrials
        FCorrected = analysisStruct.(roiMeanFCorrectedCellNames{nn})(:, tt);

        % F0 = analysisStruct.(globalF0CellNames{nn}); % Global F0
        F0 = analysisStruct.(F0CellNames{nn})(tt, 1); % Trial by trial F0
        
        dF = FCorrected - F0;
        dFF = dF/F0;
        if ismember(tt, excludeTrials) % isnan(holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, tt))
            analysisStruct.(dFCellNames{nn})(:, tt) = NaN(numFrames, 1);
            analysisStruct.(dFFCellNames{nn})(:, tt) = NaN(numFrames, 1);
        else
            analysisStruct.(dFCellNames{nn})(:, tt) = dF;
            analysisStruct.(dFFCellNames{nn})(:, tt) = dFF;
        end
    end
end

%% Break apart dFFs into stim windows and sort by hologram
counter = 0;
for nn = 1:nCells %size(vsTest_inputs, 2)
    counter = counter+1;
    disp(['Cell number: ', num2str(counter)]);

    counter2 = 0;
    for tt = 1:nTrials
        counter2 = counter2+1;
        disp(['Trial number: ', num2str(counter2)]);
        
        dFFTrace = analysisStruct.(dFFCellNames{nn})(:, tt);
        
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
            
            dFFThisHolo = dFFTrace(...
                floor((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*imagingFreq):(ceil((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*imagingFreq)+ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)));
            
            % Stim window baselining method: Same as whole sweep baselining but baseline value is set to moment of least variability right before specific target/hologram
            % Image stack in ROI before stim of this target
            dFFThisHoloPreStim = dFFThisHolo(1:preStimWindow/1000*imagingFreq);
            
            % % Baselining to prestimulus at beginning of sweep: parameters for moving window variance during pre-stim period
            % windowLimits = [1/imagingFreq, size(dFFThisHoloPreStim, 1)/imagingFreq]; % in sec The range within each sweep where variance is sampled by 10 sample windows, in ms. Will leave it as 500 ms total width, starting from 110ms after start (after the step pulse).
            % firstLimit = windowLimits(1)*imagingFreq; % Start of the first sampling window (in samples).
            % windowTime = ceil(size(dFFThisHoloPreStim, 1)/imagingFreq*1000)/1000;
            % segmentTime = 0.02; % time length of each window (s)
            % numSegments = windowTime/segmentTime; % number of sample windows within the limit
            % windowWidth = floor((windowLimits(2) - windowLimits(1))/numSegments*imagingFreq); % Span of each sample window (in samples).
            % 
            % % Calculate variances/fanofactors in the prestimulus period using rolling window
            % varBaselineThisHolo = movvar(dFFThisHoloPreStim, windowWidth); % moving window variance calculated across prestimulus 
            % fanoBaselineThisHolo = varBaselineThisHolo/(mean(dFFThisHoloPreStim));
            % 
            % % Calculate baseline value from window before this hologram stimulation based on quantile bottom 10% of fanofactor 
            % % Step 1: Define quantile threshold (e.g., bottom 10%)
            % q = 0.10;  % Change to 0.05 for bottom 5%, etc.
            % quantileCutoff = quantile(fanoBaselineThisHolo(2:end), q);
            % % Step 2: Select low points
            % [fanoLowestThisHolo, ~] = find(fanoBaselineThisHolo(2:end) < quantileCutoff);
            % dffBaselineThisHolo = mean(dFFThisHoloPreStim(fanoLowestThisHolo)); % baseline value to subtract
            
            % Baselining to just the trimmed mean of the prestim baseline period (F0)
            dFFbaselineThisHolo = trimmean(dFFThisHoloPreStim, 20); % baseline (F0) value 
            
            dFFThisHoloBaselined = dFFThisHolo - dFFbaselineThisHolo;

            if UpOrDown == '2'
                dFFThisHoloBaselined = -dFFThisHoloBaselined;
            elseif UpOrDown =='1'
                dFFThisHoloBaselined = dFFThisHoloBaselined;
            end
            
            filtdffThisHolo = filter(b, a, dFFThisHoloBaselined); % Run this if the voltage imaging trace need filtering

            if ismember(tt, excludeTrials) % isnan(holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, tt))
                %           holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{hh, 1}(:, tt) = NaN(ceil(ipi*nPulses/1000*imagingFreq)+1, 1);
                analysisStruct.(holoSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(holoSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, NaN(ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)+1, 1)];
                analysisStruct.(filtHoloSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(filtHoloSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, NaN(ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)+1, 1)];
            else
                %           holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, tt), 1}(:, tt) = dfThisHolo;
                analysisStruct.(holoSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(holoSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, dFFThisHoloBaselined];
                analysisStruct.(filtHoloSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(filtHoloSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, filtdffThisHolo];
            end        
        end % nHolos
    end % nCells
end % nTrials

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Second version below
%% Calculate df and mean df for all holos sorted by trial conditions
% Use these lines if I need to further filter the df or dff trace
input('This step splits whole trace by holograms stimmed and calculates df/f in each stim window (ctrl+c to stop!)');

% Preallocate data for each cell
% A bit cumbersome but in the final step this extra work helps reorganize data for each cell analyzed
holoSortedImagingCellNames = cell(nCells, 1);
filtHoloSortedImagingCellNames = cell(nCells, 1);
for nn = 1:nCells
    F0CellNames{nn} = ['F0AllTrials_', 'cell', num2str(nn)];
    roiMeanFCellNames{nn} = ['roiMeanF_', 'cell', num2str(nn)];
    bkgrndMeanFCellNames{nn} = ['bkgrndMeanF_', 'cell', num2str(nn)];
    subScalarCellNames{nn} = ['subScalar_', 'cell', num2str(nn)];
    roiMeanFCorrectedCellNames{nn} = ['roiMeanFCorrected_', 'cell', num2str(nn)];
    globalF0CellNames{nn} = ['globalF0_', 'cell', num2str(nn)];
    dFCellNames{nn} = ['dF_', 'cell', num2str(nn)];
    dFFCellNames{nn} = ['dFF', 'cell', num2str(nn)];
    holoSortedImagingCellNames{nn} = ['holoSortedImagingAllTrials_', 'cell', num2str(nn)];
    filtHoloSortedImagingCellNames{nn} = ['filtHoloSortedImagingAllTrials_', 'cell', num2str(nn)];

    
    analysisStruct.(roiMeanFCellNames{nn}) = [];
    analysisStruct.(bkgrndMeanFCellNames{nn}) = [];
    analysisStruct.(subScalarCellNames{nn}) = [];
    analysisStruct.(roiMeanFCorrectedCellNames{nn}) = [];
    analysisStruct.(globalF0CellNames{nn}) = [];
    analysisStruct.(dFCellNames{nn}) = [];
    analysisStruct.(dFFCellNames{nn}) = [];
    analysisStruct.(F0CellNames{nn}) = cell(nConds, 1);
    analysisStruct.(holoSortedImagingCellNames{nn}) = cell(nConds, 1);
    analysisStruct.(filtHoloSortedImagingCellNames{nn}) = cell(nConds, 1);
    
    for cc = 1:nConds
        analysisStruct.(F0CellNames{nn}){cc} = cell(nHolos(cc), 1);
        analysisStruct.(holoSortedImagingCellNames{nn}){cc} = cell(nHolos(cc), 1);
        analysisStruct.(filtHoloSortedImagingCellNames{nn}){cc} = cell(nHolos(cc), 1);
    end
end

%% F, F0, dF, dF/F0 Calculation
counter = 0;
for tt = 1:nTrials %size(vsTest_inputs, 2)
    counter = counter+1;

    % if ismember(tt, excludeTrials) % A skipping mechanism if trial is excluded (only relevant if ephys was done in this experiment)
    %     continue
    % end
    
    disp(['Trial number: ', num2str(counter)]);

    % %%%%%%%%%%%%%%%% Old Loading Method
    % % Read the multi-frame image
    % currImgPath = [ImgfolderContents(imagesIndex(tt)).folder, '/', ImgfolderContents(imagesIndex(tt)).name];
    % info = imfinfo(currImgPath);
    % numFrames = numel(info);
    % 
    % % Preallocate the image stack
    % imageStack = zeros(info(1).Height, info(1).Width, numFrames);
    % 
    % % Read each frame and store in the imageStack
    % for frameIndex = 1:numFrames
    %     imageStack(:,:,frameIndex) = imread(currImgPath, 'Index', frameIndex, 'Info', info);
    % 
    %     % % Run these lines if artifact removal is necessary, comment out if not needed
    %     % [cleanFrame] = VoltImg_mapping_removeArtifact(imageStack(:,:,frameIndex));
    %     % imageStack(:,:,frameIndex) = cleanFrame; %replace the raw frame with cleaned up frame where artifact-corrupt lines are NaN'd
    % end
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

    % Start time of first hologram stimulation in image sample
    startTimeImaging = floor(startTime*imagingFreq);    

    counter2 = 0;
    for nn = 1:nCells
        counter2 = counter2+1;
        disp(['Getting F traces and F0 for Cell: ', num2str(counter2)]);
        
        rawWholeRoiF = imageStack(roiXAllCells{nn}, roiYAllCells{nn}, :);
        roiMeanF = [];
        for ff = 1:size(rawWholeRoiF, 3)
            roiMeanF(ff, 1) = mean(mean(rawWholeRoiF(:, :, ff))); % Corrected trace this trial for this cell
        end
        
        rawWholeBkgrndF = imageStack(bkgrndRoiXAllCells{nn}, bkgrndRoiYAllCells{nn}, :);
        bkgrndMeanF = [];
        for ff = 1:size(rawWholeBkgrndF, 3)
            bkgrndMeanF(ff, 1) = mean(mean(rawWholeBkgrndF(:, :, ff))); % Generate mean trace of np/background fluorescence 
        end        
        
        % Neuropil correction subscalar (statistically sound way of deciding)
        % This is strictly decided on prestim baseline at start of trial, and technically should involve all the 
        % prestim periods before each presynaptic candidate stimmed
        baselineIndices = 1:startTimeImaging;
        
        % Robust linear fit with intercept
        b = robustfit(bkgrndMeanF(baselineIndices), roiMeanF(baselineIndices));   % b(1)=intercept, b(2)=alpha
        subScalar = b(2);
        
        % Clamp to sane range (adjust if your optics suggest otherwise)
        subScalar = min(max(subScalar, 0), 1);
    
        % Forcibly clamp subScalar to 0.9 if above calculation decides 90%+ scalar
        % to prevent oversubtracting
        if subScalar > 0.8
            subScalar = 0.8;
        end

        subScalar = 0.95;

        roiMeanFCorrected = roiMeanF - subScalar*bkgrndMeanF; % Corrected trace this trial for this cell

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Break apart the imageStack for this trial into stim windows and rearrange according to hologram sequence 
        if isempty((voltMapping.outParams.sequenceThisTrial{tt})) % Hack for 0mV trials (replacing empty holo sequence with another from a random trial)
            voltMapping.outParams.sequenceThisTrial{tt} = zeroDummySequence;
        end

        cutOffFreq = 40;   % Cutoff frequency
        [b, a] = butter(4, cutOffFreq/(imagingFreq/2));  % 4th order Butterworth filter        

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
            roiFCorrectedThisHolo = roiMeanFCorrected(...
                floor((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*imagingFreq):(ceil((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*imagingFreq)+ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)));

            % Stim window baselining method: Same as whole sweep baselining but baseline value is set to moment of least variability right before specific target/hologram
            % Image stack in ROI before stim of this target
            roiFCorrectedThisHoloPreStim = roiFCorrectedThisHolo(1:(preStimWindow/1000*imagingFreq)-1);
            
            % % Baselining to prestimulus at beginning of sweep: parameters for moving window variance during pre-stim period
            % windowLimits = [1/imagingFreq, size(dFFThisHoloPreStim, 1)/imagingFreq]; % in sec The range within each sweep where variance is sampled by 10 sample windows, in ms. Will leave it as 500 ms total width, starting from 110ms after start (after the step pulse).
            % firstLimit = windowLimits(1)*imagingFreq; % Start of the first sampling window (in samples).
            % windowTime = ceil(size(dFFThisHoloPreStim, 1)/imagingFreq*1000)/1000;
            % segmentTime = 0.02; % time length of each window (s)
            % numSegments = windowTime/segmentTime; % number of sample windows within the limit
            % windowWidth = floor((windowLimits(2) - windowLimits(1))/numSegments*imagingFreq); % Span of each sample window (in samples).
            % 
            % % Calculate variances/fanofactors in the prestimulus period using rolling window
            % varBaselineThisHolo = movvar(dFFThisHoloPreStim, windowWidth); % moving window variance calculated across prestimulus 
            % fanoBaselineThisHolo = varBaselineThisHolo/(mean(dFFThisHoloPreStim));
            % 
            % % Calculate baseline value from window before this hologram stimulation based on quantile bottom 10% of fanofactor 
            % % Step 1: Define quantile threshold (e.g., bottom 10%)
            % q = 0.10;  % Change to 0.05 for bottom 5%, etc.
            % quantileCutoff = quantile(fanoBaselineThisHolo(2:end), q);
            % % Step 2: Select low points
            % [fanoLowestThisHolo, ~] = find(fanoBaselineThisHolo(2:end) < quantileCutoff);
            % dffBaselineThisHolo = mean(dFFThisHoloPreStim(fanoLowestThisHolo)); % baseline value to subtract
            
            % Baselining to just the trimmed mean of the prestim baseline period (F0)
            f0ThisHolo = mean(roiFCorrectedThisHoloPreStim); % baseline (F0) value 
            
            dFThisHolo = roiFCorrectedThisHolo - f0ThisHolo;

            dFFThisHolo = dFThisHolo/f0ThisHolo;

            if UpOrDown == '2'
                dFFThisHolo = -dFFThisHolo;
            elseif UpOrDown =='1'
                dFFThisHolo = dFFThisHolo;
            end
            
            filtdffThisHolo = filter(b, a, dFFThisHolo); % Run this if the voltage imaging trace need filtering

            if ismember(tt, excludeTrials) % isnan(holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, tt))
                %           holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{hh, 1}(:, tt) = NaN(ceil(ipi*nPulses/1000*imagingFreq)+1, 1);
                analysisStruct.(F0CellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(F0CellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, NaN];
                analysisStruct.(holoSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(holoSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, NaN(ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)+2, 1)];
                analysisStruct.(filtHoloSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(filtHoloSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, NaN(ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)+2, 1)];
            else
                %           holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, tt), 1}(:, tt) = dfThisHolo;
                analysisStruct.(F0CellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(F0CellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, f0ThisHolo];
                analysisStruct.(holoSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(holoSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, dFFThisHolo];
                analysisStruct.(filtHoloSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(filtHoloSortedImagingCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, filtdffThisHolo];
            end        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if ismember(nn, excludeTrials) % isnan(holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, tt))
            analysisStruct.(roiMeanFCellNames{nn})(:, tt) = NaN(numFrames, 1);
            analysisStruct.(bkgrndMeanFCellNames{nn})(:, tt) = NaN(numFrames, 1);
            analysisStruct.(subScalarCellNames{nn})(tt, 1) = NaN;
            analysisStruct.(roiMeanFCorrectedCellNames{nn})(:, tt) = NaN(numFrames, 1);
        else
            analysisStruct.(roiMeanFCellNames{nn})(:, tt) = roiMeanF;
            analysisStruct.(bkgrndMeanFCellNames{nn})(:, tt) = bkgrndMeanF;
            analysisStruct.(subScalarCellNames{nn})(tt, 1) = subScalar;
            analysisStruct.(roiMeanFCorrectedCellNames{nn})(:, tt) = roiMeanFCorrected;
        end
    end
end

%% Rescaling
input('This rescales NP-subtraction. Double check what scalar method is going to be used (ctrl+c to stop!)');

% A bit cumbersome but in the final step this extra work helps reorganize data for each cell analyzed
holoSortedImagingRescaleCellNames = cell(nCells, 1);
filtHoloSortedImagingRescaleCellNames = cell(nCells, 1);
for nn = 1:nCells
    F0RescaleCellNames{nn} = ['F0RescaleAllTrials_', 'cell', num2str(nn)];
    roiMeanFRescaleCellNames{nn} = ['roiMeanFRescale_', 'cell', num2str(nn)];
    bkgrndMeanFRescaleCellNames{nn} = ['bkgrndMeanFRescale_', 'cell', num2str(nn)];
    subScalarRescaleCellNames{nn} = ['subScalarRescale_', 'cell', num2str(nn)];
    roiMeanFCorrectedRescaleCellNames{nn} = ['roiMeanFCorrectedRescale_', 'cell', num2str(nn)];
    globalF0RescaleCellNames{nn} = ['globalF0Rescale_', 'cell', num2str(nn)];
    dFRescaleCellNames{nn} = ['dFRescale_', 'cell', num2str(nn)];
    dFFRescaleCellNames{nn} = ['dFFRescale', 'cell', num2str(nn)];
    holoSortedImagingRescaleCellNames{nn} = ['holoSortedImagingAllTrialsRescale_', 'cell', num2str(nn)];
    filtHoloSortedImagingRescaleCellNames{nn} = ['filtHoloSortedImagingAllTrialsRescale_', 'cell', num2str(nn)];

    analysisStruct.(F0RescaleCellNames{nn}) = cell(nConds, 1);
    analysisStruct.(roiMeanFRescaleCellNames{nn}) = [];
    analysisStruct.(bkgrndMeanFRescaleCellNames{nn}) = [];
    analysisStruct.(subScalarRescaleCellNames{nn}) = [];
    analysisStruct.(roiMeanFCorrectedRescaleCellNames{nn}) = [];
    analysisStruct.(globalF0RescaleCellNames{nn}) = [];
    analysisStruct.(dFRescaleCellNames{nn}) = [];
    analysisStruct.(dFFRescaleCellNames{nn}) = [];
    analysisStruct.(holoSortedImagingRescaleCellNames{nn}) = cell(nConds, 1);
    analysisStruct.(filtHoloSortedImagingRescaleCellNames{nn}) = cell(nConds, 1);
    
    for cc = 1:nConds
        analysisStruct.(F0RescaleCellNames{nn}){cc} = cell(nHolos(cc), 1);
        analysisStruct.(holoSortedImagingRescaleCellNames{nn}){cc} = cell(nHolos(cc), 1);
        analysisStruct.(filtHoloSortedImagingRescaleCellNames{nn}){cc} = cell(nHolos(cc), 1);
    end
end


for nn = 1:nCells
    for tt = 1:size(analysisStruct.(roiMeanFCellNames{nn}), 2)
        roiMeanF = analysisStruct.(roiMeanFCellNames{nn})(:, tt);
        bkgrndMeanF = analysisStruct.(bkgrndMeanFCellNames{nn})(:, tt);   
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Neuropil correction subscalar (statistically sound way of deciding)
        % This is strictly decided on prestim baseline at start of trial, and technically should involve all the 
        % prestim periods before each presynaptic candidate stimmed
        baselineIndices = 1:startTimeImaging;        
        % Robust linear fit with intercept
        b = robustfit(bkgrndMeanF(baselineIndices), roiMeanF(baselineIndices));   % b(1)=intercept, b(2)=alpha
        subScalar = b(2);    
        % Clamp to sane range (adjust if your optics suggest otherwise)
        subScalar = min(max(subScalar, 0), 1);   
        % Forcibly clamp subScalar to 0.9 if above calculation decides 90%+ scalar
        % to prevent oversubtracting
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Partial override
        if subScalar > 0.8
            subScalar = 0.8;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Full override
        subScalar = 1;

        roiMeanFCorrected = roiMeanF - subScalar*bkgrndMeanF;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Break apart the imageStack for this trial into stim windows and rearrange according to hologram sequence 
        if isempty((voltMapping.outParams.sequenceThisTrial{tt})) % Hack for 0mV trials (replacing empty holo sequence with another from a random trial)
            voltMapping.outParams.sequenceThisTrial{tt} = zeroDummySequence;
        end

        cutOffFreq = 40;   % Cutoff frequency
        [b, a] = butter(4, cutOffFreq/(imagingFreq/2));  % 4th order Butterworth filter        

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
            roiFCorrectedThisHolo = roiMeanFCorrected(...
                floor((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*imagingFreq):(ceil((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*imagingFreq)+ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)));

            % Stim window baselining method: Same as whole sweep baselining but baseline value is set to moment of least variability right before specific target/hologram
            % Image stack in ROI before stim of this target
            roiFCorrectedThisHoloPreStim = roiFCorrectedThisHolo(1:(preStimWindow/1000*imagingFreq)-1);
            
            % % Baselining to prestimulus at beginning of sweep: parameters for moving window variance during pre-stim period
            % windowLimits = [1/imagingFreq, size(dFFThisHoloPreStim, 1)/imagingFreq]; % in sec The range within each sweep where variance is sampled by 10 sample windows, in ms. Will leave it as 500 ms total width, starting from 110ms after start (after the step pulse).
            % firstLimit = windowLimits(1)*imagingFreq; % Start of the first sampling window (in samples).
            % windowTime = ceil(size(dFFThisHoloPreStim, 1)/imagingFreq*1000)/1000;
            % segmentTime = 0.02; % time length of each window (s)
            % numSegments = windowTime/segmentTime; % number of sample windows within the limit
            % windowWidth = floor((windowLimits(2) - windowLimits(1))/numSegments*imagingFreq); % Span of each sample window (in samples).
            % 
            % % Calculate variances/fanofactors in the prestimulus period using rolling window
            % varBaselineThisHolo = movvar(dFFThisHoloPreStim, windowWidth); % moving window variance calculated across prestimulus 
            % fanoBaselineThisHolo = varBaselineThisHolo/(mean(dFFThisHoloPreStim));
            % 
            % % Calculate baseline value from window before this hologram stimulation based on quantile bottom 10% of fanofactor 
            % % Step 1: Define quantile threshold (e.g., bottom 10%)
            % q = 0.10;  % Change to 0.05 for bottom 5%, etc.
            % quantileCutoff = quantile(fanoBaselineThisHolo(2:end), q);
            % % Step 2: Select low points
            % [fanoLowestThisHolo, ~] = find(fanoBaselineThisHolo(2:end) < quantileCutoff);
            % dffBaselineThisHolo = mean(dFFThisHoloPreStim(fanoLowestThisHolo)); % baseline value to subtract
            
            % Baselining to just the trimmed mean of the prestim baseline period (F0)
            f0ThisHolo = mean(roiFCorrectedThisHoloPreStim); % baseline (F0) value 
            
            dFThisHolo = roiFCorrectedThisHolo - f0ThisHolo;

            dFFThisHolo = dFThisHolo/f0ThisHolo;

            if UpOrDown == '2'
                dFFThisHolo = -dFFThisHolo;
            elseif UpOrDown =='1'
                dFFThisHolo = dFFThisHolo;
            end
            
            filtdffThisHolo = filter(b, a, dFFThisHolo); % Run this if the voltage imaging trace need filtering

            if ismember(tt, excludeTrials) % isnan(holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, tt))
                %           holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{hh, 1}(:, tt) = NaN(ceil(ipi*nPulses/1000*imagingFreq)+1, 1);
                analysisStruct.(F0RescaleCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(F0RescaleCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, NaN];
                analysisStruct.(holoSortedImagingRescaleCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(holoSortedImagingRescaleCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, NaN(ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)+2, 1)];
                analysisStruct.(filtHoloSortedImagingRescaleCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(filtHoloSortedImagingRescaleCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, NaN(ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)+2, 1)];
            else
                %           holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, tt), 1}(:, tt) = dfThisHolo;
                analysisStruct.(F0RescaleCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(F0RescaleCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, f0ThisHolo];
                analysisStruct.(holoSortedImagingRescaleCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(holoSortedImagingRescaleCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, dFFThisHolo];
                analysisStruct.(filtHoloSortedImagingRescaleCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1} = [analysisStruct.(filtHoloSortedImagingRescaleCellNames{nn}){voltMapping.trialCond(tt, 1)}{holoSeqThisTrial(hh), 1}, filtdffThisHolo];
            end        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Save new version to struct
        if ismember(nn, excludeTrials) % isnan(holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, tt))
            analysisStruct.(roiMeanFRescaleCellNames{nn})(:, tt) = NaN(numFrames, 1);
            analysisStruct.(bkgrndMeanFRescaleCellNames{nn})(:, tt) = NaN(numFrames, 1);
            analysisStruct.(subScalarRescaleCellNames{nn})(tt, 1) = NaN;
            analysisStruct.(roiMeanFCorrectedRescaleCellNames{nn})(:, tt) = NaN(numFrames, 1);
        else
            analysisStruct.(roiMeanFRescaleCellNames{nn})(:, tt) = roiMeanF;
            analysisStruct.(bkgrndMeanFRescaleCellNames{nn})(:, tt) = bkgrndMeanF;
            analysisStruct.(subScalarRescaleCellNames{nn})(tt, 1) = subScalar;
            analysisStruct.(roiMeanFCorrectedRescaleCellNames{nn})(:, tt) = roiMeanFCorrected;
        end
    end
end

%% Replace with rescaled values
for nn = 1:nCells
    analysisStruct.(holoSortedImagingRescaleCellNames{nn}) = analysisStruct.(holoSortedImagingCellNames{nn});
    analysisStruct.(filtHoloSortedImagingRescaleCellNames{nn}) = analysisStruct.(filtHoloSortedImagingCellNames{nn});
end
