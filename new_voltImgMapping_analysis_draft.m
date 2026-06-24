cleanImgStack = imageStack(:, :, frameIndex);

% Initialize parameters
nlines = size(cutTemp, 1);
lineVar = zeros(nlines, 1);
lineStdDev = zeros(nlines, 1);

% Loop through each line of frame
for ii = 1:nlines
    % Calculate statistics
    lineVar(ii) = var(data(ii, :));
    lineStdDev(ii) = std(data(ii, :));
    if lineVar(ii)>7500
        cleanImgStack(ii, :) = NaN;
    end
end


% Display results
disp('Line Statistics:');
disp('Mean   Variance   Standard Deviation   Min   Max');
for ii = 1:nlines
    fprintf('Line %d: %.2f   %.2f   %.2f   %.2f   %.2f\n', ii, means(ii), variances(ii), std_devs(ii), mins(ii), maxs(ii));
end
%%

%% Calculate df and mean df for all holos sorted by trial conditions
holoSortedImagingAllTrials = cell(nConds, 1);
for cc = 1:nConds
    holoSortedImagingAllTrials{cc} = cell(nHolos(cc), 1);
end

for tt = 1:nTrials %size(vsTest_inputs, 2)
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
        [cleanFrame] = VoltImg_mapping_removeArtifact(imageStack(:,:,frameIndex));
        imageStack(:,:,frameIndex) = cleanFrame; %replace the raw frame with cleaned up frame where artifact-corrupt lines are NaN'd
    end
    
    % Break apart the imageStack for this trial into stim windows and rearrange according to hologram sequence
    df = [];
    dff = [];
    sortingIndex = [];
    sortingIndex = 1:size(holoSeqIndex{voltMapping.trialCond(tt, 1)}, 1); 
    for hh = 1:nHolos(voltMapping.trialCond(tt, 1))
        % Extract the frames associated with the current stimulation, including pre and post stimulation windows
        framesThisHolo = imageStack(:, :, ceil((voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)-preStimWindow/1000)*imagingFreq):(ceil(voltMapping.outParams.firstStimTimes{voltMapping.trialCond(tt, 1)}(hh)*imagingFreq)+ceil((ipi*nPulses+postStimWindow)/1000*imagingFreq)));
        
        % Set baseline(f0), this is based on a single mean across the preStimWindow, so it may be a flawed approach!
        roiBaselineAllPixels = framesThisHolo(roiX, roiY, 1:ceil(preStimWindow/1000*imagingFreq));
        roiBaselineMean = nanmean(roiPixelsBaseline, 'all');
        
        dfThisHolo = [];
        dffThisHolo = [];
        for ff = 1:size(framesThisHolo, 3)
            % ROI pixels for this frame
            currFrameRoi = framesThisHolo(roiX, roiY, ff);
            currFrameRoiMean = nanmean(currFrameRoi, 'all');
            
            % Calculate df
            intensityChange = currFrameRoiMean - roiBaselineMean;
            dfThisHolo = [dfThisHolo; intensityChange];
            dffThisHolo = [dffThisHolo; intensityChange/roiBaselineMean];
        end
        
        if UpOrDown == '2'
            dffThisHolo = -dffThisHolo;
        elseif UpOrDown =='1'
            dffThisHolo = dffThisHolo;
        end
        
        % Re-normalize the dff data
        dffThisHolo = dffThisHolo - mean(dffThisHolo(1:ceil(preStimWindow/1000*imagingFreq)));
        
        if ismember(tt, excludeTrials) % isnan(holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, tt))
%           holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{hh, 1}(:, tt) = NaN(ceil(ipi*nPulses/1000*imagingFreq)+1, 1);
            holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{sortingIndex(hh), 1} = [holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{sortingIndex(hh), 1}, NaN(ceil((ipi*nPulses+(preStimWindow+postStimWindow))/1000*imagingFreq)+1, 1)];
        else
%           holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{holoSeqIndex{voltMapping.trialCond(tt, 1)}(hh, tt), 1}(:, tt) = dfThisHolo;
            holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{sortingIndex(hh), 1} = [holoSortedImagingAllTrials{voltMapping.trialCond(tt, 1)}{sortingIndex(hh), 1}, dffThisHolo];
        end
    end
end