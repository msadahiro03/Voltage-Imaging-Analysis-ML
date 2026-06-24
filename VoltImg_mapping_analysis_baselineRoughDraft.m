%%
tt = 1;

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

%%
startTimeImaging = floor(startTime*imagingFreq);
postStimSeqTime = voltMapping.holoStimParams.postStimSeqTime;
postStimSeqTimeImaging = ceil(postStimSeqTime/1000*imagingFreq);

baselinePreImageStack = imageStack(roiX, roiY, 1:startTimeImaging);
baselinePostImageStack = imageStack(roiX, roiY, (size(imageStack, 3)-postStimSeqTimeImaging):size(imageStack, 3));

roiBaselinePreImageStack = [];
roiBaselinePostImageStack = [];
for ff = 1:size(baselinePreImageStack, 3)
    roiBaselinePreImageStack(ff, 1) = mean(mean(baselinePreImageStack(:, :, ff)));
end
for ff = 1:size(baselinePostImageStack, 3)
    roiBaselinePostImageStack(ff, 1) = mean(mean(baselinePostImageStack(:, :, ff)));
end


windowLimits = [1/imagingFreq, size(roiBaselineImageStack, 1)/imagingFreq]; % in sec The range within each sweep where variance is sampled by 10 sample windows, in ms. Will leave it as 500 ms total width, starting from 110ms after start (after the step pulse).
firstLimit = windowLimits(1)*imagingFreq; % Start of the first sampling window (in samples).
windowTime = ceil(size(roiBaselineImageStack, 1)/imagingFreq*1000)/1000;
segmentTime = 0.050; % time length of each window (s)
numSegments = windowTime/segmentTime; % number of sample windows within the limit
windowWidth = floor((windowLimits(2) - windowLimits(1))/numSegments*imagingFreq); % Span of each sample window (in samples). Will leave as 50 ms for now.

% Rolling window approach
varBaseline = movvar(roiBaselinePreImageStack, windowWidth);
fanoBaseline = varBaseline/(mean(roiBaselinePreImageStack));

% Find moments of var/fano below threshold
% Current threshold set at standard deviation from mean of entire period before stim begins
varThreshold = mean(varBaseline) - std(varBaseline);
fanoThreshold = mean(fanoBaseline)-std(fanoBaseline);

[vLowest, ~] = find(varBaseline < varThreshold);
[fanoLowest, ~] = find(fanoBaseline < fanoThreshold);

baselineFluor = mean(roiBaselinePreImageStack(fanoLowest));

%%

for ff = 1:size(imageStack(roiX, roiY, :), 3)
    roiImageStack(ff, 1) = mean(mean(imageStack(roiX, roiY, ff)));
end

figure(4); plot(roiImageStack);
figure(5); plot(roiImageStack - baselineFluor);

%%
figure(1); plot(roiBaselinePreImageStack);
figure(2); plot(varBaseline); hold on; yline(mean(varBaseline)); yline(mean(varBaseline)-std(varBaseline), 'color', 'red'); yline(mean(varBaseline)+std(varBaseline), 'color', 'red'); hold off;
figure(3); plot(fanoBaseline); hold on; yline(mean(fanoBaseline)); yline(mean(fanoBaseline)-std(fanoBaseline), 'color', 'red'); yline(mean(fanoBaseline)+std(fanoBaseline), 'color', 'red'); hold off;
