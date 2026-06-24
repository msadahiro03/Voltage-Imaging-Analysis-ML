function [deNoisedFilesPath, deNoisedFilesContents, deNoisedFilesIndex] = voltMapping_deNoiser(ImgfolderContents, imagesIndex, imagingFreq, startTime)

% Baseline
for bb = 1:floor(imagingFreq*startTime)

end
for ii = 1:length(imagesIndex)
    % Read the multi-frame image
    currImgPath = [ImgfolderContents(imagesIndex(tt)).folder, '/', ImgfolderContents(imagesIndex(tt)).name];
    info = imfinfo(currImgPath);
    numFrames = numel(info);
    
    samplingLines = [100:420];
    
    baselineImageStack = zeros(info(1).Height, info(1).Width, floor(imagingFreq*startTime));
    for bb = 1:floor(imagingFreq*startTime)
        baselineImageStack(:, :, bb) = imread(currImgPath, 'Index', bb, 'Info', info);
    end
    meanBaselineImage = mean(baselineImageStack(:, 100:420), 3);
    
    for ff = 39
        % Preallocate the image stack
        imageStack = zeros(info(1).Height, info(1).Width, 1);
        imageStack(:, :, 1) = imread(currImgPath, 'Index', ff, 'Info', info);
        figure(ff); clf; imagesc(imageStack(:, 100:420));
    end
    baseSubImage = [];
    baseSubImage = meanBaselineImage-imageStack(:, 100:420);
    figure(200); clf; imagesc(baseSubImage);
  

        for nn = 1:length(samplingLines)
            meanBaselineImage(:, samplingLines(nn))-imageStack;
        end
end
temp = mean(imageStack(:, 100:420), 2);
figure(101010); imagesc(temp);

ceil(imagingFreq*startTime);
%% Calculate time index of each stim pulse that occurs in a trial
pulseTimes = [];
holoStimTimes = [];
for hh = 1:nHolos
    for nn = 1:nPulses
        pulseTimes(:, nn) = voltMapping.outParams.firstStimTimes{1}(hh)+ (nn-1)*ipi/1000;
    end
    holoStimTimes = [holoStimTimes, pulseTimes];
end

removeFrameNums = [removeFrameNums, ceil(imagingFreq*startTime) 


%%
for tt = 1:nTrials %size(vsTest_inputs, 2)
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
        dfAllTrialsCondSorted{voltMapping.trialCond(tt, 1), 1} = [dfAllTrialsCondSorted{voltMapping.trialCond(tt, 1), 1}, -df];
        dfAllTrialsNonSorted{tt} = -df;
    elseif UpOrDown =='1'
        dfAllTrialsCondSorted{voltMapping.trialCond(tt, 1), 1} = [dfAllTrialsCondSorted{voltMapping.trialCond(tt, 1), 1}, df];
        dfAllTrialsNonSorted{tt} = df;
    end

end

voltMapping.dfAllTrialsCondSorted = dfAllTrialsCondSorted;
voltMapping.dfAllTrialsNonSorted = dfAllTrialsNonSorted;
