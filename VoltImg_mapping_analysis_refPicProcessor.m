% Prompt user to select a file from a folder on the specified path
initialPath = '/Volumes/phoenixinthesky/Masato/Voltage Imaging Data/voltMapping/voltMapping_ReferencePics';
if ~isfolder(initialPath)
    initialPath = pwd; % fallback to current folder if path doesn't exist
end

% Let user pick a file (not just a folder)
[selFile, selPath] = uigetfile({'*.tif;*.tiff;*.mat;*.*','Image or data files (*.*)'}, ...
    'Select a file from the folder', initialPath);
if isequal(selFile,0)
    error('No file selected. Operation cancelled by user.');
end

% Also create a path to the folder where the picked file is saved
% Ensure selPath is a full absolute path; uigetfile may return relative paths.
if isempty(selPath)
    selPath = pwd;
else
    selPath = char(selPath);
    if ~isfolder(selPath)
        % Try to resolve relative path
        selPath = fullfile(pwd, selPath);
        if ~isfolder(selPath)
            error('Selected file path could not be resolved: %s', selPath);
        end
    end
end
% Full path to the selected file
selectedFullPath = fullfile(selPath, selFile);

% Load the selected file into the workspace depending on its extension
[~, ~, ext] = fileparts(selFile);
switch lower(ext)
    case {'.mat'}
        % Load variables from .mat file into struct to avoid overwriting workspace unintentionally
        loadedData = load(selectedFullPath);
        % If the MAT file contains a variable named 'imageStack' and/or 'info', prefer those for downstream code
        if isfield(loadedData, 'imageStack')
            imageStack = loadedData.imageStack;
        else
            % Try to guess common image variables
            fn = fieldnames(loadedData);
            if ~isempty(fn)
                imageStack = loadedData.(fn{1});
            end
        end
        if isfield(loadedData, 'info')
            info = loadedData.info;
        end
    case {'.tif', '.tiff'}
        % Read multipage TIFF into a 3D array
        tiffInfo = imfinfo(selectedFullPath);
        numFrames = numel(tiffInfo);
        h = tiffInfo(1).Height;
        w = tiffInfo(1).Width;
        % Preallocate based on bit depth
        if isfield(tiffInfo(1), 'BitsPerSample') && tiffInfo(1).BitsPerSample == 16
            imageStack = zeros(h, w, numFrames, 'uint16');
        else
            imageStack = zeros(h, w, numFrames, 'uint8');
        end
        for k = 1:numFrames
            imageStack(:, :, k) = imread(selectedFullPath, k);
        end
        % create a minimal info struct for downstream usage if not present
        info(1).Height = h;
        info(1).Width = w;
    otherwise
        error('Unsupported file type: %s', ext);
end

% Notify user that the file has been loaded
disp(['Selected file: ', selectedFullPath]);

% Split interleaved channels: odd -> green, even -> red
nFrames = size(imageStack, 3);
if nFrames < 2
    error('Not enough frames to split interleaved channels.');
end
% Determine number of frames per channel
nGreen = ceil(nFrames/2);
nRed = floor(nFrames/2);

% Preallocate channel stacks with same class as original
greenStack = zeros(size(imageStack,1), size(imageStack,2), nGreen, class(imageStack));
redStack   = zeros(size(imageStack,1), size(imageStack,2), nRed,   class(imageStack));

% Assign frames: 1,3,5,... -> green; 2,4,6,... -> red
greenIdx = 1:2:nFrames;
redIdx   = 2:2:nFrames;
greenStack(:,:,1:numel(greenIdx)) = imageStack(:,:,greenIdx);
redStack(:,:,1:numel(redIdx))     = imageStack(:,:,redIdx);

% Replace imageStack with greenStack for downstream processing if desired,
% and keep redStack available as imageStack_red.
imageStack_red = redStack;
imageStack = greenStack;

% Set ImgsFilePath to the selected file's folder and get folder contents
ImgsFilePath = selPath;
ImgfolderContents = dir(ImgsFilePath);
disp(['Selected file: ', fullfile(selPath, selFile)]);

% Setup save directory and prepare output filename for motion-corrected images
saveDirectory = ImgsFilePath; % same folder where the selected file resides

% Build motion-corrected filename: original name with '_mc' before extension
[~, rawName, rawExt] = fileparts(selFile);
mcFilename = [rawName, '_mc', rawExt];

% Full path for saving the motion-corrected version of the selected file
mcFullPath = fullfile(saveDirectory, mcFilename);

% Display intended save path
disp(['Motion-corrected file will be saved as: ', mcFullPath]);

% Create subfolder for motion-corrected TIFFs (separate for green and red)
mcTiffFolder = fullfile(saveDirectory, 'Motion_Corrected_Tiffs');
mcTiffFolderGreen = fullfile(mcTiffFolder, 'Green');
mcTiffFolderRed   = fullfile(mcTiffFolder, 'Red');
if ~exist(mcTiffFolderGreen, 'dir'); mkdir(mcTiffFolderGreen); end
if ~exist(mcTiffFolderRed,   'dir'); mkdir(mcTiffFolderRed);   end

% Avoid hidden files and non image files (kept for compatibility but not used further here)
fileNames = {};
fileType = '.tif';
for ii = 1:length(ImgfolderContents)
    if ~ImgfolderContents(ii).isdir && ~startsWith(ImgfolderContents(ii).name, '.') && endsWith(ImgfolderContents(ii).name, fileType)
        fileNames{end+1, 1} = ImgfolderContents(ii).name; %#ok<SAGROW>
    end
end
imagesIndex = 1:numel(fileNames);

% Add NoRMCorre path if needed (adjust path to your installation)
normcorrePath = '/Users/masatosadahiro/Library/CloudStorage/GoogleDrive-masato.sadahiro@gmail.com/Other computers/Work PC/MATLAB/Scripts/NoRMCorre-master';
if exist(normcorrePath, 'dir')
    addpath(normcorrePath);
else
    warning('NoRMCorre path not found, please adjust normcorrePath variable if motion correction is required.');
end

%%
% Prepare processing for each color channel separately: imageStack (green) and imageStack_red (red)
channelNames = {'Green','Red'};
stacks = {imageStack, imageStack_red};

% Preallocate outputs containers
imageStack_mc_all = cell(1,2);
imageStackMean_all = cell(1,2);

for ch = 1:2
    Y = stacks{ch};
    if isempty(Y)
        imageStack_mc_all{ch} = [];
        imageStackMean_all{ch} = [];
        continue;
    end

    % Ensure single precision for NoRMCorre
    Y = single(Y);
    [d1,d2,nt] = size(Y);

    % Choose reasonable NoRMCorre parameters (these can be adjusted or exposed as inputs)
    % Choose NoRMCorre parameters with a heuristic for bin_width selection.
    % Heuristic: bin_width should be large enough to ensure stable
    % cross-correlation estimates but small enough to capture slow drift.
    % Use image size and number of frames to set a default, and allow
    % adjustment based on measured frame-to-frame variability.
    %
    % Rules used:
    %  - baseBin = round(max(5, min(50, sqrt(d1*d2)/50)));
    %    (depends on image area; keeps bin sizes reasonable for small/large images)
    %  - if many frames, increase bin_width so each bin contains at least ~20 frames:
    %       bin_width = max(baseBin, round(nt / max(1, floor(nt/20))));
    %  - if frame-to-frame median absolute difference is large, reduce bin_width to track faster changes.
    %
    % This provides a quantitative starting point; the user can still override
    % by changing options_rigid.bin_width after this point.

    % compute a base bin width from image area
    baseBin = round(max(5, min(50, sqrt(double(d1)*double(d2))/50)));

    % try to ensure ~20 frames per bin when possible
    targetFramesPerBin = 20;
    estBins = max(1, floor(nt / targetFramesPerBin));
    binFromFrames = max(1, round(nt / estBins));

    % measure median frame-to-frame change (robust proxy for motion magnitude)
    if nt > 1
        % sample up to first 100 frame differences to limit cost
        nSample = min(100, nt-1);
        idx = round(linspace(1, nt-1, nSample));
        diffs = zeros(1, nSample, 'like', Y);
        for kk = 1:nSample
            diffs(kk) = median(abs(Y(:,:,idx(kk)+1) - Y(:,:,idx(kk))), 'all');
        end
        medDiff = median(diffs);
    else
        medDiff = 0;
    end

    % adjust bin width: if median difference is large => decrease bin_width (track faster)
    if medDiff > 0
        % normalize medDiff relative to data dynamic range (approx)
        dynRange = double(max(Y(:)) - min(Y(:)));
        if dynRange <= 0
            dynRange = 1;
        end
        relDiff = double(medDiff) / dynRange;
    else
        relDiff = 0;
    end

    % form final bin_width: combine heuristics, clamp to [3, 200]
    bin_width = round(median([baseBin, binFromFrames, max(3, round(baseBin*(1 - 0.5*relDiff)))]));
    bin_width = max(3, min(200, bin_width)); 

    % Expose other reasonable defaults
    max_shift = 4;   % maximum allowed rigid shift (pixels)
    us_fac = 50;     % upsampling factor for subpixel registration

    options_rigid = NoRMCorreSetParms('d1',d1,'d2',d2, ...
        'bin_width',bin_width,'max_shift',max_shift,'us_fac',us_fac);

    % Informative display of chosen parameters
    disp(['NoRMCorre parameters for ', channelNames{ch}, ' channel: bin_width=', num2str(bin_width), ...
          ', max_shift=', num2str(max_shift), ', us_fac=', num2str(us_fac)]);

    % Run rigid motion correction
    try
        [M_rigid, ~] = normcorre(Y, options_rigid);
    catch ME
        warning('NoRMCorre failed for %s channel: %s. Returning uncorrected stack for this channel.', channelNames{ch}, ME.message);
        M_rigid = Y;
    end

    % Store motion-corrected stack and its mean image
    imageStack_mc_all{ch} = M_rigid;
    imageStackMean_all{ch} = mean(M_rigid, 3);

    % Save motion-corrected stack as multi-page TIFF in corresponding channel folder
    if ~isempty(M_rigid)
        % Convert to uint16 for saving (scale per-channel)
        mcMin = min(M_rigid(:));
        mcMax = max(M_rigid(:));
        if mcMax > mcMin
            M_uint16 = uint16( (M_rigid - mcMin) ./ (mcMax - mcMin) * double(intmax('uint16')) );
        else
            M_uint16 = uint16(zeros(size(M_rigid)));
        end

        outFolder = mcTiffFolderGreen;
        if ch == 2
            outFolder = mcTiffFolderRed;
        end
        % Use the original base name and append channel and _mc
        outBase = [rawName, '_', lower(channelNames{ch}), '_mc.tif'];
        outPath = fullfile(outFolder, outBase);

        % Write multipage TIFF
        for f = 1:size(M_uint16,3)
            if f == 1
                imwrite(M_uint16(:,:,f), outPath, 'Compression','none','WriteMode','overwrite');
            else
                imwrite(M_uint16(:,:,f), outPath, 'Compression','none','WriteMode','append');
            end
        end
        disp(['Saved motion-corrected ', channelNames{ch}, ' stack to: ', outPath]);
    end
end

% Assign outputs back into workspace variables expected downstream
imageStack_mc = imageStack_mc_all{1};       % motion-corrected green
imageStack_mc_red = imageStack_mc_all{2};   % motion-corrected red
imageStackMean = imageStackMean_all{1};     % mean green
imageStackMean_red = imageStackMean_all{2}; % mean red

% Save mean images (green and red) as single-page TIFFs in the folder where the original file was selected
if exist('rawPath','var') && ~isempty(rawPath)
    baseFolder = fileparts(rawPath);
elseif exist('selPath','var') && ~isempty(selPath)
    baseFolder = fileparts(selPath);
else
    baseFolder = pwd;
end

% Prepare filenames
meanGreenName = fullfile(baseFolder, [rawName, '_mean_green.tif']);
meanRedName   = fullfile(baseFolder, [rawName, '_mean_red.tif']);

% Save green mean if available
if exist('imageStackMean','var') && ~isempty(imageStackMean)
    img = imageStackMean;
    % convert to uint16 scaling across the image
    mn = min(img(:)); mx = max(img(:));
    if mx > mn
        Iout = uint16( (double(img) - double(mn)) ./ double(mx - mn) * double(intmax('uint16')) );
    else
        Iout = uint16(zeros(size(img)));
    end
    imwrite(Iout, meanGreenName, 'Compression','none','WriteMode','overwrite');
    disp(['Saved mean green image to: ', meanGreenName]);
end

% Save red mean if available
if exist('imageStackMean_red','var') && ~isempty(imageStackMean_red)
    img = imageStackMean_red;
    mn = min(img(:)); mx = max(img(:));
    if mx > mn
        Iout = uint16( (double(img) - double(mn)) ./ double(mx - mn) * double(intmax('uint16')) );
    else
        Iout = uint16(zeros(size(img)));
    end
    imwrite(Iout, meanRedName, 'Compression','none','WriteMode','overwrite');
    disp(['Saved mean red image to: ', meanRedName]);
end

% Display the two mean images (green and red) side-by-side if available
figure('Name','Mean Images (Green | Red)','NumberTitle','off');
if ~isempty(imageStackMean)
    imgG = mat2gray(imageStackMean); % normalize for display
else
    imgG = zeros(size(imageStackMean_red));
end
if ~isempty(imageStackMean_red)
    imgR = mat2gray(imageStackMean_red);
else
    imgR = zeros(size(imageStackMean));
end

% Ensure same size for montage: pad if necessary
szG = size(imgG);
szR = size(imgR);
if ~isequal(szG, szR)
    newH = max(szG(1), szR(1));
    newW = max(szG(2), szR(2));
    tmpG = zeros(newH, newW);
    tmpR = zeros(newH, newW);
    tmpG(1:szG(1),1:szG(2)) = imgG;
    tmpR(1:szR(1),1:szR(2)) = imgR;
    imgG = tmpG; imgR = tmpR;
end

% Display side-by-side
subplot(1,2,1);
imshow(imgG, []);
title('Mean Green');

subplot(1,2,2);
imshow(imgR, []);
title('Mean Red');

% Compare original green channel (before motion correction) with motion-corrected version
if exist('greenStack','var') && ~isempty(greenStack)
    origG = mean(greenStack, 3);
else
    origG = [];
end

mcG = [];
if exist('imageStack_mc','var') && ~isempty(imageStack_mc)
    mcG = mean(imageStack_mc, 3);
end

if ~isempty(origG) && ~isempty(mcG)
    % Normalize for display
    orig_disp = mat2gray(origG);
    mc_disp   = mat2gray(mcG);
    % Ensure same size
    szO = size(orig_disp); szM = size(mc_disp);
    if ~isequal(szO, szM)
        newH = max(szO(1), szM(1));
        newW = max(szO(2), szM(2));
        tmpO = zeros(newH, newW); tmpM = zeros(newH, newW);
        tmpO(1:szO(1),1:szO(2)) = orig_disp;
        tmpM(1:szM(1),1:szM(2)) = mc_disp;
        orig_disp = tmpO; mc_disp = tmpM;
    end

    % Show side-by-side with difference map
    figure('Name','Green Channel: Original | Motion-Corrected | Difference','NumberTitle','off','Position',[100 100 1200 400]);
    subplot(1,3,1);
    imshow(orig_disp, []);
    title('Original Green (mean)');

    subplot(1,3,2);
    imshow(mc_disp, []);
    title('Motion-Corrected Green (mean)');

    subplot(1,3,3);
    diffimg = imabsdiff(orig_disp, mc_disp);
    imshow(diffimg, []);
    title('Absolute Difference (normalized)');
else
    warning('Original or motion-corrected green stack not available; skipping comparison figure.');
end

% Compare original red channel (before motion correction) with motion-corrected version
if exist('redStack','var') && ~isempty(redStack)
    origR = mean(redStack, 3);
else
    origR = [];
end

mcR = [];
if exist('imageStack_mc_red','var') && ~isempty(imageStack_mc_red)
    mcR = mean(imageStack_mc_red, 3);
elseif exist('imageStack_mc_all','var') && exist('imageStack_mc','var') && ~isempty(imageStack_mc) && size(imageStack_mc_all,2) >= 2
    % fallback if variable naming differs; try to retrieve from imageStack_mc_all
    try
        mcR = mean(imageStack_mc_all{2}, 3);
    catch
        mcR = [];
    end
end

if ~isempty(origR) && ~isempty(mcR)
    % Normalize for display
    orig_disp = mat2gray(origR);
    mc_disp   = mat2gray(mcR);
    % Ensure same size
    szO = size(orig_disp); szM = size(mc_disp);
    if ~isequal(szO, szM)
        newH = max(szO(1), szM(1));
        newW = max(szO(2), szM(2));
        tmpO = zeros(newH, newW); tmpM = zeros(newH, newW);
        tmpO(1:szO(1),1:szO(2)) = orig_disp;
        tmpM(1:szM(1),1:szM(2)) = mc_disp;
        orig_disp = tmpO; mc_disp = tmpM;
    end

    % Show side-by-side with difference map
    figure('Name','Red Channel: Original | Motion-Corrected | Difference','NumberTitle','off','Position',[100 100 1200 400]);
    subplot(1,3,1);
    imshow(orig_disp, []);
    title('Original Red (mean)');

    subplot(1,3,2);
    imshow(mc_disp, []);
    title('Motion-Corrected Red (mean)');

    subplot(1,3,3);
    diffimg = imabsdiff(orig_disp, mc_disp);
    imshow(diffimg, []);
    title('Absolute Difference (normalized)');
else
    warning('Original or motion-corrected red stack not available; skipping comparison figure.');
end