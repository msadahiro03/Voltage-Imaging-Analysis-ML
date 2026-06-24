%% tiffStack_frameAverageVideo
% Read a multi-page (stacked) TIFF, average every nFrameAvg successive frames
% into one image, and write an MP4 (or other VideoWriter format).
%
% Large stacks: frames are read in chunks only (no full stack in memory).

%% ---- User parameters ----
tiffPath  = '/path/to/your_stack.tif';   % input multi-page TIFF
outPath   = '/path/to/output_avg10.mp4'; % output video file
nFrameAvg = 10;                         % number of raw frames per output frame
fps       = 30;                         % playback frame rate of the output video

% Optional: only use every dirStep-th directory (e.g. 2 for green channel in interleaved 2-color TIFFs)
dirStep = 1;

%% ---- Count TIFF directories ----
t = Tiff(tiffPath, 'r');
nDir = 1;
while true
    try
        t.setDirectory(nDir);
        nDir = nDir + 1;
    catch
        nDir = nDir - 1;
        break;
    end
end
if nDir < 1
    t.close();
    error('No readable directories in: %s', tiffPath);
end

dirList = 1:dirStep:nDir;
nFrames = numel(dirList);
nOut = ceil(nFrames / nFrameAvg);
fprintf('TIFF: %d directories, using %d frames (%d-step), -> %d output frames @ %g fps\n', ...
    nDir, nFrames, dirStep, nOut, fps);

%% ---- Video writer ----
[outDir, ~, ~] = fileparts(outPath);
if ~isempty(outDir) && ~isfolder(outDir)
    mkdir(outDir);
end

v = VideoWriter(outPath, 'MPEG-4');
v.FrameRate = fps;
open(v);

%% ---- First frame: size / channels ----
t.setDirectory(dirList(1));
first = t.read();
frameSize = size(first);
lastPct = -1;

for outIdx = 1:nOut
    i0 = (outIdx - 1) * nFrameAvg + 1;
    i1 = min(i0 + nFrameAvg - 1, nFrames);
    chunkLen = i1 - i0 + 1;

    acc = zeros(frameSize, 'single');
    for k = 0:(chunkLen - 1)
        t.setDirectory(dirList(i0 + k));
        acc = acc + single(t.read());
    end
    avgFrame = acc / chunkLen;

    lo = min(avgFrame(:));
    hi = max(avgFrame(:));
    u8 = uint8(255 * (single(avgFrame - lo) / single(hi - lo + eps)));
    if ndims(u8) == 2 || size(u8, 3) == 1
        vidFrame = repmat(u8, [1 1 3]);
    else
        vidFrame = u8;
    end
    writeVideo(v, vidFrame);

    pct = floor(100 * outIdx / nOut);
    if pct >= lastPct + 5 || outIdx == nOut
        fprintf('  wrote %d / %d output frames (%d%%)\n', outIdx, nOut, pct);
        lastPct = pct;
    end
end

t.close();
close(v);
fprintf('Done: %s\n', outPath);
