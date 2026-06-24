function [meanImgZ, corrImgZ, Yhp] = voltImg_volpySummaryImages(Y, fps, varargin)
%VOLTIMG_VOLPYSUMMARYIMAGES Mean and local-correlation summary images (VolPy-style).
%   [meanImgZ, corrImgZ, Yhp] = voltImg_volpySummaryImages(Y, fps)
%   Y: H x W x T (single/double). fps: imaging rate in Hz.
%
%   Mean image: temporal mean of Y, z-scored across pixels (VolPy paper).
%   Correlation image: on high-passed movie (fc = 1/3 Hz, 3rd-order Butterworth),
%   each pixel gets the mean Pearson r to its 8 in-plane neighbors (Smith/Hausser style).
%
%   Name-value pairs:
%     'fcHp' - high-pass for correlation branch (default 1/3 Hz)

p = inputParser;
addParameter(p, 'fcHp', 1/3, @(x) isnumeric(x) && isscalar(x));
parse(p, varargin{:});
fcHp = p.Results.fcHp;

Y = single(Y);
[H, W, T] = size(Y);
if T < 8
    error('voltImg_volpySummaryImages: need at least 8 frames (got %d).', T);
end

Yhp = voltImg_temporalFilterVolume(Y, fps, fcHp, 'high', 3);

mu = mean(Yhp, 3);
meanImgZ = (mu - mean(mu(:), 'omitnan')) ./ max(std(mu(:), 0, 'omitnan'), eps);

sigT = std(Yhp, 0, 3);
sigT = max(sigT, eps);
Zn = (Yhp - mean(Yhp, 3)) ./ sigT;

acc = zeros(H, W, 'single');
cnt = zeros(H, W, 'single');
for di = -1:1
    for dj = -1:1
        if di == 0 && dj == 0
            continue
        end
        iDest = max(1, 1 - di):min(H, H - di);
        jDest = max(1, 1 - dj):min(W, W - dj);
        iSrc = max(1, 1 + di):min(H, H + di);
        jSrc = max(1, 1 + dj):min(W, W + dj);
        if isempty(iDest) || isempty(jDest)
            continue
        end
        P = Zn(iDest, jDest, :);
        Q = Zn(iSrc, jSrc, :);
        rloc = mean(P .* Q, 3, 'omitnan');
        acc(iDest, jDest) = acc(iDest, jDest) + rloc;
        cnt(iDest, jDest) = cnt(iDest, jDest) + 1;
    end
end
corrImg = acc ./ max(cnt, 1);
corrImgZ = (corrImg - mean(corrImg(:), 'omitnan')) ./ max(std(corrImg(:), 0, 'omitnan'), eps);
end
