function roiMeanF = VoltImg_roiMeanPerFrameExcludeBadRows(imageStack, roiX, roiY, badRowMask)
%VOLTIMG_ROIMEANPERFRAMEEXCLUDEBADROWS  Mean fluorescence per frame over ROI pixels, excluding bad rows.
%
%   roiMeanF = VoltImg_roiMeanPerFrameExcludeBadRows(imageStack, roiX, roiY, badRowMask)
%
%   roiX, roiY — same-length column vectors of row/column indices into imageStack(:,:,ff).
%   badRowMask — H x numFrames logical; ROI pixels whose row is flagged for that frame are omitted.

[H, W, numFrames] = size(imageStack);
roiMeanF = zeros(numFrames, 1);

if isempty(roiX)
    roiMeanF(:) = NaN;
    return
end

useMask = nargin >= 4 && ~isempty(badRowMask);
if useMask && ~isequal(size(badRowMask), [H, numFrames])
    error('VoltImg_roiMeanPerFrameExcludeBadRows:SizeMismatch', ...
        'badRowMask must be H x numFrames matching imageStack.');
end

lin = sub2ind([H, W], roiX(:), roiY(:));

for ff = 1:numFrames
    slice = imageStack(:, :, ff);
    pv = slice(lin);
    if useMask
        br = badRowMask(roiX(:), ff);
        pv = pv(~br);
    end
    if isempty(pv)
        roiMeanF(ff) = NaN;
    else
        roiMeanF(ff) = mean(pv, 'omitnan');
    end
end

end
