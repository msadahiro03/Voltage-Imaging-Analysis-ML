function roiMeanF = VoltImg_roiMeanPerFrameExcludeBadRows_claude(imageStack, roiX, roiY, badRowMask)
%VOLTIMG_ROIMEANPERFRAMEEXCLUDEBADROWS_CLAUDE  Vectorized mean F per frame over ROI pixels, excluding bad rows.
%
%   roiMeanF = VoltImg_roiMeanPerFrameExcludeBadRows_claude(imageStack, roiX, roiY, badRowMask)
%
%   roiX, roiY — same-length column vectors of row/column indices into imageStack(:,:,ff).
%   badRowMask — H x numFrames logical; ROI pixels whose row is flagged for that frame are set to NaN.
%
%   Drop-in replacement for VoltImg_roiMeanPerFrameExcludeBadRows. Computes the same
%   per-frame mean over the same pixel set, but vectorized via reshape + linear indexing.

[H, W, numFrames] = size(imageStack);

if isempty(roiX)
    roiMeanF = nan(numFrames, 1);
    return
end

useMask = nargin >= 4 && ~isempty(badRowMask);
if useMask && ~isequal(size(badRowMask), [H, numFrames])
    error('VoltImg_roiMeanPerFrameExcludeBadRows_claude:SizeMismatch', ...
        'badRowMask must be H x numFrames matching imageStack.');
end

lin = sub2ind([H, W], roiX(:), roiY(:));

stack2D = reshape(imageStack, H*W, numFrames);
pv = double(stack2D(lin, :));               % nROIpix x numFrames

if useMask
    rowMask = badRowMask(roiX(:), :);       % nROIpix x numFrames
    pv(rowMask) = NaN;
end

allNaN = all(isnan(pv), 1);
roiMeanF = mean(pv, 1, 'omitnan').';
roiMeanF(allNaN) = NaN;

end
