function [imageStackOut, stackStats] = VoltImg_applyLaserRowArtifactToStack(imageStack, gateColFirst, gateColLast, threshMode, threshParam, mcMode)
%VOLTIMG_APPLYLASERROWARTIFACTTOSTACK  Per-frame laser row detection on a 3-D stack (H x W x T).
%
%   [imageStackOut, stackStats] = VoltImg_applyLaserRowArtifactToStack( ...
%       imageStack, gateColFirst, gateColLast, threshMode, threshParam, mcMode)
%
%   mcMode:
%     'fill_for_mc' — replace flagged rows with the median of all pixels in non-flagged
%                     rows (finite values for NoRMCorre and uint16 save path).
%     'nan'         — set flagged rows to NaN (use only if normcorre and downstream
%                     aggregation support NaN; not compatible with uint16 MC TIFF).
%
%   stackStats.threshUsed(t), .nBadRows(t) optional for QC.

if nargin < 6 || isempty(mcMode)
    mcMode = 'fill_for_mc';
end

mcMode = lower(char(mcMode));
if ~ismember(mcMode, {'fill_for_mc', 'nan'})
    error('VoltImg_applyLaserRowArtifactToStack:BadMode', 'mcMode must be fill_for_mc or nan.');
end

[H, W, nKeep] = size(imageStack);
if strcmp(mcMode, 'nan')
    imageStackOut = zeros(H, W, nKeep, 'single');
else
    imageStackOut = zeros(H, W, nKeep, 'like', imageStack);
end
stackStats.threshUsed = zeros(nKeep, 1);
stackStats.nBadRows = zeros(nKeep, 1);

for ki = 1:nKeep
    fr = imageStack(:, :, ki);
    applyNan = strcmp(mcMode, 'nan');
    frDetect = fr;
    if applyNan
        frDetect = single(fr);
    end
    [frOut, badRows, ~, tu] = VoltImg_mapping_removeArtifact_v2( ...
        frDetect, gateColFirst, gateColLast, threshMode, threshParam, applyNan);
    stackStats.threshUsed(ki) = double(tu);
    stackStats.nBadRows(ki) = sum(badRows);

    if strcmp(mcMode, 'fill_for_mc') && any(badRows)
        goodMask = ~badRows;
        if ~any(goodMask)
            mv = median(fr(:), 'omitnan');
        else
            mv = median(fr(goodMask, :), 'all', 'omitnan');
        end
        if isnan(mv) || isinf(mv)
            mv = 0;
        end
        frOut = fr;
        mv = cast(mv, 'like', frOut(1));
        frOut(badRows, :) = mv;
    elseif strcmp(mcMode, 'nan')
        % already single with NaNs from v2
    else
        frOut = fr;
    end

    imageStackOut(:, :, ki) = frOut;
end

end
