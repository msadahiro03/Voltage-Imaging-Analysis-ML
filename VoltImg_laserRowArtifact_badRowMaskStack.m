function badRowMask = VoltImg_laserRowArtifact_badRowMaskStack(imageStack, gateColFirst, gateColLast, threshMode, threshParam)
%VOLTIMG_LASERROWARTIFACT_BADROWMASKSTACK  H x T logical: rows flagged by same rule as removeArtifact_v2.
%
%   badRowMask = VoltImg_laserRowArtifact_badRowMaskStack(imageStack, gateColFirst, ...
%       gateColLast, threshMode, threshParam)
%
%   imageStack may be uint16/single/double; detection uses single per frame.

[H, ~, nKeep] = size(imageStack);
badRowMask = false(H, nKeep);

for ki = 1:nKeep
    fr = single(imageStack(:, :, ki));
    [~, badRows, ~, ~] = VoltImg_mapping_removeArtifact_v2(fr, gateColFirst, gateColLast, ...
        threshMode, threshParam, false);
    badRowMask(:, ki) = badRows;
end

end
