function [cleanFrame, badRows, lineVar, threshUsed] = VoltImg_mapping_removeArtifact_v2(currFrame, gateColFirst, gateColLast, threshMode, threshParam, applyNan)
%VOLTIMG_MAPPING_REMOVEARTIFACT_V2  Row-wise laser streak detection from line variance in a column band.
%
%   [cleanFrame, badRows, lineVar, threshUsed] = ...
%       VoltImg_mapping_removeArtifact_v2(currFrame, gateColFirst, gateColLast, ...
%       threshMode, threshParam, applyNan)
%
%   For each row, variance is computed across columns gateColFirst:gateColLast
%   (inclusive). Rows exceeding a threshold are flagged.
%
%   threshMode:
%     'fixed'      — threshParam is the variance cutoff (scalar).
%     'mad'        — threshParam is k; threshold = median(lineVar) + k*mad(lineVar,1).
%     'percentile' — threshParam is P; threshold = prctile(lineVar, P); bad if lineVar > T.
%
%   applyNan — if true, flagged rows in cleanFrame are set to NaN; otherwise cleanFrame
%              is a copy of currFrame (use badRows for external fill).

if nargin < 6
    applyNan = false;
end

[H, W] = size(currFrame);
g1 = max(1, min(gateColFirst, gateColLast));
g2 = min(W, max(gateColFirst, gateColLast));
if g2 < g1
    error('VoltImg_mapping_removeArtifact_v2:InvalidGate', 'Invalid gate column range.');
end

% var/median/mad/prctile require floating point in MATLAB.
statArea = single(currFrame(:, g1:g2));
lineVar = var(statArea, 0, 2);

switch lower(char(threshMode))
    case 'fixed'
        threshUsed = threshParam;
        badRows = lineVar > threshUsed;
    case 'mad'
        k = threshParam;
        medLv = median(lineVar, 'omitnan');
        madLv = mad(lineVar, 1);
        threshUsed = medLv + k * madLv;
        badRows = lineVar > threshUsed;
    case 'percentile'
        P = threshParam;
        threshUsed = prctile(lineVar, P);
        badRows = lineVar > threshUsed;
    otherwise
        error('VoltImg_mapping_removeArtifact_v2:BadMode', 'threshMode must be fixed, mad, or percentile.');
end

cleanFrame = currFrame;
if applyNan
    if ~(isa(cleanFrame, 'single') || isa(cleanFrame, 'double'))
        cleanFrame = single(cleanFrame);
    end
    cleanFrame(badRows, :) = NaN;
end

badRows = logical(badRows);

end
