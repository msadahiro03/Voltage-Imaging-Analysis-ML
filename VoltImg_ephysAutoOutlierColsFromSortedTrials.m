function [outlierColsByCond, stats] = VoltImg_ephysAutoOutlierColsFromSortedTrials(holoInc, nConds, nHolos, varargin)
%VoltImg_ephysAutoOutlierColsFromSortedTrials  Flag outlier trial columns in hologram-sorted ephys.
%
% Each element holoInc{cc}{hh} is nSamples x nTrials (one column per trial for
% that condition, in ascending global trial index among trials of that
% condition). Acquisition glitches often appear as traces far from the
% per-time median; columns are scored by RMSE to the trial-wise median trace
% and flagged with a robust median/MAD rule (MATLAB isoutlier, 'median').
%
% [outlierColsByCond, stats] = VoltImg_ephysAutoOutlierColsFromSortedTrials(holoInc, nConds, nHolos)
% ... = VoltImg_ephysAutoOutlierColsFromSortedTrials(..., 'MinCols', 6, 'MedianThresholdFactor', 4)
%
% outlierColsByCond{cc} is 1 x nCols logical. A column is true if ANY hologram
% slot for that condition flags it (CombineRule 'anyHolo').

p = inputParser;
addParameter(p, 'MinCols', 6, @(x) isnumeric(x) && isscalar(x) && x >= 3);
addParameter(p, 'MedianThresholdFactor', 4, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'CombineRule', 'anyHolo', @(s) ischar(s) || isstring(s));
parse(p, varargin{:});
minCols = p.Results.MinCols;
medianTF = p.Results.MedianThresholdFactor;
combineRule = char(p.Results.CombineRule);

outlierColsByCond = cell(nConds, 1);
stats = struct();
stats.rmseByCondHolo = cell(nConds, 1);
stats.outlierByCondHolo = cell(nConds, 1);

for cc = 1:nConds
    nH = nHolos(cc);
    nCols = [];
    for hh = 1:nH
        M = holoInc{cc}{hh, 1};
        if isempty(M)
            continue
        end
        nCols = size(M, 2);
        break
    end
    if isempty(nCols)
        outlierColsByCond{cc} = logical.empty(1, 0);
        stats.rmseByCondHolo{cc} = cell(nH, 1);
        stats.outlierByCondHolo{cc} = cell(nH, 1);
        continue
    end

    votes = false(1, nCols);
    stats.rmseByCondHolo{cc} = cell(nH, 1);
    stats.outlierByCondHolo{cc} = cell(nH, 1);

    for hh = 1:nH
        M = holoInc{cc}{hh, 1};
        if isempty(M)
            continue
        end
        if size(M, 2) ~= nCols
            error('ephysAutoOutlier: condition %d hologram %d has %d columns; expected %d.', ...
                cc, hh, size(M, 2), nCols);
        end
        nc = size(M, 2);
        if nc < minCols
            stats.rmseByCondHolo{cc}{hh} = [];
            stats.outlierByCondHolo{cc}{hh} = false(1, nc);
            continue
        end

        ref = median(M, 2, 'omitnan');
        R = M - ref;
        rmse = sqrt(mean(R .* R, 1, 'omitnan'));
        stats.rmseByCondHolo{cc}{hh} = rmse;

        tf = isoutlier(rmse, 'median', 'ThresholdFactor', medianTF);
        stats.outlierByCondHolo{cc}{hh} = tf;

        switch combineRule
            case 'anyHolo'
                votes = votes | tf;
            otherwise
                error('ephysAutoOutlier: unknown CombineRule ''%s''.', combineRule);
        end
    end
    outlierColsByCond{cc} = votes;
end

end
