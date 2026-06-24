function [maxDvTrialMask, maxDvFrameCap] = VoltImg_mapping_maxDvStackSamplingPlan(nImgTrials, imagingFreq, eligibleTrialTT)
%VOLTIMG_MAPPING_MAXDVSTACKSAMPLINGPLAN  Plan maxDvStack reference mean: 50% random trials + 4s frame cap.
%
%   [maxDvTrialMask, maxDvFrameCap] = VoltImg_mapping_maxDvStackSamplingPlan(nImgTrials, imagingFreq, eligibleTrialTT)
%
%   maxDvFrameCap = max(1, floor(imagingFreq*4)) — only first cap frames of each MC stack
%   contribute to the per-trial mean plane used for maxDvStack (full MC movies are still saved).
%
%   eligibleTrialTT — row or column vector of trial indices tt in 1..nImgTrials that are allowed
%   to contribute (e.g. non-ephys-excluded). A uniform random percentage subset of these is selected
%   (at least one trial if eligibleTrialTT is non-empty).

maxDvFrameCap = max(1, floor(double(imagingFreq) * 4));
maxDvTrialMask = false(1, nImgTrials);
eligibleTrialTT = eligibleTrialTT(:).';
eligibleTrialTT = eligibleTrialTT(eligibleTrialTT >= 1 & eligibleTrialTT <= nImgTrials);
eligibleTrialTT = unique(eligibleTrialTT, 'stable');
nEl = numel(eligibleTrialTT);
if nEl == 0
    return;
end

pickPercentage = 100/100;

nPick = max(1, ceil(nEl * pickPercentage));
picked = eligibleTrialTT(randperm(nEl, nPick));
maxDvTrialMask(picked) = true;

end
