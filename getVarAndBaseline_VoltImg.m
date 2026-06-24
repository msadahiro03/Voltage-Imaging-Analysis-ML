% This code is based on Hillel's "smart_zero" code
% This code will take a select span from each sweep/trial and break it down into sampling windows, 
% then analyzes the variance in each sampling window. 
% The window with the lowest variance then becomes the new baseline value that gets subtracted from the whole trace.
function [ TracesZeroedSorted ] = getVarAndBaseline_VoltImg(TracesSorted, imagingFreq, roiBaselineImageStack)

windowLimits = [1/imagingFreq, size(roiBaselineImageStack, 1)/imagingFreq]; % in sec The range within each sweep where variance is sampled by 10 sample windows, in ms. Will leave it as 500 ms total width, starting from 110ms after start (after the step pulse).
firstLimit = windowLimits(1)*imagingFreq; % Start of the first sampling window (in samples).
windowTime = ceil(size(roiBaselineImageStack, 1)/imagingFreq*1000)/1000;
segmentTime = 0.100; % time length of each window (s)
numSegments = windowTime/segmentTime; % number of sample windows within the limit
windowWidth = floor((windowLimits(2) - windowLimits(1))/numSegments*imagingFreq); % Span of each sample window (in samples). Will leave as 50 ms for now.

% Rolling window approach
v = movvar(roiBaselineImageStack, windowWidth);


% Break baseline into 10 "windowWidth" segments and analyze the variance for each window
varVector = zeros(1,numSegments); % Vector to collect variances from each sampling window in a select trial
for ss = 1:numSegments
    varVector(ss) = var(roiBaselineImageStack((firstLimit+(ss-1)*windowWidth):(firstLimit+(ss*windowWidth)))); % excerpt of trace corresponding to the segment
end
[~, minVarSegNum] = min(varVector);

% Find segment number that has the lowest variance
minVarSegNum = cell(visStimTypes, 1);
for vv = 1:visStimTypes
    minVarSegPerTrial = [];
    minVarSegAllTrials = zeros(1, size(varVectorSorted{vv}, 2));
    for ii = 1:size(varVectorSorted{vv}, 2)
        [~, minVarSegPerTrial] = min(varVectorSorted{vv, 1}{1, ii});
        minVarSegAllTrials(1, ii) = minVarSegPerTrial; 
    end
    minVarSegNum{vv} = minVarSegAllTrials;
end        

% Compute offset as the mean value of the least variable segment
% Subtract offset from trial to zero the trial trace
baselineOffsetSorted = cell(visStimTypes, 1);
TracesZeroedSorted = cell(visStimTypes, 1);
for vv = 1:visStimTypes
    baselineOffsetAllTrials = [];
    TraceZeroedAllTrials = cell(1, size(TracesSorted{vv}, 2));
    temp = [];
    for ii = 1:size(TracesSorted{vv}, 2)
        if minVarSegNum{vv, 1}(1, ii) == 17
            temp = TracesSorted{vv}(1:1000, ii);
        else
            temp = TracesSorted{vv}(round(firstLimit+(minVarSegNum{vv, 1}(1, ii)-1)*windowWidth+1):round((firstLimit+(minVarSegNum{vv, 1}(1, ii)*windowWidth))), ii);
        end
        baselineOffsetAllTrials(1, ii) = mean(temp); % offset calculated by mean of trace segment
        TraceZeroedAllTrials{1, ii} = TracesSorted{vv}(:, ii) - baselineOffsetAllTrials(1, ii);
    end
    TraceZeroedAllTrials = horzcat(TraceZeroedAllTrials{:}); % un-nest one level
    baselineOffsetSorted{vv} = baselineOffsetAllTrials;
    TracesZeroedSorted{vv, 1} = TraceZeroedAllTrials; % Main output, traces across all trials with respective baselines subtracted
end    

% % zero by setting max or min of baseline to zero
% if mean(thisTrace<0)
%     zeroed_trace=thisTrace-max(thisTrace(Exp_Defaults.Fs*0.15:Exp_Defaults.Fs*0.5));
% else
%     zeroed_trace=thisTrace-min(thisTrace(Exp_Defaults.Fs*0.15:Exp_Defaults.Fs*0.5));
% end
end