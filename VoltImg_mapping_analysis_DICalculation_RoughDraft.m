%% VoltMapping ConnMatrix Analysis
clearvars -except CCbsAllSessions D4ImagingAllCells nCells cellID CCwsAllSessions.D2

%% Load Analysis Results
% Run this section after loading the specific cell analysis file, then
% re-run the above sections to regenerate figures

names = fieldnames(voltMapping);
for i = 1:numel(names)
    assignin('caller', names{i}, voltMapping.(names{i}));
end
nConds = length(outParams.power);

names = fieldnames(ephys);
for i = 1:numel(names)
    assignin('caller', names{i}, ephys.(names{i}));
end

Fs = voltMapping.daqParams.Fs;
imagingFreq = voltMapping.imagingFreq;
trialTime = voltMapping.daqParams.maxSweepLengthSec;
nTrials = length(voltMapping.trialCond); % ALTERNATIVELY "length(find(cellfun(@isempty, ExpStruct.inputs)==0))". Instead of "length(ExpStruct.inputs)" this puts out true number of trials successfully recorded
powers = voltMapping.outParams.power; % ALTERNATIVELY "unique(ExpStruct.trialCond)" what powers were used
nConds = length(voltMapping.outParams.power); % ALTERNATIVELY "length(unique(ExpStruct.trialCond))" total number of powers used
nHolos = voltMapping.holoStimParams.nHolos; % number of holograms in grid
pulseDurs = unique(voltMapping.outParams.pulseDur);
nPulses = unique(voltMapping.outParams.nPulses);
ipi = voltMapping.outParams.ipi;
totalPulses = nHolos*nPulses;
% SpotCoordinates = SortedData.holoRequest.targets;
startTime = (voltMapping.holoStimParams.startTime)/1000;
imagesIndex = voltMapping.imagesIndex;
UpOrDown = voltMapping.UpOrDown;
ephysFilePath = voltMapping.ephysFilePath;
ImgsFilePath = voltMapping.ImgsFilePath;
cellID = voltMapping.cellID;

cutOffFreq = 50;   % Cutoff frequency
[b, a] = butter(4, cutOffFreq/(imagingFreq/2));

%% Pearson correlation coefficient (between sessions)
% nIter = 500; % number of times random halves will be compared. Choose from 100~500
% CCwsAllCells = cell(nCells, 1);
% for nn = 1:nCells
%     CCwsAllConds = cell(nConds, 1);
%     clear exclHoloSortedImagingAllTrials
%     exclHoloSortedImagingAllTrials = voltMapping.(cellID{nn}).exclHoloSortedImagingAllTrials;
%     for cc = 1:nConds
%         CCwsAllHolos = zeros(nHolos(cc), 1);
%         for hh = 1:nHolos(cc)
%             nTrials = size(exclHoloSortedImagingAllTrials{cc}{hh}, 2); halfTrials = floor(nTrials/2);
%             vals = nan(nIter, 1);
%             for ii = 1:nIter
%                 randIdx = randperm(nTrials);
%                 firstHalf = nanmean(exclHoloSortedImagingAllTrials{cc}{hh}(:, randIdx(1:halfTrials)), 2);
%                 secondHalf = nanmean(exclHoloSortedImagingAllTrials{cc}{hh}(:, randIdx(halfTrials+1:2*halfTrials)), 2);
%                 oddHalf = nanmean(exclHoloSortedImagingAllTrials{cc}{hh}(:, 2:2:end), 2);
%                 evenHalf = nanmean(exclHoloSortedImagingAllTrials{cc}{hh}(:, 1:2:end), 2);
%                 CCThisIter = corr(firstHalf, secondHalf);
%                 CCAllIter(ii) = max(CCThisIter, 0);
%             end
%             CCwsAllHolos(hh) = mean(CCAllIter, 'omitnan');
%         end
%         CCwsAllConds{cc} = CCwsAllHolos;
%     end
%     CCwsAllCells{nn} = CCwsAllConds;
% end

%%
nHolos = 40;

D4ImagingAllCells = cell(nCells, 1);
for nn = 1:nCells
    clear D1ImagingThisCell
    D4ImagingThisCell = voltMapping.(cellID{nn}).exclHoloSortedImagingAllTrials;
    D4ImagingAllCells{nn, 1} = D4ImagingThisCell;
end

D10ImagingAllCells = cell(nCells, 1);
for nn = 1:nCells
    clear D1ImagingThisCell
    D10ImagingThisCell = voltMapping.(cellID{nn}).exclHoloSortedImagingAllTrials;
    D10ImagingAllCells{nn, 1} = D10ImagingThisCell;
end

D4ImagingCell1 = D4ImagingAllCells{1,1}{2,1};
D10ImagingCell1 = D10ImagingAllCells{1,1}{1,1};
D4ImagingCell2 = D4ImagingAllCells{2,1}{2,1};
D10ImagingCell2 = D10ImagingAllCells{2,1}{1,1};
D4ImagingCell3 = D4ImagingAllCells{3,1}{2,1};
D10ImagingCell3 = D10ImagingAllCells{3,1}{1,1};
D4ImagingCell4 = D4ImagingAllCells{4,1}{2,1};
D10ImagingCell4 = D10ImagingAllCells{4,1}{1,1};
D4ImagingCell5 = D4ImagingAllCells{5,1}{2,1};
D10ImagingCell5 = D10ImagingAllCells{5,1}{1,1};

% Cell1
D4vD10Cell5 = [];
for hh = 1:nHolos
    nIter = 500;
    nT1 = size(D4ImagingCell5{hh}, 2);
    nT2 = size(D10ImagingCell5{hh}, 2);
    half1 = floor(nT1/2);
    half2 = floor(nT2/2);
    
    vals = nan(nIter,1);
    
    for k = 1:nIter
        % random split day 1
        idx1 = randperm(nT1);
        A1 = nanmean(D4ImagingCell5{hh}(:, idx1(1:half1)), 2);
        B1 = nanmean(D4ImagingCell5{hh}(:, idx1(half1+1:2*half1)), 2);
    
        % random split day N
        idx2 = randperm(nT2);
        A2 = nanmean(D10ImagingCell5{hh}(:, idx2(1:half2)), 2);
        B2 = nanmean(D10ImagingCell5{hh}(:, idx2(half2+1:2*half2)),2 );
    
        % correlate matching halves
        c1 = corr(A1(:), A2(:));
        c2 = corr(B1(:), B2(:));
    
        vals(k) = mean([max(c1,0), max(c2,0)]); % rectify negatives
    end
    
    D4vD10Cell5(hh, 1) = mean(vals,'omitnan');
end

CCwsD4Cell1 = CCwsAllSessions.D2{1,1}{2,1};
CCwsD4Cell2 = CCwsAllSessions.D2{2,1}{2,1};
CCwsD4Cell3 = CCwsAllSessions.D2{3,1}{2,1};
CCwsD4Cell4 = CCwsAllSessions.D2{4,1}{2,1};
CCwsD4Cell5 = CCwsAllSessions.D2{5,1}{2,1};

for hh = 1:nHolos
    RDID4vD10_Cell1(hh,1) = (CCwsD4Cell1(hh)-D4vD10Cell4(hh))/(CCwsD4Cell1(hh)+D4vD10Cell4(hh));
    RDID4vD10_Cell2(hh,1) = (CCwsD4Cell2(hh)-D4vD10Cell4(hh))/(CCwsD4Cell2(hh)+D4vD10Cell4(hh));
    RDID4vD10_Cell3(hh,1) = (CCwsD4Cell3(hh)-D4vD10Cell4(hh))/(CCwsD4Cell3(hh)+D4vD10Cell4(hh));
    RDID4vD10_Cell4(hh,1) = (CCwsD4Cell4(hh)-D4vD10Cell4(hh))/(CCwsD4Cell4(hh)+D4vD10Cell4(hh));
    RDID4vD10_Cell5(hh,1) = (CCwsD4Cell5(hh)-D4vD10Cell4(hh))/(CCwsD4Cell5(hh)+D4vD10Cell4(hh));
end

%% Plot DI distributions across multiple sessions
% Concatenate sessions of one cell into one array
RDID4vD7ThisCell = [RDID4vD10_Cell1;...
    RDID4vD10_Cell2;...
    RDID4vD10_Cell3;...
    RDID4vD10_Cell4;...
    RDID4vD10_Cell5];
labels = [ones(size(RDID4vD10_Cell1));
          2*ones(size(RDID4vD10_Cell2));
          3*ones(size(RDID4vD10_Cell3));
          4*ones(size(RDID4vD10_Cell4));
          5*ones(size(RDID4vD10_Cell5))];

% Create box plot
figure(103); 
set(gcf, 'Position',  [100, 100, 560, 420])
clf
hold on;
boxplot(RDID4vD7ThisCell, labels, 'Colors',[0 0 0], 'Symbol','');
swarmchart(labels, RDID4vD7ThisCell, 12, 'k','filled','MarkerFaceAlpha',0.3);
set(gca,'XTick',1:6,'XTickLabel',{'Cell1','Cell2','Cell3','Cell4','Cell5'});
set(gca, 'fontsize', 12);
ylabel('Drift Index');
title('Drift Across Cell and presynaptic candidates: D4 vs D10');
ylim([-1 1]); 
box on; 
% grid on;

% Optionally overlay mean ± SEM for each session
RDID4vD7ThisCell2 = [RDID4vD10_Cell1,...
    RDID4vD10_Cell2,...
    RDID4vD10_Cell3,...
    RDID4vD10_Cell4,...
    RDID4vD10_Cell5];

for ii = 1:size(RDID4vD7ThisCell2, 2)
    m = mean(RDID4vD7ThisCell2(:, ii),'omitnan');
    s = std(RDID4vD7ThisCell2(:, ii),'omitnan')/sqrt(numel(RDID4vD7ThisCell2(:, ii)));
    plot(ii, m, 'ko','MarkerFaceColor','k');
    line([ii ii], [m-s, m+s], 'Color','r','LineWidth',1.5);
end
hold off;