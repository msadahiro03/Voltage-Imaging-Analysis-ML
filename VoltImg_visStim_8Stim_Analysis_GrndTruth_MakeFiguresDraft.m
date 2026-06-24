count = 0;
EightVisStimImagingAllTrials_Cell1_subBackground = EightVisStimImagingAllTrials_Cell1 - EightVisStimImagingAllTrials_Background;
for tt = 1:nTrials
    disp(count+tt);
    fig1 = figure(1);
    fig1.Name = ['Cell 1: trial ' num2str(tt)];
    plot(linspace(0, length(EightVisStimImagingAllTrials_Cell1(:, tt))/imagingFreq, length(EightVisStimImagingAllTrials_Cell1(:, tt))), ...
        EightVisStimImagingAllTrials_Cell1(:, tt), 'color', [0, 0.4470, 0.7410]);
    xline(0.5, 'Color', [0.5 0.5 0.5], 'LineWidth', 3);
    xlabel('Time (s)');
    ylabel('dF/F');
    ylim([-0.35 0.45]);

    fig2 = figure(2);
    fig2.Name = ['Cell 2: trial ' num2str(tt)];
    plot(linspace(0, length(EightVisStimImagingAllTrials_Cell2(:, tt))/imagingFreq, length(EightVisStimImagingAllTrials_Cell2(:, tt))), ...
        EightVisStimImagingAllTrials_Cell2(:, tt), 'color', [0.8500, 0.3250, 0.0980]);
    xline(0.5, 'Color', [0.5 0.5 0.5], 'LineWidth', 3);
    xlabel('Time (s)');
    ylabel('dF/F');
    ylim([-0.35 0.45]);

    fig3 = figure(3);
    fig3.Name = ['Cell 3: trial ' num2str(tt)];
    plot(linspace(0, length(EightVisStimImagingAllTrials_Cell3(:, tt))/imagingFreq, length(EightVisStimImagingAllTrials_Cell3(:, tt))), ...
        EightVisStimImagingAllTrials_Cell3(:, tt), 'color', [0.9290, 0.6940, 0.1250]);
    xline(0.5, 'Color', [0.5 0.5 0.5], 'LineWidth', 3);
    xlabel('Time (s)');
    ylabel('dF/F');
    ylim([-0.35 0.45]);

    fig4 = figure(4);
    fig4.Name = ['Background: trial ' num2str(tt)];
    plot(linspace(0, length(EightVisStimImagingAllTrials_Background(:, tt))/imagingFreq, length(EightVisStimImagingAllTrials_Background(:, tt))), ...
        EightVisStimImagingAllTrials_Background(:, tt), 'color', [0 0 0]);
    xline(0.5, 'Color', [0.5 0.5 0.5], 'LineWidth', 3);
    xlabel('Time (s)');
    ylabel('dF/F');
    ylim([-0.35 0.45]);

    fig5 = figure(5);
    fig5.Name = ['Cell1_subtracted: trial ' num2str(tt)];
    plot(linspace(0, length(EightVisStimImagingAllTrials_Cell1_subBackground(:, tt))/imagingFreq, length(EightVisStimImagingAllTrials_Cell1_subBackground(:, tt))), ...
        EightVisStimImagingAllTrials_Cell1_subBackground(:, tt), 'color', [0, 0.4470, 0.7410]);
    xline(0.5, 'Color', [0.5 0.5 0.5], 'LineWidth', 3);
    xlabel('Time (s)');
    ylabel('dF/F');
    ylim([-0.35 0.45]);

    fig6 = figure(6);
    fig6.Name = ['All cells: trial ' num2str(tt)];
    plot(linspace(0, length(EightVisStimImagingAllTrials_Cell1(:, tt))/imagingFreq, length(EightVisStimImagingAllTrials_Cell1(:, tt))), ...
        EightVisStimImagingAllTrials_Cell1(:, tt), 'color', [0, 0.4470, 0.7410]);
    hold on;
    plot(linspace(0, length(EightVisStimImagingAllTrials_Cell2(:, tt))/imagingFreq, length(EightVisStimImagingAllTrials_Cell2(:, tt))), ...
        EightVisStimImagingAllTrials_Cell2(:, tt), 'color', [0.8500, 0.3250, 0.0980]);
    plot(linspace(0, length(EightVisStimImagingAllTrials_Cell3(:, tt))/imagingFreq, length(EightVisStimImagingAllTrials_Cell3(:, tt))), ...
        EightVisStimImagingAllTrials_Cell3(:, tt), 'color', [0.9290, 0.6940, 0.1250]);
    xline(0.5, 'Color', [0.5 0.5 0.5], 'LineWidth', 3);
    xlabel('Time (s)');
    ylabel('dF/F');
    ylim([-0.35 0.45]);
    hold off;
    pause
end

%%
[numRows, numCols] = size(EightVisStimImagingAllTrials_Cell1);
VisStimImagingSorted_8Vis_Cell1 = cell(nStims, 1);
% Loop through the columns and place them in the correct group
for tt = 1:numCols
    VisStimImagingSorted_8Vis_Cell1{stimulus_sequence(tt)} = [EightVisStimImagingAllTrials_Cell1{stimulus_sequence(tt)} EightVisStimImagingAllTrials_Cell1(:, tt)];
end
    
% Now merge the 8 groups into 4 groups (1+5, 2+6, 3+7, 4+8), combining 2 directions into 1
numFinalGroups = 4;
VisStimImagingSorted_4Vis_Cell1 = cell(numFinalGroups, 1);
for i = 1:numFinalGroups
    % Combine groups i and i+4
    VisStimImagingSorted_4Vis_Cell1{i} = [VisStimImagingSorted_8Vis_Cell1{i}, VisStimImagingSorted_8Vis_Cell1{i+4}];
end
    
for vv = 1:length(VisStimImagingSorted_8Vis_Cell1)
    meanVisStim_8Vis(:, vv) = nanmean(VisStimImagingSorted_8Vis_Cell1{vv}, 2);
end

for vv = 1:length(VisStimImagingSorted_4Vis_Cell1)
    meanVisStim_4Vis(:, vv) = nanmean(VisStimImagingSorted_4Vis_Cell1{vv}, 2);
end
    
figure(99);
clf
plot(linspace(0, length(meanfiltVisStim_4Vis(:, 1))/imagingFreq, length(meanfiltVisStim_4Vis(:, 1))), meanfiltVisStim_4Vis(:, 1), 'LineWidth', 1.5);
hold on;
plot(linspace(0, length(meanfiltVisStim_4Vis(:, 2))/imagingFreq, length(meanfiltVisStim_4Vis(:, 2))), meanfiltVisStim_4Vis(:, 2), 'LineWidth', 1.5);
plot(linspace(0, length(meanfiltVisStim_4Vis(:, 3))/imagingFreq, length(meanfiltVisStim_4Vis(:, 3))), meanfiltVisStim_4Vis(:, 3), 'LineWidth', 1.5);
plot(linspace(0, length(meanfiltVisStim_4Vis(:, 4))/imagingFreq, length(meanfiltVisStim_4Vis(:, 4))), meanfiltVisStim_4Vis(:, 4), 'LineWidth' , 1.5);
xlabel('Time (s)');
ylabel('dF/F');
legend('0', '45', '90', '135');
hold off

%% Trials sorted by orientation with Cell1 Data
[numRows, numCols] = size(EightVisStimImagingAllTrials_Cell1);
VisStimImagingSorted_8Vis_Cell1 = cell(nStims, 1);
% Loop through the columns and place them in the correct group
for tt = 1:numCols
    VisStimImagingSorted_8Vis_Cell1{stimulus_sequence(tt)} = [VisStimImagingSorted_8Vis_Cell1{stimulus_sequence(tt)} EightVisStimImagingAllTrials_Cell1(:, tt)];
end
    
% Now merge the 8 groups into 4 groups (1+5, 2+6, 3+7, 4+8), combining 2 directions into 1
numFinalGroups = 4;
VisStimImagingSorted_4Vis_Cell1 = cell(numFinalGroups, 1);
for i = 1:numFinalGroups
    % Combine groups i and i+4
    VisStimImagingSorted_4Vis_Cell1{i} = [VisStimImagingSorted_8Vis_Cell1{i}, VisStimImagingSorted_8Vis_Cell1{i+4}];
end
    
for vv = 1:length(VisStimImagingSorted_8Vis_Cell1)
    meanVisStim_8Vis(:, vv) = nanmean(VisStimImagingSorted_8Vis_Cell1{vv}, 2);
end

for vv = 1:length(VisStimImagingSorted_4Vis_Cell1)
    meanVisStim_4Vis(:, vv) = nanmean(VisStimImagingSorted_4Vis_Cell1{vv}, 2);
end
    
figure(1000); clf
plot(linspace(0, length(meanVisStim_4Vis(:, 1))/imagingFreq, length(meanVisStim_4Vis(:, 1))), meanVisStim_4Vis(:, 1), 'LineWidth', 1);
xlabel('Time (s)');
ylabel('dF/F');
xline(0.5, 'Color', [0.5 0.5 0.5], 'LineWidth', 3);
ylim([-0.05 0.03])

figure(1045); clf
plot(linspace(0, length(meanVisStim_4Vis(:, 2))/imagingFreq, length(meanVisStim_4Vis(:, 2))), meanVisStim_4Vis(:, 2), 'LineWidth', 1);
xlabel('Time (s)');
ylabel('dF/F');
xline(0.5, 'Color', [0.5 0.5 0.5], 'LineWidth', 3);
ylim([-0.05 0.03])

figure(1090); clf
plot(linspace(0, length(meanVisStim_4Vis(:, 3))/imagingFreq, length(meanVisStim_4Vis(:, 3))), meanVisStim_4Vis(:, 3), 'LineWidth', 1);
xlabel('Time (s)');
ylabel('dF/F');
xline(0.5, 'Color', [0.5 0.5 0.5], 'LineWidth', 3);
ylim([-0.05 0.03])

figure(1135); clf
plot(linspace(0, length(meanVisStim_4Vis(:, 4))/imagingFreq, length(meanVisStim_4Vis(:, 4))), meanVisStim_4Vis(:, 4), 'LineWidth' , 1);
xlabel('Time (s)');
ylabel('dF/F');
xline(0.5, 'Color', [0.5 0.5 0.5], 'LineWidth', 3);
ylim([-0.05 0.03])

%% Trials sorted by orientation with background subtracted Cell1 Data
[numRows, numCols] = size(EightVisStimImagingAllTrials_Cell1_subBackground);
VisStimImagingSorted_8Vis_Cell1_sub = cell(nStims, 1);
% Loop through the columns and place them in the correct group
for tt = 1:numCols
    VisStimImagingSorted_8Vis_Cell1_sub{stimulus_sequence(tt)} = [VisStimImagingSorted_8Vis_Cell1_sub{stimulus_sequence(tt)} EightVisStimImagingAllTrials_Cell1_subBackground(:, tt)];
end
    
% Now merge the 8 groups into 4 groups (1+5, 2+6, 3+7, 4+8), combining 2 directions into 1
numFinalGroups = 4;
VisStimImagingSorted_4Vis_Cell1_sub = cell(numFinalGroups, 1);
for i = 1:numFinalGroups
    % Combine groups i and i+4
    VisStimImagingSorted_4Vis_Cell1_sub{i} = [VisStimImagingSorted_8Vis_Cell1_sub{i}, VisStimImagingSorted_8Vis_Cell1_sub{i+4}];
end
    
for vv = 1:length(VisStimImagingSorted_8Vis_Cell1_sub)
    meanVisStim_8Vis_sub(:, vv) = nanmean(VisStimImagingSorted_8Vis_Cell1_sub{vv}, 2);
end

for vv = 1:length(VisStimImagingSorted_4Vis_Cell1_sub)
    meanVisStim_4Vis_sub(:, vv) = nanmean(VisStimImagingSorted_4Vis_Cell1_sub{vv}, 2);
end
    
figure(2000); clf
plot(linspace(0, length(meanVisStim_4Vis_sub(:, 1))/imagingFreq, length(meanVisStim_4Vis_sub(:, 1))), meanVisStim_4Vis_sub(:, 1), 'LineWidth', 1);
xlabel('Time (s)');
ylabel('dF/F');
ylim([-0.04 0.04])
xline(0.5, 'Color', [0.5 0.5 0.5], 'LineWidth', 3);

figure(2045); clf
plot(linspace(0, length(meanVisStim_4Vis_sub(:, 2))/imagingFreq, length(meanVisStim_4Vis_sub(:, 2))), meanVisStim_4Vis_sub(:, 2), 'LineWidth', 1);
xlabel('Time (s)');
ylabel('dF/F');
ylim([-0.04 0.04])
xline(0.5, 'Color', [0.5 0.5 0.5], 'LineWidth', 3);

figure(2090); clf
plot(linspace(0, length(meanVisStim_4Vis_sub(:, 3))/imagingFreq, length(meanVisStim_4Vis_sub(:, 3))), meanVisStim_4Vis_sub(:, 3), 'LineWidth', 1);
xlabel('Time (s)');
ylabel('dF/F');
ylim([-0.04 0.04])
xline(0.5, 'Color', [0.5 0.5 0.5], 'LineWidth', 3);

figure(2135); clf
plot(linspace(0, length(meanVisStim_4Vis_sub(:, 4))/imagingFreq, length(meanVisStim_4Vis_sub(:, 4))), meanVisStim_4Vis_sub(:, 4), 'LineWidth' , 1);
xlabel('Time (s)');
ylabel('dF/F');
ylim([-0.04 0.04])
xline(0.5, 'Color', [0.5 0.5 0.5], 'LineWidth', 3);