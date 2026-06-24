%%
clear all;

%% Load files and setup
% Step 1: Read the ephys file
UpOrDown = '2';
analysisExperimentPath = char(uigetdir('/Volumes/ExData2/Voltage Imaging/VoltImg_slice_test/Analysis Results'));
analysisExperimentDir = dir(analysisExperimentPath);
analysisExperimentNames = {analysisExperimentDir.name}; % Extract names into a cell array
analysisExperimentDir = analysisExperimentDir(~startsWith(analysisExperimentNames, '.'));

% ASAP7y is 1 Jedi-2P is 2
firstPulseVoltages = cell(2, 1);
firstPulseDffs = cell(2, 1);
meanFirstPulseVoltages = cell(2, 1);
meanFirstPulseDffs = cell(2,1);

allPulsesVoltage = cell(2, 1);
allPulsesDff = cell(2, 1);

for gg = 1:size(analysisExperimentDir, 1) % For each GEVI...
    analysisFilePath = [analysisExperimentDir(gg).folder, '/', analysisExperimentDir(gg).name];
    analysisFileDir = dir(analysisFilePath);
    analysisFileNames = {analysisFileDir.name}; % Extract names into a cell array
    analysisFileDir = analysisFileDir(~startsWith(analysisFileNames, '.'));

    for nn = 1:size(analysisFileDir, 1) % For each cell...
        % Load Analysis Results
        load([analysisFileDir(nn).folder, '/', analysisFileDir(nn).name]);
    
        names = fieldnames(voltImgTest_Analysis);
        for i = 1:numel(names)
            assignin('caller', names{i}, voltImgTest_Analysis.(names{i}));
        end
        nConds = length(dfAllConds);
        
        names = fieldnames(pulseParams);
        for i = 1:numel(names)
            assignin('caller', names{i}, pulseParams.(names{i}));
        end
        
        names = fieldnames(ephysData);
        for i = 1:numel(names)
            assignin('caller', names{i}, ephysData.(names{i}));
        end
        
        % Get first pulse voltage and dff
        for vv = 1:length(dvToTest)
            firstPulseVoltages{gg}(nn, vv) = max(meanTraceAllConds(Fs*pulseStart/1000:(Fs*pulseStart/1000+Fs/pulseFreq), vv)) - mean(meanTraceAllConds(Fs*(pulseStart-50)/1000:Fs*pulseStart/1000, vv));
            firstPulseDffs{gg}(nn, vv) = max(meanFiltDffAllConds(floor(imagingFreq*pulseStart/1000):floor(imagingFreq*pulseStart/1000+imagingFreq/pulseFreq), vv))*100;
            
            for pp = 1:nPulses
                allPulsesVoltage{gg}{nn, 1}(vv, pp) = max(meanTraceAllConds(Fs*pulseStart/1000+(Fs/pulseFreq*(pp-1)):(Fs*pulseStart/1000+Fs/pulseFreq+(Fs/pulseFreq*(pp-1))), vv)) - mean(meanTraceAllConds(Fs*(pulseStart-50)/1000:Fs*pulseStart/1000, vv));
                allPulsesDff{gg}{nn, 1}(vv, pp) = max(meanFiltDffAllConds(floor(imagingFreq*pulseStart/1000+(imagingFreq/pulseFreq*(pp-1))):floor((imagingFreq*pulseStart/1000+imagingFreq/pulseFreq+(imagingFreq/pulseFreq*(pp-1)))), vv))*100;
            end
        end
        
    end
    meanFirstPulseVoltages{gg} = mean(firstPulseVoltages{gg}, 1);
    meanFirstPulseDffs{gg} = mean(firstPulseDffs{gg}, 1);
    end
%%
figure(100); clf
for gg = 1:size(analysisExperimentDir, 1)
    scatter(meanFirstPulseVoltages{gg}, meanFirstPulseDffs{gg});
    hold on
end

%compile all readings together per GEVI
compiledPulseVoltages = cell(size(analysisExperimentDir, 1), 1);
compiledPulseDffs = cell(size(analysisExperimentDir, 1), 1);
for gg = 1:size(analysisExperimentDir, 1)
    analysisFilePath = [analysisExperimentDir(gg).folder, '/', analysisExperimentDir(gg).name];
    analysisFileDir = dir(analysisFilePath);
    analysisFileNames = {analysisFileDir.name}; % Extract names into a cell array
    analysisFileDir = analysisFileDir(~startsWith(analysisFileNames, '.'));
    for nn = 1:size(analysisFileDir, 1) % For each cell...
        compiledPulseVoltages{gg} = [compiledPulseVoltages{gg}; allPulsesVoltage{gg}{nn}(1:4, :)];
        compiledPulseDffs{gg} = [compiledPulseDffs{gg}; allPulsesDff{gg}{nn}(1:4, :)];
    end
end

figure(200); clf
for gg = 1:size(analysisExperimentDir, 1)
    scatter(compiledPulseVoltages{gg}, compiledPulseDffs{gg}, 'green', 'filled');
    hold on;
    xlabel('dV from Vm (mV)')
    ylabel('dF/F')
    set(gca,'fontsize', 18, 'linewidth', 2);
    xlim([0 20]);
    ylim([0 20]);
    axis square
end
hold off

figure(201); clf
for gg = 1:size(analysisExperimentDir, 1)
    scatter(compiledPulseVoltages{gg}, compiledPulseDffs{gg}, 'green', 'filled');
    hold on;
    xlabel('dV from Vm (mV)')
    ylabel('dF/F')
    set(gca,'fontsize', 18, 'linewidth', 2);
    xlim([0 5]);
    ylim([0 5]);
    axis square
end
hold off

figure(202); clf
scatter(compiledPulseVoltages{1}, compiledPulseDffs{1}, 'red', 'filled');
hold on;
scatter(compiledPulseVoltages{2}, compiledPulseDffs{2}, 'blue', 'filled');
axis square
ylabel('dF/F');
xlabel('dV from Vm');
gca;
set(gca,'fontsize', 18, 'linewidth', 2);


figure(203); clf
scatter(compiledPulseVoltages{1}, compiledPulseDffs{1}, 'red', 'filled');
hold on;
scatter(compiledPulseVoltages{2}, compiledPulseDffs{2}, 'blue', 'filled');
axis square
xlim([0 10]);
ylim([0 4])
% ylabel('dF/F');
xlabel('dV from Vm');
gca;
set(gca,'fontsize', 18, 'linewidth', 2);



