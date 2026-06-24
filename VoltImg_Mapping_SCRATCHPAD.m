for cc = 1:nConds
    for tt = 1:size(condSortedInputs{cc}, 2)
    figure(1001); clf;
    plot(sweepThisTrial);
    hold on;
    plot(ExpStruct.outParams.trialHoloSeqIds{cc}(:, tt))
    pause
    end
end

%% Frequency analysis for imaging
% Sampling frequency
fs = 330.22;  % Sampling frequency in Hz

% Compute and plot the Power Spectral Density (PSD) using Welch's method
figure;
pwelch(data,[],[],[],fs);  % pwelch(data, window, noverlap, nfft, fs)
title('Power Spectral Density of Data');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

%% Frequency analysis for ephys
fs = 20000;  % Sampling frequency in Hz

data = thisHoloSweep;

% Compute and plot the Power Spectral Density (PSD) using Welch's method
figure;
pwelch(data,[],[],[],fs);  % pwelch(data, window, noverlap, nfft, fs)
title('Power Spectral Density of Data');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

%% Filter ephys
cutOffFreq = 30;   % Cutoff frequency
[blp, alp] = butter(4, cutOffFreq/(Fs/2), 'low');  %

dataFilt = filtfilt(blp, alp, data);
figure(101010); clf;
plot(data);
hold on
plot(dataFilt);

%%
data = mappingInputsBaselined{tt};
dataFilt = filtfilt(blp, alp, data);
figure(101010); clf;
% plot(data);
hold on
plot(dataFilt);

%%
roiAllTrials = cell(nTrials, 1);
counter = 0;
for tt = 1:nTrials %size(vsTest_inputs, 2)
    counter = counter+1;
    disp(['Trial number: ', num2str(counter)]);
    
    % Read the multi-frame image
    currImgPath = [ImgfolderContents(imagesIndex(tt)).folder, '/', ImgfolderContents(imagesIndex(tt)).name];
    info = imfinfo(currImgPath);
    numFrames = numel(info);

    % Preallocate the image stack
    imageStack = zeros(length(roiX), length(roiY), numFrames);
    
    % Read each frame and store in the stack
    for frameIndex = 1:numFrames
        thisFrame = imread(currImgPath, 'Index', frameIndex, 'Info', info);
        roiThisFrame = thisFrame(roiX, roiY);
        
%         imageStack(:, :, frameIndex) = roiThisFrame;
        meanRoiThisTrial(frameIndex, 1) = mean(roiThisFrame, 'all');
    end
    
%     roiAllTrials{tt} = imageStack;
    meanRoiAllTrials(:, tt) = meanRoiThisTrial;
end

%% Piece I wrote for downsampling (prxy decrease framerate) the mapping data
downSamplingFactor = 2;

filtHoloSortedImagingAllTrials_downsampled = cell(nConds, 1);
for cc = 1:nConds
    filtHoloSortedImagingAllTrials_downsampled{cc} = cell(nHolos(cc), 1);
    for hh = 1:nHolos(cc)
        for tt = 1:size(filtHoloSortedImagingAllTrials{cc}{hh}, 2)
            filtHoloSortedImagingAllTrials_downsampled{cc}{hh}(:, tt) = filtHoloSortedImagingAllTrials{cc}{hh}(1:downSamplingFactor:end, tt);
        end
    end
end

filtHoloSortedImagingMean_downsampled = cell(nConds, 1);
for cc = 1:nConds
    filtHoloSortedImagingMean_downsampled{cc} = cell(nHolos(cc), 1);
end

imagingFreq_downsampled_downsampled = imagingFreq/downSamplingFactor;

for cc = 1:nConds
    for hh = 1:nHolos(cc)
        filtHoloSortedImagingMean_downsampled{cc}{hh} = nanmean(filtHoloSortedImagingAllTrials_downsampled{cc}{hh}, 2);        
        % Baseline the mean holo traces
%         holoSortedImagingMean{cc}{hh} = holoSortedImagingMean{cc}{hh} - mean(holoSortedImagingMean{cc}{hh}(1:preStimWindow/1000*imagingFreq_downsampled));
%         filtHoloSortedImagingMean{cc}{hh} = filtHoloSortedImagingMean{cc}{hh} - mean(filtHoloSortedImagingMean{cc}{hh}(1:preStimWindow/1000*imagingFreq_downsampled));
      
        figure(30);
        clf
        for nn = 1:length(nPulseCoordsImaging)
            xline(nPulseCoordsImaging(nn)/imagingFreq, '--', 'LineWidth', 1.5, 'color', [.8 .8 .8]);
        end
        hold on
        plot(linspace(0, length(filtHoloSortedImagingMean_downsampled{cc}{hh})/imagingFreq_downsampled_downsampled, length(filtHoloSortedImagingMean_downsampled{cc}{hh})), filtHoloSortedImagingMean_downsampled{cc}{hh}, 'LineWidth', 1.5);
    %     axis off
        hold off
        xlabel('Time(s)')
        ylabel('df/f')
        % pause
    end
end

%% Align downsampled imaging with ephys traces
for cc = 1:nConds
    for hh = 1:nHolos(cc)
        figure(cc*1000+hh);
%         set(gcf, 'Position',  [100, 100, 600, 400])
        clf
        hold on
%         fill([linspace(0, size(CIDfAllConds{cc}{hh}, 1)/imagingFreq_downsampled, size(CIDfAllConds{cc}{hh}, 1)), fliplr(linspace(0, size(CIDfAllConds{cc}{hh}, 1)/imagingFreq_downsampled, size(CIDfAllConds{cc}{hh}, 1)))],...
%         [CIDfAllConds{cc}{hh}(:, 1)', fliplr(CIDfAllConds{cc}{hh}(:, 2)')], [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
%         % plot CI lowerbound
%         plot(linspace(0, size(CIDfAllConds{cc}{hh}, 1)/imagingFreq_downsampled, size(CIDfAllConds{cc}{hh}, 1)), CIDfAllConds{cc}{hh}(:, 1), '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]); 
%         % plot CI upperbound
%         plot(linspace(0, size(CIDfAllConds{cc}{hh}, 1)/imagingFreq_downsampled, size(CIDfAllConds{cc}{hh}, 1)), CIDfAllConds{cc}{hh}(:, 2), '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);   

        % fill([linspace(0, size(filtCIDffAllConds{cc}{hh}, 1)/imagingFreq_downsampled, size(filtCIDffAllConds{cc}{hh}, 1)), fliplr(linspace(0, size(filtCIDffAllConds{cc}{hh}, 1)/imagingFreq_downsampled, size(filtCIDffAllConds{cc}{hh}, 1)))],...
        % [filtCIDffAllConds{cc}{hh}(:, 1)'*100, fliplr(filtCIDffAllConds{cc}{hh}(:, 2)'*100)], [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
        % % plot CI lowerbound
        % plot(linspace(0, size(filtCIDffAllConds{cc}{hh}, 1)/imagingFreq_downsampled, size(filtCIDffAllConds{cc}{hh}, 1)), filtCIDffAllConds{cc}{hh}(:, 1)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]); 
        % % plot CI upperbound
        % plot(linspace(0, size(filtCIDffAllConds{cc}{hh}, 1)/imagingFreq_downsampled, size(filtCIDffAllConds{cc}{hh}, 1)), filtCIDffAllConds{cc}{hh}(:, 2)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);   
        
%         % plot ephys and voltage traces
%         [ax, hl, h2] = plotyy(linspace(0, size(holoSortedImagingMean{cc}(:, hh), 1)/imagingFreq_downsampled, size(holoSortedImagingMean{cc}(:, hh), 1)), holoSortedImagingMean{cc}(:, hh), linspace(0, size(holoSortedDataMean{cc}(:, hh), 1)/Fs, size(holoSortedDataMean{cc}(:, hh), 1)), holoSortedDataMean{cc}(:, hh));
% %         set(ax, 'XAxisLocation', 'origin');
% %         ylim([min(holoSortedImagingMean{cc}(:, hh)), max(holoSortedImagingMean{cc}(:, hh))]);
% %         ylim([-1, 2])
%         ax(2).YLim = ax(1).YLim;
%         set(ax(1), 'ytick', floor(min(holoSortedImagingMean{cc}(:, hh))):1:ceil(max(holoSortedImagingMean{cc}(:, hh))));
%         set(ax, {'ycolor'},{'g';'k'});
%         set(ax, 'FontSize', 14);
%         ylabel(ax(1), 'dF/F');
%         ylabel(ax(2), 'dV');
%         xlabel(ax(2), 'Time (s)');
%         set(hl, 'Color', 'g');
%         set(h2, 'Color', [0.7 0.7 0.7]);
%         set(hl, 'LineWidth', 2);    

        % plot voltage imaging trace only
%         yyaxis right
%         plot(linspace(0, size(holoSortedImagingMean{cc}{hh}, 1)/imagingFreq_downsampled, size(holoSortedImagingMean{cc}{hh}, 1)), holoSortedImagingMean{cc}{hh}, '-', 'linewidth', 2, 'color', 'g');
%         plot(linspace(0, size(holoSortedImagingMean{cc}{hh}, 1)/imagingFreq_downsampled, size(holoSortedImagingMean{cc}{hh}, 1)), filter(b, a, holoSortedImagingMean{cc}{hh})*100, '-', 'linewidth', 2, 'color', 'g');
        plot(linspace(0, size(filtHoloSortedImagingMean_downsampled{cc}{hh}, 1)/imagingFreq_downsampled, size(filtHoloSortedImagingMean_downsampled{cc}{hh}, 1)), filtHoloSortedImagingMean_downsampled{cc}{hh}*100, '-', 'linewidth', 2, 'color', 'g');
        ylabel('dF/F (%)');
        xlabel('Time (s)');
%         ylim([-1 max(holoSortedImagingMean{cc}{hh})*100])
%         ylim([min(holoSortedImagingMean{cc}{hh})*100 max(holoSortedImagingMean{cc}{hh})*100])
        ylim([-0.5 1.5])
        ax.YColor = [0 1 0];
        xticks([0:0.05:size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);

        % axis off
%         
        % plot ephys trace only
%         yyaxis left
%     %     axes('Position',[.70 .12 .2 .2]);
%     %     box on
%         plot(linspace(0, size(holoSortedDataMean{cc}(:, hh), 1)/Fs, size(holoSortedDataMean{cc}(:, hh), 1)), holoSortedDataMean{cc}(:, hh), 'linewidth', 1.5, 'color', [0 0 0]);
%         gca;
%         set(gca,'xtick',[], 'fontsize', 18);
% %         ylim([min(holoSortedDataMean{cc}(:, hh)), max(holoSortedDataMean{cc}(:, hh))]);
%         ylim([-0.5 2]);
%         ylabel('dV');
%         ax.YColor = [0, 0, 0];
        axis off
%         xlim([0 size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);
%         xticks([0:0.05:size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);
%         
        % show line at dff = 0
%         yline(0, '-', 'LineWidth', 1.5, 'color', [0.9 0.9 0.9]);
        for nn = 1:length(nPulseCoords)
            xline(nPulseCoords(nn)/Fs, '-', 'LineWidth', 10, 'color', [1 0 0]);
        end
        
        plot([0.01; 0.035], [0.5; 0.5], '-k', 'LineWidth', 2);
        plot([0.01; 0.01], [0.5; 1], '-k', 'LineWidth', 2);
        hold off
     pause
    end
end
