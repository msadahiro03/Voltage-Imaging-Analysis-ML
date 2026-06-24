filtCIDffAllConds_D1 = filtCIDffAllConds;
filtHoloSortedImagingMean_D1 = filtHoloSortedImagingMean;
holoSortedDataMean_D1 = holoSortedDataMean;

filtCIDffAllConds_D2 = filtCIDffAllConds;
filtHoloSortedImagingMean_D2 = filtHoloSortedImagingMean;
holoSortedDataMean_D2 = holoSortedDataMean;

filtCIDffAllConds_D3 = filtCIDffAllConds;
filtHoloSortedImagingMean_D3 = filtHoloSortedImagingMean;
holoSortedDataMean_D3 = holoSortedDataMean;

%%
%% D1 Align imaging with ephys traces
for cc = 1:nConds
    for hh = 1:nHolos(cc)
        figure(10000+1000*cc+hh);
%         set(gcf, 'Position',  [100, 100, 600, 400])
        clf
        hold on

        fill([linspace(0, size(filtCIDffAllConds_D1{cc}{hh}, 1)/imagingFreq, size(filtCIDffAllConds_D1{cc}{hh}, 1)), fliplr(linspace(0, size(filtCIDffAllConds_D1{cc}{hh}, 1)/imagingFreq, size(filtCIDffAllConds_D1{cc}{hh}, 1)))],...
        [filtCIDffAllConds_D1{cc}{hh}(:, 1)'*100, fliplr(filtCIDffAllConds_D1{cc}{hh}(:, 2)'*100)], [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
        % plot CI lowerbound
        plot(linspace(0, size(filtCIDffAllConds_D1{cc}{hh}, 1)/imagingFreq, size(filtCIDffAllConds_D1{cc}{hh}, 1)), filtCIDffAllConds_D1{cc}{hh}(:, 1)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]); 
        % plot CI upperbound
        plot(linspace(0, size(filtCIDffAllConds_D1{cc}{hh}, 1)/imagingFreq, size(filtCIDffAllConds_D1{cc}{hh}, 1)), filtCIDffAllConds_D1{cc}{hh}(:, 2)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);   
        
        plot(linspace(0, size(filtHoloSortedImagingMean_D1{cc}{hh}, 1)/imagingFreq, size(filtHoloSortedImagingMean_D1{cc}{hh}, 1)), filtHoloSortedImagingMean_D1{cc}{hh}*100, '-', 'linewidth', 2, 'color', 'g');
        ylabel('dF/F (%)');
        xlabel('Time (s)');
%         ylim([-1 max(holoSortedImagingMean{cc}{hh})*100])
%         ylim([min(holoSortedImagingMean{cc}{hh})*100 max(holoSortedImagingMean{cc}{hh})*100])
        ylim([-0.5 2])
        ax.YColor = [0 1 0];
        xticks([0:0.05:size(holoSortedDataMean_D1{cc}(:, hh), 1)/Fs]);

        axis off

        for nn = 1:length(nPulseCoords)
            xline(nPulseCoords(nn)/Fs, '-', 'LineWidth', 10, 'color', [1 0 0]);
        end
        
        plot([0.01; 0.035], [0.5; 0.5], '-k', 'LineWidth', 2);
        plot([0.01; 0.01], [0.5; 1], '-k', 'LineWidth', 2);
        hold off
     
    end
end

%% D2 Align imaging with ephys traces
for cc = 1:nConds
    for hh = 1:nHolos(cc)
        figure(20000+1000*cc+hh);
%         set(gcf, 'Position',  [100, 100, 600, 400])
        clf
        hold on

        fill([linspace(0, size(filtCIDffAllConds_D2{cc}{hh}, 1)/imagingFreq, size(filtCIDffAllConds_D2{cc}{hh}, 1)), fliplr(linspace(0, size(filtCIDffAllConds_D2{cc}{hh}, 1)/imagingFreq, size(filtCIDffAllConds_D2{cc}{hh}, 1)))],...
        [filtCIDffAllConds_D2{cc}{hh}(:, 1)'*100, fliplr(filtCIDffAllConds_D2{cc}{hh}(:, 2)'*100)], [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
        % plot CI lowerbound
        plot(linspace(0, size(filtCIDffAllConds_D2{cc}{hh}, 1)/imagingFreq, size(filtCIDffAllConds_D2{cc}{hh}, 1)), filtCIDffAllConds_D2{cc}{hh}(:, 1)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]); 
        % plot CI upperbound
        plot(linspace(0, size(filtCIDffAllConds_D2{cc}{hh}, 1)/imagingFreq, size(filtCIDffAllConds_D2{cc}{hh}, 1)), filtCIDffAllConds_D2{cc}{hh}(:, 2)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);   
        
        plot(linspace(0, size(filtHoloSortedImagingMean_D2{cc}{hh}, 1)/imagingFreq, size(filtHoloSortedImagingMean_D2{cc}{hh}, 1)), filtHoloSortedImagingMean_D2{cc}{hh}*100, '-', 'linewidth', 2, 'color', 'g');
        ylabel('dF/F (%)');
        xlabel('Time (s)');
%         ylim([-1 max(holoSortedImagingMean{cc}{hh})*100])
%         ylim([min(holoSortedImagingMean{cc}{hh})*100 max(holoSortedImagingMean{cc}{hh})*100])
        ylim([-0.5 2])
        ax.YColor = [0 1 0];
        xticks([0:0.05:size(holoSortedDataMean_D2{cc}(:, hh), 1)/Fs]);

        axis off

        for nn = 1:length(nPulseCoords)
            xline(nPulseCoords(nn)/Fs, '-', 'LineWidth', 10, 'color', [1 0 0]);
        end
        
        plot([0.01; 0.035], [0.5; 0.5], '-k', 'LineWidth', 2);
        plot([0.01; 0.01], [0.5; 1], '-k', 'LineWidth', 2);
        hold off
     
    end
end

%% D3 Align imaging with ephys traces
for cc = 1:nConds
    for hh = 1:nHolos(cc)
        figure(30000+1000*cc+hh);
%         set(gcf, 'Position',  [100, 100, 600, 400])
        clf
        hold on

        fill([linspace(0, size(filtCIDffAllConds_D3{cc}{hh}, 1)/imagingFreq, size(filtCIDffAllConds_D3{cc}{hh}, 1)), fliplr(linspace(0, size(filtCIDffAllConds_D3{cc}{hh}, 1)/imagingFreq, size(filtCIDffAllConds_D3{cc}{hh}, 1)))],...
        [filtCIDffAllConds_D3{cc}{hh}(:, 1)'*100, fliplr(filtCIDffAllConds_D3{cc}{hh}(:, 2)'*100)], [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
        % plot CI lowerbound
        plot(linspace(0, size(filtCIDffAllConds_D3{cc}{hh}, 1)/imagingFreq, size(filtCIDffAllConds_D3{cc}{hh}, 1)), filtCIDffAllConds_D3{cc}{hh}(:, 1)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]); 
        % plot CI upperbound
        plot(linspace(0, size(filtCIDffAllConds_D3{cc}{hh}, 1)/imagingFreq, size(filtCIDffAllConds_D3{cc}{hh}, 1)), filtCIDffAllConds_D3{cc}{hh}(:, 2)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);   
        
        plot(linspace(0, size(filtHoloSortedImagingMean_D3{cc}{hh}, 1)/imagingFreq, size(filtHoloSortedImagingMean_D3{cc}{hh}, 1)), filtHoloSortedImagingMean_D3{cc}{hh}*100, '-', 'linewidth', 2, 'color', 'g');
        ylabel('dF/F (%)');
        xlabel('Time (s)');
%         ylim([-1 max(holoSortedImagingMean{cc}{hh})*100])
%         ylim([min(holoSortedImagingMean{cc}{hh})*100 max(holoSortedImagingMean{cc}{hh})*100])
        ylim([-0.5 2])
        ax.YColor = [0 1 0];
        xticks([0:0.05:size(holoSortedDataMean_D3{cc}(:, hh), 1)/Fs]);

        axis off

        for nn = 1:length(nPulseCoords)
            xline(nPulseCoords(nn)/Fs, '-', 'LineWidth', 10, 'color', [1 0 0]);
        end
        
        plot([0.01; 0.035], [0.5; 0.5], '-k', 'LineWidth', 2);
        plot([0.01; 0.01], [0.5; 1], '-k', 'LineWidth', 2);
        hold off
     
    end
end