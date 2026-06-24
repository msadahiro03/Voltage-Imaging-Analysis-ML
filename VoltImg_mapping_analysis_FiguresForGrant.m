%% Align Exclusion-applied imaging with ephys traces
nn = double(input('which cell number? '));

exclCIDffAllConds = voltMapping.(cellID{nn}).exclCIDffAllConds;
exclHoloSortedImagingMean = voltMapping.(cellID{nn}).exclHoloSortedImagingMean;
exclFiltCIDffAllConds = voltMapping.(cellID{nn}).exclFiltCIDffAllConds;
exclFiltHoloSortedImagingMean = voltMapping.(cellID{nn}).exclFiltHoloSortedImagingMean;

for cc = 1:nConds
    for hh = 1:nHolos(cc)
        figure(cc*1000+hh);
        set(gcf, 'Position',  [100, 100, 500, 400])
        clf

        fill([linspace(0, size(exclFiltCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclFiltCIDffAllConds{cc}{hh}, 1)), fliplr(linspace(0, size(exclFiltCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclFiltCIDffAllConds{cc}{hh}, 1)))],...
            [exclFiltCIDffAllConds{cc}{hh}(:, 1)'*100, fliplr(exclFiltCIDffAllConds{cc}{hh}(:, 2)'*100)], [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
        hold on;
        % plot CI lowerbound
        plot(linspace(0, size(exclFiltCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclFiltCIDffAllConds{cc}{hh}, 1)), exclFiltCIDffAllConds{cc}{hh}(:, 1)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);

        % plot CI upperbound
        plot(linspace(0, size(exclFiltCIDffAllConds{cc}{hh}, 1)/imagingFreq, size(exclFiltCIDffAllConds{cc}{hh}, 1)), exclFiltCIDffAllConds{cc}{hh}(:, 2)*100, '--', 'linewidth', 1, 'color', [0.7 0.7 0.7]);

        %         % plot ephys and voltage traces
        %         [ax, hl, h2] = plotyy(linspace(0, size(holoSortedImagingMean{cc}(:, hh), 1)/imagingFreq, size(holoSortedImagingMean{cc}(:, hh), 1)), holoSortedImagingMean{cc}(:, hh), linspace(0, size(holoSortedDataMean{cc}(:, hh), 1)/Fs, size(holoSortedDataMean{cc}(:, hh), 1)), holoSortedDataMean{cc}(:, hh));
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
        % yyaxis right
        plot(linspace(0, size(exclFiltHoloSortedImagingMean{cc}{hh}, 1)/imagingFreq, size(exclFiltHoloSortedImagingMean{cc}{hh}, 1)), exclFiltHoloSortedImagingMean{cc}{hh}*100, '-', 'linewidth', 2, 'color', '#00A14B');
        % ylabel('dF/F (%)');
        xlabel('Time (s)');
        ylim([-1.5 3.5])
        xlim([0.01, 0.16])
        ax.YColor = [0 1 0];
        set(gca, 'fontsize', 16);
        set(gca, 'fontsize', 16);
        % xticks([0:0.05:size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);

        % plot ephys trace only
        %         yyaxis left
        %     %     axes('Position',[.70 .12 .2 .2]);
        %     %     box on
        %         plot(linspace(0, size(holoSortedDataMean{cc}(:, hh), 1)/Fs, size(holoSortedDataMean{cc}(:, hh), 1)), holoSortedDataMean{cc}(:, hh), 'linewidth', 1.5, 'color', [0 0 0]);
        % gca;
        % set(gca,'xtick',[], 'fontsize', 18);
        % %         ylim([min(holoSortedDataMean{cc}(:, hh)), max(holoSortedDataMean{cc}(:, hh))]);
        %         ylim([-0.5 2]);
        %         ylabel('dV');
        %         ax.YColor = [0, 0, 0];
        % axis off
        %         xlim([0 size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);
        %         xticks([0:0.05:size(holoSortedDataMean{cc}(:, hh), 1)/Fs]);

        for nn = 1:length(nPulseCoords)
            xline(nPulseCoords(nn)/Fs, '-', 'LineWidth', 2, 'color', [1 0 0 0]);
        end
        plot([0.01; 0.035], [-1.5; -1.5], '-k', 'LineWidth', 2);
        plot([0.01; 0.01], [-1.5; -1], '-k', 'LineWidth', 2);
        hold off
        box off
        axis off
        pause
    end
end
   