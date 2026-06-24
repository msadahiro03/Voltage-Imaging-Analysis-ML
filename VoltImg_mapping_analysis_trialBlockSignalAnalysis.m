%% Code for looking at signal of response across random blocks of trials
% Need to select a good hologram first. So far this version of the code will not go through all the holograms.
% therefore set cc = 1, and manually select hologram (aka 'hh').

%%
selectHolos = [3 10 13 16 24];
imagingTrialBlocks_4selectHolos = cell(length(selectHolos), 1);


% for cc = 2
    for hh = length(selectHolos)
        nTrialsThisHolo = size(exclFiltHoloSortedImagingAllTrials{cc}{selectHolos(hh)}, 2);
        
        trialsPerBlock = [nTrialsThisHolo, ...
            floor(nTrialsThisHolo/2), ...
            floor(nTrialsThisHolo/4), ...
            floor(nTrialsThisHolo/8), ...
            floor(nTrialsThisHolo/16), 1];
        imagingTrialBlocks_4selectHolos{hh} = cell(length(trialsPerBlock), 1);
        for tt = 1:length(trialsPerBlock)
            remainingTrialsThisBlock = 1:nTrialsThisHolo;
            trialIdxThisBlock = cell(floor(nTrialsThisHolo/trialsPerBlock(tt)), 1);
            for ii = 1:length(trialIdxThisBlock)
                idx = randperm(length(remainingTrialsThisBlock), trialsPerBlock(tt));
                trialIdxThisBlock{ii} = remainingTrialsThisBlock(idx);
                remainingTrialsThisBlock(idx) = [];
            end
            
            for bb = 1:length(trialIdxThisBlock)
                % if trialsPerBlock(tt) == 1
                %     imagingTrialBlocks_4selectHolos{hh}{tt, 1} = nanmean(holoSortedImagingAllTrials{cc}{selectHolos(hh)});
                % else
                imagingTrialBlocks_4selectHolos{hh}{tt, 1} = [imagingTrialBlocks_4selectHolos{hh}{tt, 1}, nanmean(exclFiltHoloSortedImagingAllTrials{cc}{selectHolos(hh)}(:, trialIdxThisBlock{bb}), 2)];
                % end
            end
        end
    end
% end  

%%
for hh = length(selectHolos)
    for tt = 1:length(trialsPerBlock)
        figure(tt)
        for bb = 1:size(imagingTrialBlocks_4selectHolos{hh}{tt}, 2)
            plot(linspace(0, size(imagingTrialBlocks_4selectHolos{hh}{tt}(:, bb), 1)/imagingFreq, size(imagingTrialBlocks_4selectHolos{hh}{tt}(:, bb), 1)), imagingTrialBlocks_4selectHolos{hh}{tt}(:, bb)*100, '-', 'linewidth', 2, 'color', 'g');
            ylabel('dF/F (%)');
            xlabel('Time (s)');
    %         ylim([-1 max(holoSortedImagingMean{cc}{hh})*100])
    %         ylim([min(holoSortedImagingMean{cc}{hh})*100 max(holoSortedImagingMean{cc}{hh})*100])
            ylim([-1.5 3])
            ax.YColor = [0 1 0];
            hold on;
        end
    end
end


    imagingTrialBlocks_4selectHolos{hh}
