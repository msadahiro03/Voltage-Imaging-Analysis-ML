
% This code is for subsampling all the trials for all presynaptic targets
% for one condition that must be manually selected. Then, the second
% section runs x number if iterations of random subsample blocks to
% determine detection rate (aka test power of given number of trials). This
% can statistically justify how many trials I need in any given experiment
% to claim I sampled enough.
% The variable 'testData' is used to designate the condition.

%%
nValidTrials = [];
for hh = 1:nHolos(2)
    testData = [];
    testData = voltMapping.Cell1.exclFiltHoloSortedImagingAllTrials_commonF0{2, 1}{hh, 1};
    nValidTrials(hh, 1) = sum(~isnan(testData(1, :)));
end

subSampPermutations = 4;
figure(5000);
clf
set(gcf, 'Position',  [100, 100, 1300, 1000])
hold on;
for pp = 1:subSampPermutations
    subSampConnMatrix = [];
    
    for hh = 1:nHolos(2)
        testData = [];
        testData = voltMapping.Cell1.exclFiltHoloSortedImagingAllTrials_commonF0{2, 1}{hh, 1};
        validData = [];
        invalidTrials = length(find(isnan(testData(1, :))));
        validTrials = length(find(~isnan(testData(1, :))));
        
        % Remove excluded (nan'd) trials
        for tt = 1:size(testData, 2)
            if ~isnan(testData(1, tt))
                validData = [validData, testData(:, tt)];
            end
        end
        
        % subSampFactors = [1, 0.75, 0.5, 0.25, 0.1]; 
        % subSampIndex = floor(validTrials*subSampFactors);
        % subSampIndex = [validTrials, 800, 700, 600, 500, 400, 300, 200];
        % subSampIndexTop = floor(min(nValidTrials)/100)*100;
        % subSampIndexTop = validTrials-100+50;
        subSampIndexTop = floor(min(nValidTrials)/75)*75;
        subSampIndex = [validTrials, subSampIndexTop : -25 : 25]; 
        
        preStimIdx  = 1:floor(preStimWindow/1000*imagingFreq);
        postStimIdx = ceil(preStimWindow/1000*imagingFreq):(ceil(preStimWindow/1000*imagingFreq+ipi*nPulses/1000*imagingFreq));
    
        subSampMeans = [];
        subSampCI_bounds = cell(length(subSampIndex), 1);
        confidence_level = 0.95;
    
        for ff = 1:length(subSampIndex)
            randTrials = randperm(validTrials, subSampIndex(ff));
            subSampMeans(:, ff) = mean(validData(:, randTrials), 2);
            
            preStimVals  = mean(validData(preStimIdx,  randTrials), 1);   % 1 × 589
            postStimVals = mean(validData(postStimIdx, randTrials), 1);   % 1 × 589
            [h,p,ci,stats] = ttest(preStimVals, postStimVals);
            subSampPValues(1, ff) = p;
    
            if p<=0.05
                subSampConnMatrix(hh, ff) = 1;
            else if p>0.05
                subSampConnMatrix(hh, ff) = 0;
            end
            end
    
            subSamp_std_errors = std(validData(:, randTrials), 0, 2) / sqrt(length(randTrials));
            ci95 = 1.96*subSamp_std_errors;
            t_score = tinv((1 + confidence_level) / 2, length(randTrials) - 1);
            margin_of_error = t_score * subSamp_std_errors;
            lower_bounds = subSampMeans(:, ff) - margin_of_error;
            upper_bounds = subSampMeans(:, ff) + margin_of_error;
        
            subSampCI_bounds{ff}(:, 1) = lower_bounds;
            subSampCI_bounds{ff}(:, 2) = upper_bounds;
        end
        
        figure(10000+hh);
        set(gcf, 'Position',  [100, 100, 1800, 350])
        clf
        for ff = 1:length(subSampIndex)
            subplot(1, length(subSampIndex), length(subSampIndex)-ff+1)

            fill([linspace(0, size(subSampCI_bounds{ff}, 1)/imagingFreq, size(subSampCI_bounds{ff}, 1)), fliplr(linspace(0, size(subSampCI_bounds{ff}, 1)/imagingFreq, size(subSampCI_bounds{ff}, 1)))],...
            [subSampCI_bounds{ff}(:, 1)'*100, fliplr(subSampCI_bounds{ff}(:, 2)'*100)], [0.95, 0.95, 0.95], 'EdgeColor', [0.95, 0.95, 0.95]);
            hold on

            for nn = 1:length(nPulseCoords)
                xline(nPulseCoords(nn)/Fs, '-', 'LineWidth', 2, 'color', [1 0 0 0.1]);
            end

            plot(linspace(0, size(subSampMeans(:, ff), 1)/imagingFreq, size(subSampMeans(:, ff), 1)), subSampMeans(:, ff)*100, '-', 'linewidth', 2, 'color', 'g');
            ylabel('dF/F (%)');
            xlabel('Time (s)');
            hold on
            ylim([-0.5 3.5])
            ax.YColor = [0 1 0];
            if ff == 1
                title(['All ' num2str(subSampIndex(ff)), ' Trials;', ' p = ', num2str(subSampPValues(ff))]);     
            else title([num2str(subSampIndex(ff)), ' Trials subsample;', ' p = ', num2str(subSampPValues(ff))]);
            end
            hold off
        end
        % pause;
    end

    % Significance of average response as a function of subsampling
    respSigMatrix = fliplr(subSampConnMatrix);
    respSubSampIndex = fliplr(subSampIndex);
    % figure(5000+pp);
    % set(gcf, 'Position',  [100+(pp-1)*350, 100, 350, 1800])
    figure(5000);
    subplot(1, subSampPermutations, pp)

    % imagesc(respSigMatrix(:, 1:(size(respSigMatrix, 2)-1)))
    imagesc(respSigMatrix);
    caxis([0, 1]);
    axis image
    colormap(gray(2));
    set(gca,'XTick',[]);
    % xticklabels(respSubSampIndex(1:(size(respSubSampIndex, 2)-1)));
    if pp == 1
        % xlabel(['Random subsamples'], 'FontSize', 20);
        ylabel(['Candidate Presynaptic Target'], 'FontSize', 30);
    end
end

%% Power Analysis: Subsample Detection Rate
powerHatAllHolos = [];
pQuantAllHolos = cell(nHolos(2), 1);
targetSubsampAllHolos = zeros(nHolos(2), 1);

for hh = 1:nHolos(2)
    testData = [];
    testData = voltMapping.Cell1.exclFiltHoloSortedImagingAllTrials_commonF0{2, 1}{hh, 1};
    validData = [];
    invalidTrials = length(find(isnan(testData(1, :))));
    validTrials = length(find(~isnan(testData(1, :))));

    % Remove excluded (nan'd) trials
    for tt = 1:size(testData, 2)
        if ~isnan(testData(1, tt))
            validData = [validData, testData(:, tt)];
        end
    end

    % Define epochs (edit if needed, but should carry over from previous section)
    % preStimIdx  = 1:floor(preStimWindow/1000*imagingFreq);
    % postStimIdx = ceil(preStimWindow/1000*imagingFreq):(ceil(preStimWindow/1000*imagingFreq+ipi(2)*nPulses(2)/1000*imagingFreq));
    
    % Per-trial reduction (paired observations)
    preStimVals  = mean(validData(preStimIdx,  :), 1);   % 1 x nTrials
    postStimVals = mean(validData(postStimIdx, :), 1);   % 1 x nTrials
    
    nTrials = numel(preStimVals);
    if numel(postStimVals) ~= nTrials
        error('preVals and postVals sizes do not match.');
    end
    
    % Empirical power settings
    alpha = 0.05;
    
    % Trial counts to evaluate (edit as desired)
    % subSampIndexTop = floor(min(nValidTrials)/100)*100;
    % subSampIndexTop = validTrials-100+50;
    subSampIndexTop = floor(min(nValidTrials)/75)*75;
    subSampIndex = [validTrials, subSampIndexTop : -25 : 25];
    % subSampIndex = fliplr([validTrials, (floor(min(nValidTrials)/100)*100) : -25 : 25]);

    subSampIndex = fliplr(subSampIndex);
    
    nIter = 1000;        % number of Iterations per N
    rng(1);             % reproducible results
    
    powerHat = nan(size(subSampIndex));
    pQuant   = nan(numel(subSampIndex), 3);  % optional: p-value quantiles [0.1 0.5 0.9]
    
    % Repeated iterations of t-test of subsample blocks to test detection rate
    for ss = 1:numel(subSampIndex)
        subSamp = subSampIndex(ss);
    
        sig = false(nIter,1);
        pvals = nan(nIter,1);
    
        for ii = 1:nIter
            idx = randperm(nTrials, subSamp);
    
            % Paired t-test across trials on epoch-averaged values
            [h,p,ci,stats] = ttest(preStimVals(idx), postStimVals(idx), 'Alpha', alpha);
            pvals(ii) = p;
            sig(ii) = (p < alpha);
        end
    
        powerHat(ss) = mean(sig); 
        pQuant(ss,:) = quantile(pvals, [0.1 0.5 0.9]);
    
        fprintf('N=%d: power=%.3f\n', subSamp, powerHat(ss));
    end
    
    % Power curve
    figure(7000+hh); clf; hold on;
    plot(subSampIndex, powerHat, 'LineWidth', 2, 'Color', 'k');
    yline(0.80, '--', '0.80', 'LabelHorizontalAlignment','left');
    yline(0.90, '--', '0.90', 'LabelHorizontalAlignment','left');
    yline(0.95, '--', '0.95', 'LabelHorizontalAlignment','left', 'Color', 'r');
    xlabel('Number of trials (subsample size)');
    ylabel(sprintf('Empirical power (P(p < %.2f))', alpha));
    title(sprintf('Empirical power curve (paired t-test, %d iters per)', nIter));
    ylim([0 1]);
    set(gca, 'FontSize', 16);
    box on;
    
    % P-value distribution (10th, median, 90th percentile) as a function of subsample blocks
    figure(8000+hh); clf; hold on;
    plot(subSampIndex, pQuant(:,2), 'LineWidth', 2, 'Color', 'k');           % median p
    plot(subSampIndex, pQuant(:,1), '--', 'LineWidth', 1, 'Color', [0.75, 0.75, 0.75]);     % 10th percentile p
    plot(subSampIndex, pQuant(:,3), '--', 'LineWidth', 1, 'Color', [0.75, 0.75, 0.75]);     % 90th percentile p
    set(gca, 'YScale', 'log');
    yline(alpha, '--', sprintf('\\alpha=%.2f', alpha), 'Color', 'r');
    xlabel('Number of trials (subsample size)');
    ylabel('p-value (log scale)');
    title('p-value distribution across subsamples (M, P10, P90)');
    set(gca, 'FontSize', 16);
    box on;
    
    % Extrapolation of trials meeting 0.95 power
    thresholdPow = 0.95;
    % Enforce monotone non-decreasing power (removes Monte Carlo wiggles)
    powMono = cummax(powerHat);
    % Find first index at/above target
    indexPostThresh = find(powMono >= thresholdPow, 1, 'first');
    
    if isempty(indexPostThresh)
        warning('Target power not reached within total trials. Insufficient trials or response is nonsignificant');
        targetSubsamp = NaN;
    else
        if indexPostThresh == 1
            targetSubsamp = subSampIndex(1);
        else
            % Linear interpolation between the two surrounding points
            subsampPreThresh = subSampIndex(indexPostThresh-1);  
            powerPreThresh = powMono(indexPostThresh-1);
            subsampPostThresh = subSampIndex(indexPostThresh);    
            powerPostThresh = powMono(indexPostThresh);
    
            % Avoid division by zero if P1==P2
            if powerPostThresh == powerPreThresh
                targetSubsamp = subsampPostThresh;
            else
                targetSubsamp = subsampPreThresh + (thresholdPow - powerPreThresh) * (subsampPostThresh - subsampPreThresh) / (powerPostThresh - powerPreThresh);
            end
        end
    end
    
    targetSubsampAllHolos(hh, 1) = targetSubsamp;
    fprintf('Estimated N for power %.2f: %.1f trials\n', thresholdPow, targetSubsamp);

    powerHatAllHolos(hh, :) = powerHat;
    pQuantAllHolos{hh, 1} = pQuant;

    % Compute and report effect size on full dataset
    diffVals = postStimVals - preStimVals;
    cohens_d_paired = mean(diffVals) / std(diffVals);
    
    [~, pFull, ciFull, statsFull] = ttest(preStimVals, postStimVals, 'Alpha', alpha);
    
    fprintf('\nFull dataset results:\n');
    fprintf('Paired t-test: t(%d)=%.3f, p=%.3g\n', statsFull.df, statsFull.tstat, pFull);
    fprintf('Mean diff (post-pre)=%.6g, 95%% CI=[%.6g, %.6g]\n', mean(diffVals), ciFull(1), ciFull(2));
    fprintf('Cohen''s d (paired)=%.3f\n', cohens_d_paired);
end

meanTargetTrials = nanmean(targetSubsampAllHolos);

% Matrix of power (average of all iterations) across the subsample blocks
for hh = 1:nHolos(2)
    binaryPowerHat(hh, find(powerHatAllHolos(hh, :) <= 0.95)) = 0;
    binaryPowerHat(hh, find(powerHatAllHolos(hh, :) >= 0.95)) = 1;
end

figure(9000); clf;
set(gcf, 'Position',  [350, 100, 250, 1800])
imagesc(binaryPowerHat(:, 1:size(binaryPowerHat, 2)-1));
caxis([0, 1]);
axis image
colormap(gray(2));
set(gca,'XTick',[]);
set(gca,'XGrid', ' on', 'YGrid', 'on', 'GridColor', [0.1 0.1 0.1]);
xlabel(['Random subsamples'], 'FontSize', 20);
ylabel(['Candidate Presynaptic Target'], 'FontSize', 30);
% xticks(1.5:1:size(binaryPowerHat, 2));
% yticks(1.5:1:size(binaryPowerHat, 1));
xticklabels([]);
