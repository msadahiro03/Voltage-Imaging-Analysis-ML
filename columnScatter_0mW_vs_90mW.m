% Column scatter: paired rows (connector lines); combined means as black bars; circles for data.

% MATLAB variable names cannot start with a digit; use mW0 / mW90 aliases.

mW0 = [0.01 0.0 0.07 0.06 0.05 0.04 0.04 0.02 0.02];
mW90 = [0.125 0.10 0.25 0.25 0.5 0.6 0.49 0.32 0.08];
mW0ThisCell = [0.0500; 0.0500; 0.0750; 0.0500];
mW90ThisCell = [0.4500; 0.2000; 0.5750; 0.3500];

x0 = 1;
x90 = 2;

nPop = numel(mW0);
nCell = numel(mW0ThisCell);
assert(numel(mW90) == nPop, 'mW0 and mW90 must be the same length (paired rows).');
assert(numel(mW90ThisCell) == nCell, ...
    'mW0ThisCell and mW90ThisCell must be the same length (paired rows).');

% Light horizontal jitter so overlapping points stay visible in one column
rng(1);
j0pop = x0 + 0.06 * (rand(1, nPop) - 0.5);
j90pop = x90 + 0.06 * (rand(1, nPop) - 0.5);
j0cell = x0 + 0.06 * (rand(1, nCell) - 0.5);
j90cell = x90 + 0.06 * (rand(1, nCell) - 0.5);

mean0 = mean([mW0(:); mW0ThisCell(:)]);
mean90 = mean([mW90(:); mW90ThisCell(:)]);

% Publication-style figure (light mode, single-column width)
figW = 3.6;  % inches
figH = 3.1;
figure('Color', 'w', 'Units', 'inches', ...
    'Position', [1 1 figW figH], 'PaperPositionMode', 'auto');

hold on;

colPop = [0.25 0.25 0.25];
colCell = [0.15 0.55 0.22];

meanBarHalfW = 0.14;  % half-width of mean marker in x data units

% Paired observations: line k connects row k of mW0 to row k of mW90 (same for this cell)
for k = 1:nPop
    plot([j0pop(k), j90pop(k)], [mW0(k), mW90(k)], '-', ...
        'Color', [0.78 0.78 0.78], 'LineWidth', 0.85, 'Clipping', 'on', ...
        'HandleVisibility', 'off');
end
for k = 1:nCell
    plot([j0cell(k), j90cell(k)], [mW0ThisCell(k), mW90ThisCell(k)], '-', ...
        'Color', [0.72 0.86 0.74], 'LineWidth', 0.9, 'Clipping', 'on', ...
        'HandleVisibility', 'off');
end

scatter(j0pop(:), mW0(:), 42, colPop, 'filled', 'Marker', 'o', ...
    'MarkerFaceAlpha', 0.75, 'LineWidth', 0.4, 'MarkerEdgeColor', [0.15 0.15 0.15], ...
    'DisplayName', 'Population (paired)');
scatter(j90pop(:), mW90(:), 42, colPop, 'filled', 'Marker', 'o', ...
    'MarkerFaceAlpha', 0.75, 'LineWidth', 0.4, 'MarkerEdgeColor', [0.15 0.15 0.15], ...
    'HandleVisibility', 'off');

scatter(j0cell(:), mW0ThisCell(:), 42, colCell, 'filled', 'Marker', 'o', ...
    'LineWidth', 0.5, 'MarkerEdgeColor', [0 0.38 0.12], ...
    'DisplayName', 'This cell (paired)');
scatter(j90cell(:), mW90ThisCell(:), 42, colCell, 'filled', 'Marker', 'o', ...
    'LineWidth', 0.5, 'MarkerEdgeColor', [0 0.38 0.12], ...
    'HandleVisibility', 'off');

plot([x0 - meanBarHalfW, x0 + meanBarHalfW], [mean0, mean0], 'k-', ...
    'LineWidth', 2.8, 'Clipping', 'on', 'DisplayName', 'Mean (combined)');
plot([x90 - meanBarHalfW, x90 + meanBarHalfW], [mean90, mean90], 'k-', ...
    'LineWidth', 2.8, 'Clipping', 'on', 'HandleVisibility', 'off');

hold off;

xlim([0.45 2.55]);
ylim tight;
yL = ylim;
ylim([yL(1) - 0.03 * diff(yL), yL(2) + 0.03 * diff(yL)]);

ax = gca;
set(ax, 'XTick', [x0 x90], 'XTickLabel', {'0 mW', '90 mW'}, ...
    'TickDir', 'out', 'TickLength', [0.022 0.022], ...
    'LineWidth', 1.1, 'FontSize', 11, 'FontName', 'Arial', ...
    'XColor', [0 0 0], 'YColor', [0 0 0], 'Color', 'w', ...
    'Box', 'off', 'Layer', 'top');
grid off;

ylabel('Value', 'FontSize', 12, 'FontName', 'Arial', 'Color', [0 0 0]);
title('Population vs. this cell', 'FontSize', 12, 'FontName', 'Arial', ...
    'FontWeight', 'normal', 'Color', [0 0 0]);

lgd = legend('Location', 'southoutside', 'Orientation', 'horizontal', ...
    'Box', 'off', 'FontSize', 9, 'FontName', 'Arial', 'TextColor', [0 0 0]);
lgd.ItemTokenSize = [18, 10];

set(gcf, 'InvertHardcopy', 'off');
