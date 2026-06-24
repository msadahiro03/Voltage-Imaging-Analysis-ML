% Single-column scatter: y = [11 30 9]
fig = figure;
set(fig, 'Units', 'pixels', 'Position', [100 100 100 250], 'Color', 'w');

y = [11; 30; 9];
x = ones(numel(y), 1);
m = mean(y);
sem_val = std(y) / sqrt(numel(y)); % sample SEM

scatter(x, y, 48, 'k', 'filled', 'MarkerEdgeColor', 'k');
hold on;
errorbar(1, m, sem_val, 'k', 'LineStyle', 'none', 'LineWidth', 1.2, ...
    'CapSize', 10);

xlim([0.5 1.5]);
set(gca, 'XTick', []);
ypad = 0.1 * (max([y; m + sem_val]) + eps);
ylim([0, max([y; m + sem_val]) + ypad]);

ax = gca;
set(ax, 'Color', 'w', 'Box', 'off', 'LineWidth', 1.5, ...
    'XColor', 'k', 'YColor', 'k', ...
    'GridColor', [0.15 0.15 0.15], 'MinorGridColor', [0.5 0.5 0.5]);
set(ax, 'FontSize', 14); % tick labels
ylabel('Connectivity (%)', 'FontSize', 14, 'Color', 'k');

dy = diff(ax.YLim);
text(1, m + sem_val + 0.04 * dy, ...
    sprintf('mean = %.1f%%\n(n = %d)', m, numel(y)), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'FontSize', 14, 'Color', 'k', 'Interpreter', 'none');
