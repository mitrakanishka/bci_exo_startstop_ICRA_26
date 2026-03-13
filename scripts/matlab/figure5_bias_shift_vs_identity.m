function figure5_bias_shift_vs_identity()
% Paper Figure 5: task-based recentering bias shift vs identity.

clc;
script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(script_dir));
addpath(fullfile(script_dir, '_shared'));

of_dir = fullfile(repo_root, 'BCI_Harmony_ExperimentalData', 'epoched_data_of');
on_dir = fullfile(repo_root, 'BCI_Harmony_ExperimentalData', 'epoched_data_on');
outdir = fullfile(repo_root, 'figures');
subject_csv = fullfile(repo_root, 'fig_data', 'fig5_bias_task_vs_identity.csv');
summary_csv = fullfile(repo_root, 'fig_data', 'fig5_summary_stats.csv');

[bias_tbl, ~] = compute_recentering_tables(of_dir, on_dir, 8, 0.30, 1e-3, 0.90);
ensure_parent(subject_csv);
writetable(bias_tbl, subject_csv);
summary_tbl = build_summary(bias_tbl);
writetable(summary_tbl, summary_csv);

fig = figure('Color', 'w', 'Position', [100 100 1100 540], 'Renderer', 'painters');
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Task Recentering vs Identity (\Delta Median Margin)', 'FontSize', 16, 'FontWeight', 'bold');

phase_order = {"ONSET", "OFFSET"};
metric_order = {"Delta_pos", "Delta_neg", "Delta_sep"};
colors.ONSET = [0.1804 0.4902 0.1961; 0.1176 0.5333 0.8980; 0.7882 0.7882 0.7882];
colors.OFFSET = [0.8039 0.3608 0.3608; 0.4941 0.3412 0.7608; 0.7882 0.7882 0.7882];
labels.ONSET = {'\Delta_{pos}\newline(Start-MI)','\Delta_{neg}\newline(REST)','\Delta_{sep}'};
labels.OFFSET = {'\Delta_{pos}\newline(Stop-MI)','\Delta_{neg}\newline(Maintain-MI)','\Delta_{sep}'};

for pi = 1:numel(phase_order)
    phase = phase_order{pi};
    ax = nexttile(tl, pi); hold(ax, 'on');
    ax.Toolbar.Visible = 'off';
    phase_rows = summary_tbl(summary_tbl.phase == phase, :);
    phase_rows = sort_summary_rows(phase_rows, metric_order);
    x = 1:numel(metric_order);
    means = phase_rows.mean';
    lo = phase_rows.ci_low';
    hi = phase_rows.ci_high';
    err = [means - lo; hi - means];
    cols = colors.(phase);

    bh = bar(ax, x, means, 0.8, 'FaceColor', 'flat', 'EdgeColor', 'black', 'LineWidth', 0.5);
    bh.CData = cols;
    errorbar(ax, x, means, err(1, :), err(2, :), 'k.', 'LineWidth', 1.0, 'CapSize', 10);
    yline(ax, 0.0, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.0);

    phase_vals = bias_tbl.value(bias_tbl.phase == phase & bias_tbl.metric == "Delta_sep");
    phase_vals = phase_vals(isfinite(phase_vals));
    res = paired_wilcoxon_exact(zeros(size(phase_vals)), phase_vals);
    if isfinite(res.pvalue)
        ptxt = sprintf('\\Delta_{sep} vs S_I: p=%0.4f', res.pvalue);
    else
        ptxt = sprintf('\\Delta_{sep} vs S_I: n=%d', res.n);
    end

    text(ax, 0.5, 0.04, ptxt, 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 15, 'BackgroundColor', 'white', 'Margin', 6);

    ax.XTick = x;
    ax.XTickLabel = labels.(phase);
    ax.TickLabelInterpreter = 'tex';
    ax.FontSize = 14;
    title(ax, lower_title(phase), 'FontSize', 16, 'FontWeight', 'bold');
    if pi == 1
        ylabel(ax, '\Delta vs S_I (median margin)', 'Interpreter', 'tex', 'FontSize', 18, 'FontWeight', 'bold');
    end
    grid(ax, 'on');
    ax.YGrid = 'on';
    ax.XGrid = 'off';
    ax.GridLineStyle = ':';
    box(ax, 'off');
end

ensure_parent(fullfile(outdir, 'dummy'));
exportgraphics(fig, fullfile(outdir, 'fig5_bias_shift_vs_identity.png'), 'Resolution', 300);
exportgraphics(fig, fullfile(outdir, 'fig5_bias_shift_vs_identity.pdf'), 'ContentType', 'vector');
close(fig);
end

function summary_tbl = build_summary(bias_tbl)
phase_order = {"ONSET", "OFFSET"};
metric_order = {"Delta_pos", "Delta_neg", "Delta_sep"};
phase = strings(0,1);
metric = strings(0,1);
mean_val = [];
ci_low = [];
ci_high = [];
n = [];

for pi = 1:numel(phase_order)
    for mi = 1:numel(metric_order)
        vals = bias_tbl.value(bias_tbl.phase == phase_order{pi} & bias_tbl.metric == metric_order{mi});
        vals = vals(isfinite(vals));
        phase(end + 1, 1) = string(phase_order{pi}); %#ok<AGROW>
        metric(end + 1, 1) = string(metric_order{mi}); %#ok<AGROW>
        if isempty(vals)
            mean_val(end + 1, 1) = NaN; %#ok<AGROW>
            ci_low(end + 1, 1) = NaN; %#ok<AGROW>
            ci_high(end + 1, 1) = NaN; %#ok<AGROW>
            n(end + 1, 1) = 0; %#ok<AGROW>
        else
            seed = (phase_order{pi} == "ONSET") * 7 + (phase_order{pi} == "OFFSET") * 17 + (mi - 1);
            [lo, hi] = bootstrap_ci(vals, 10000, seed);
            mean_val(end + 1, 1) = mean(vals); %#ok<AGROW>
            ci_low(end + 1, 1) = lo; %#ok<AGROW>
            ci_high(end + 1, 1) = hi; %#ok<AGROW>
            n(end + 1, 1) = numel(vals); %#ok<AGROW>
        end
    end
end
summary_tbl = table(phase, metric, mean_val, ci_low, ci_high, n, ...
    'VariableNames', {'phase', 'metric', 'mean', 'ci_low', 'ci_high', 'n'});
end

function T = sort_summary_rows(T, metric_order)
idx = zeros(height(T), 1);
for i = 1:numel(metric_order)
    idx(T.metric == metric_order{i}) = i;
end
[~, order] = sort(idx);
T = T(order, :);
end

function [lo, hi] = bootstrap_ci(vals, n_boot, seed)
vals = vals(isfinite(vals));
if isempty(vals)
    lo = NaN;
    hi = NaN;
    return;
end
rng(seed);
sample_idx = randi(numel(vals), numel(vals), n_boot);
means = mean(vals(sample_idx), 1);
lo = percentile_linear(means, 2.5);
hi = percentile_linear(means, 97.5);
end

function value = percentile_linear(x, pct)
x = sort(x(:));
if isempty(x)
    value = NaN;
    return;
end
idx = 1 + (numel(x) - 1) * (pct / 100);
lo = floor(idx);
hi = ceil(idx);
if lo == hi
    value = x(lo);
else
    value = x(lo) + (idx - lo) * (x(hi) - x(lo));
end
end

function out = lower_title(phase)
if phase == "ONSET"
    out = 'Onset';
else
    out = 'Offset';
end
end

function ensure_parent(pathstr)
folder = fileparts(pathstr);
if exist(folder, 'dir') ~= 7
    mkdir(folder);
end
end
