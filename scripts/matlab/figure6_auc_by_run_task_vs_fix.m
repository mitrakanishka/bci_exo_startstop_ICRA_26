function figure6_auc_by_run_task_vs_fix()
% Paper Figure 6: AUC by run for TASK vs fixation-based recentering.

clc;
script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(script_dir));
addpath(fullfile(script_dir, '_shared'));

of_dir = fullfile(repo_root, 'BCI_Harmony_ExperimentalData', 'epoched_data_of');
on_dir = fullfile(repo_root, 'BCI_Harmony_ExperimentalData', 'epoched_data_on');
outdir = fullfile(repo_root, 'figures');
subject_csv = fullfile(repo_root, 'fig_data', 'fig6_auc_subject_run.csv');
summary_csv = fullfile(repo_root, 'fig_data', 'fig6_run_level_summary.csv');

[~, auc_tbl] = compute_recentering_tables(of_dir, on_dir, 8, 0.30, 1e-3, 0.90);
ensure_parent(subject_csv);
writetable(auc_tbl, subject_csv);
[summary_tbl, pvals] = build_summary(auc_tbl);
writetable(summary_tbl, summary_csv);

fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 7.2 8.0], 'Renderer', 'painters');
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'AUC: Task vs Fixation-Based Recentering', 'FontSize', 20, 'FontWeight', 'bold');

phase_order = {"ONSET", "OFFSET"};
scheme_order = {"TASK", "FIX"};
colors.TASK = [0.9020 0.5255 0.1647];
colors.FIX = [0.1608 0.4745 1.0000];
runs = 1:8;
labels = arrayfun(@run_label, runs, 'UniformOutput', false);
axes_list = gobjects(numel(phase_order), 1);
legend_handles = gobjects(3, 1);

for pi = 1:numel(phase_order)
    phase = phase_order{pi};
    ax = nexttile(tl, pi); hold(ax, 'on');
    ax.Toolbar.Visible = 'off';
    axes_list(pi) = ax;
    phase_rows = summary_tbl(summary_tbl.phase == phase, :);
    for si = 1:numel(scheme_order)
        scheme = scheme_order{si};
        rows = phase_rows(phase_rows.scheme == scheme, :);
        rows = sortrows(rows, 'run');
        x = rows.run';
        y = rows.mean';
        lo = rows.ci_low';
        hi = rows.ci_high';
        fill(ax, [x, fliplr(x)], [lo, fliplr(hi)], colors.(scheme), ...
            'FaceAlpha', 0.14, 'EdgeColor', 'none');
        h = plot(ax, x, y, '-o', 'Color', colors.(scheme), 'LineWidth', 2.0, ...
            'MarkerSize', 3.2, 'DisplayName', scheme);
        if pi == 1
            legend_handles(si) = h;
        end
    end
    hchance = yline(ax, 0.5, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.0, 'DisplayName', 'Chance');
    if pi == 1
        legend_handles(3) = hchance;
    end
    title(ax, lower_title(phase), 'FontSize', 13, 'FontWeight', 'bold');
    ylabel(ax, 'AUC', 'FontSize', 12, 'FontWeight', 'bold');
    xlim(ax, [1 8]);
    ylim(ax, [0.30 1.01]);
    ax.XTick = runs;
    ax.XTickLabel = labels;
    ax.FontSize = 11;
    text(ax, 0.98, 0.04, p_text(pvals.(phase)), 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
        'FontSize', 11, 'BackgroundColor', 'white', 'Margin', 4);
    grid(ax, 'on');
    ax.GridAlpha = 0.2;
    box(ax, 'off');
    if pi == 1
        legend(ax, legend_handles, {'TASK', 'FIX', 'Chance'}, 'Location', 'southwest', 'Box', 'on', 'FontSize', 11);
    end
end
xlabel(axes_list(end), 'Session Number, Run Number', 'FontSize', 12, 'FontWeight', 'bold');

ensure_parent(fullfile(outdir, 'dummy'));
exportgraphics(fig, fullfile(outdir, 'fig6_auc_by_run_task_vs_fix.png'), 'Resolution', 300);
exportgraphics(fig, fullfile(outdir, 'fig6_auc_by_run_task_vs_fix.pdf'), 'ContentType', 'vector');
close(fig);
end

function [summary_tbl, pvals] = build_summary(auc_tbl)
phase_order = {"ONSET", "OFFSET"};
scheme_order = {"TASK", "FIX"};
runs = 1:8;
phase = strings(0,1);
scheme = strings(0,1);
run = [];
mean_val = [];
sd_val = [];
ci_low = [];
ci_high = [];
n = [];
pvals = struct('ONSET', NaN, 'OFFSET', NaN);

for pi = 1:numel(phase_order)
    curr_phase = phase_order{pi};
    for si = 1:numel(scheme_order)
        curr_scheme = scheme_order{si};
        for ri = 1:numel(runs)
            vals = auc_tbl.auc(auc_tbl.phase == curr_phase & auc_tbl.scheme == curr_scheme & auc_tbl.run == runs(ri));
            vals = vals(isfinite(vals));
            phase(end + 1, 1) = string(curr_phase); %#ok<AGROW>
            scheme(end + 1, 1) = string(curr_scheme); %#ok<AGROW>
            run(end + 1, 1) = runs(ri); %#ok<AGROW>
            if isempty(vals)
                mean_val(end + 1, 1) = NaN; %#ok<AGROW>
                sd_val(end + 1, 1) = NaN; %#ok<AGROW>
                ci_low(end + 1, 1) = NaN; %#ok<AGROW>
                ci_high(end + 1, 1) = NaN; %#ok<AGROW>
                n(end + 1, 1) = 0; %#ok<AGROW>
            else
                m = mean(vals);
                sd = std(vals, 1);
                half = 1.96 * sd / sqrt(numel(vals));
                mean_val(end + 1, 1) = m; %#ok<AGROW>
                sd_val(end + 1, 1) = sd; %#ok<AGROW>
                ci_low(end + 1, 1) = m - half; %#ok<AGROW>
                ci_high(end + 1, 1) = m + half; %#ok<AGROW>
                n(end + 1, 1) = numel(vals); %#ok<AGROW>
            end
        end
    end

    subs = unique(auc_tbl.subj(auc_tbl.phase == curr_phase));
    task_vals = [];
    fix_vals = [];
    for si = 1:numel(subs)
        task = auc_tbl.auc(auc_tbl.phase == curr_phase & auc_tbl.subj == subs(si) & auc_tbl.scheme == "TASK");
        fix = auc_tbl.auc(auc_tbl.phase == curr_phase & auc_tbl.subj == subs(si) & auc_tbl.scheme == "FIX");
        task = task(isfinite(task));
        fix = fix(isfinite(fix));
        if isempty(task) || isempty(fix)
            continue;
        end
        task_vals(end + 1, 1) = mean(task); %#ok<AGROW>
        fix_vals(end + 1, 1) = mean(fix); %#ok<AGROW>
    end
    res = paired_wilcoxon_exact(task_vals, fix_vals);
    pvals.(curr_phase) = res.pvalue;
end

summary_tbl = table(phase, scheme, run, mean_val, sd_val, ci_low, ci_high, n, ...
    'VariableNames', {'phase', 'scheme', 'run', 'mean', 'sd', 'ci_low', 'ci_high', 'n'});
end

function label = run_label(run_idx)
if run_idx <= 4
    label = sprintf('S2,R%d', run_idx);
else
    label = sprintf('S3,R%d', run_idx - 4);
end
end

function out = lower_title(phase)
if phase == "ONSET"
    out = 'Onset';
else
    out = 'Offset';
end
end

function txt = p_text(p)
if isfinite(p)
    txt = sprintf('p=%0.4f', p);
else
    txt = 'p=NA';
end
end

function ensure_parent(pathstr)
folder = fileparts(pathstr);
if exist(folder, 'dir') ~= 7
    mkdir(folder);
end
end
