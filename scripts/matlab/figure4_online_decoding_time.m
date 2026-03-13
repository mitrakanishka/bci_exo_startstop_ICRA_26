function figure4_online_decoding_time()
% Paper Figure 4: online decoding-time distributions (onset/offset).

clc;
script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(script_dir));
addpath(fullfile(script_dir, '_shared'));

log_root = fullfile(repo_root, 'BCI_Harmony_ExperimentalData', 'online_python_log');
outdir = fullfile(repo_root, 'figures');
out_csv = fullfile(repo_root, 'fig_data', 'fig4_session2_vs3_stats.csv');

T = load_all_online_logs(log_root);
T = T(ismember(T.session, [2 3]), :);
RT = build_rt_table(T);
Stats = paired_session_stats(RT);
ensure_parent(out_csv);
writetable(Stats, out_csv);

fig = figure('Color', 'w', 'Position', [100 100 820 430], 'Renderer', 'painters');
ax = axes(fig); hold(ax, 'on');
ax.Toolbar.Visible = 'off';
centers = containers.Map({'Onset', 'Offset'}, [1.0, 2.0]);
offsets = containers.Map({'2', '3'}, [-0.13, 0.13]);
colors = containers.Map({'2', '3'}, {[0.2980 0.4470 0.6900], [0.3333 0.6588 0.4078]});
width = 0.34;
phases = {'Onset', 'Offset'};
sessions = [2 3];

for pi = 1:numel(phases)
    phase = phases{pi};
    for sj = 1:numel(sessions)
        sess = sessions(sj);
        vals_real = RT.rt_s(strcmp(RT.phase, phase) & RT.session == sess);
        if isempty(vals_real)
            continue;
        end
        pos = centers(phase) + offsets(num2str(sess));
        vals_violin = vals_real;
        if strcmp(phase, 'Offset')
            vals_violin = [vals_real; 0];
        end
        draw_violin(ax, vals_violin, pos, width, colors(num2str(sess)));
        rng(42 + pi + sj);
        jitter = (rand(size(vals_real)) - 0.5) * 0.10;
        scatter(ax, pos + jitter, vals_real, 11, colors(num2str(sess)), 'filled', ...
            'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.0);
        m = mean(vals_real);
        span = violin_span(vals_violin, pos, width, m);
        plot(ax, span, [m m], 'k-', 'LineWidth', 1.0);
    end
end

yline(ax, 0, 'Color', 'none');
ax.XLim = [0.55 2.45];
ax.YLim = [0 6];
ax.YTick = 0:1:6;
ax.XTick = [1 2];
ax.XTickLabel = phases;
set(ax, 'FontSize', 11);
for i = 1:numel(ax.XTickLabel)
    ax.XAxis.TickLabel{i} = phases{i};
end
xlabel(ax, '');
ylabel(ax, 'Decoding Time (s)', 'FontSize', 13, 'FontWeight', 'bold');
title(ax, 'Online Decoding Time', 'FontSize', 16, 'FontWeight', 'bold');
grid(ax, 'on');
ax.GridLineStyle = '-';
ax.GridAlpha = 0.32;
box(ax, 'off');

leg1 = patch(NaN, NaN, colors('2'), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
leg2 = patch(NaN, NaN, colors('3'), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
legend(ax, [leg1 leg2], {'Session 2', 'Session 3'}, 'Location', 'southeast', 'Box', 'on');

ensure_parent(fullfile(outdir, 'dummy'));
exportgraphics(fig, fullfile(outdir, 'fig4_online_decoding_time.png'), 'Resolution', 300);
exportgraphics(fig, fullfile(outdir, 'fig4_online_decoding_time.pdf'), 'ContentType', 'vector');
close(fig);
end

function RT = build_rt_table(T)
Ton = T(T.onset_label == "hit", {'subject', 'session', 'onset_rt_s'});
Ton = renamevars(Ton, 'onset_rt_s', 'rt_s');
Ton.phase = repmat("Onset", height(Ton), 1);

Toff = T(T.offset_attempted & T.offset_label == "hit", {'subject', 'session', 'offset_rt_s'});
Toff = renamevars(Toff, 'offset_rt_s', 'rt_s');
Toff.rt_s = Toff.rt_s + 2.0;
Toff.phase = repmat("Offset", height(Toff), 1);

RT = [Ton; Toff];
RT = RT(isfinite(RT.rt_s) & RT.rt_s <= 6.0, :);
end

function Stats = paired_session_stats(RT)
phases = {"Onset", "Offset"};
rows = table();
for i = 1:numel(phases)
    phase = phases{i};
    Tp = RT(RT.phase == phase, :);
    subs = unique(Tp.subject);
    pairs = nan(numel(subs), 2);
    keep = false(numel(subs), 1);
    for si = 1:numel(subs)
        for sess = 2:3
            vals = Tp.rt_s(Tp.subject == subs(si) & Tp.session == sess);
            if ~isempty(vals)
                pairs(si, sess - 1) = median(vals);
            end
        end
        keep(si) = all(isfinite(pairs(si, :)));
    end
    pairs = pairs(keep, :);
    if size(pairs, 1) < 2
        rows = [rows; table(string(phase), 0, NaN, NaN, NaN, NaN, ...
            'VariableNames', {'phase', 'n', 'w_stat', 'p_value', 'mean_delta_s', 'sd_delta_s'})]; %#ok<AGROW>
        continue;
    end
    delta = pairs(:, 2) - pairs(:, 1);
    res = paired_wilcoxon_exact(pairs(:, 1), pairs(:, 2));
    rows = [rows; table(string(phase), res.n, res.statistic, res.pvalue, mean(delta), std(delta, 0), ...
        'VariableNames', {'phase', 'n', 'w_stat', 'p_value', 'mean_delta_s', 'sd_delta_s'})]; %#ok<AGROW>
end
Stats = rows;
end

function draw_violin(ax, vals, pos, width, color)
[ygrid, dens] = simple_density(vals, linspace(0, 6, 300));
if max(dens) <= 0
    return;
end
half = (dens / max(dens)) * (width / 2);
patch(ax, [pos - half, fliplr(pos + half)], [ygrid, fliplr(ygrid)], color, ...
    'FaceAlpha', 0.35, 'EdgeColor', 'none');
end

function span = violin_span(vals, pos, width, y)
[ygrid, dens] = simple_density(vals, linspace(0, 6, 300));
if max(dens) <= 0
    span = [pos - width * 0.24, pos + width * 0.24];
    return;
end
half = (dens / max(dens)) * (width / 2);
h = interp1(ygrid, half, y, 'linear', width * 0.24);
span = [pos - h, pos + h];
end

function [ygrid, dens] = simple_density(vals, ygrid)
vals = vals(:);
vals = vals(isfinite(vals));
if isempty(vals)
    dens = zeros(size(ygrid));
    return;
end
sd = std(vals);
if sd <= 0 || ~isfinite(sd)
    sd = 0.25;
end
bw = 1.06 * sd * numel(vals)^(-1/5);
if ~isfinite(bw) || bw <= 0
    bw = 0.25;
end
z = (ygrid(:) - vals.') / bw;
dens = mean(exp(-0.5 * z.^2) / (sqrt(2 * pi) * bw), 2).';
end

function ensure_parent(pathstr)
folder = fileparts(pathstr);
if exist(folder, 'dir') ~= 7
    mkdir(folder);
end
end
