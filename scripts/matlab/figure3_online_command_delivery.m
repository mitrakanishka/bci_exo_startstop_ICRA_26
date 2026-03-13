function figure3_online_command_delivery()
% Paper Figure 3: online command-delivery accuracy (onset/offset).

clc;
script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(script_dir));
addpath(fullfile(script_dir, '_shared'));

log_root = fullfile(repo_root, 'BCI_Harmony_ExperimentalData', 'online_python_log');
outdir = fullfile(repo_root, 'figures');
out_csv = fullfile(repo_root, 'fig_data', 'fig3_group_composition.csv');

T = load_all_online_logs(log_root);
sessions = [2 3];
outcomes = ["hit", "miss", "timeout"];
phases = ["Onset", "Offset"];
colors = [0.3333 0.6588 0.4078; 0.8824 0.4863 0.0196; 0.2980 0.4706 0.6588];

Ton = subject_session_props(T(ismember(T.session, sessions) & ismember(T.onset_label, outcomes), :), ...
    'onset_label', outcomes, sessions);
Ton.phase = repmat("Onset", height(Ton), 1);
Toff_in = T(ismember(T.session, sessions) & T.offset_attempted & ismember(T.offset_label, outcomes), :);
Toff = subject_session_props(Toff_in, 'offset_label', outcomes, sessions);
Toff.phase = repmat("Offset", height(Toff), 1);
summary = [Ton; Toff];
summary = fill_missing_rows(summary, sessions, phases, outcomes);
summary = summary(:, {'session', 'phase', 'outcome', 'mean_prop'});
ensure_parent(out_csv);
writetable(summary, out_csv);

x = [1 2 4 5];
plot_sessions = [2; 3; 2; 3];
plot_phases = ["Onset"; "Onset"; "Offset"; "Offset"];
Y = zeros(4, 3);
for i = 1:4
    for j = 1:numel(outcomes)
        mask = summary.session == plot_sessions(i) & summary.phase == plot_phases(i) & summary.outcome == outcomes(j);
        Y(i, j) = summary.mean_prop(find(mask, 1));
    end
end

fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 5.8 3.4], 'Renderer', 'painters');
ax = axes(fig); hold(ax, 'on');
ax.Toolbar.Visible = 'off';
ax.Clipping = 'off';
ax.Position = [0.12 0.22 0.78 0.68];

bh = bar(ax, x, Y, 0.9, 'stacked', 'LineWidth', 0.5, 'EdgeColor', 'black');
for j = 1:numel(outcomes)
    bh(j).FaceColor = colors(j, :);
end

for i = 1:numel(x)
    if Y(i, 1) > 0.06
        text(ax, x(i), Y(i, 1) / 2, sprintf('%0.0f%%', Y(i, 1) * 100), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 13, 'FontWeight', 'bold', 'Color', 'w');
    end
    bottoms = cumsum(Y(i, 1:2));
    for j = 1:numel(bottoms)
        if bottoms(j) > 0
            plot(ax, [x(i) - 0.45, x(i) + 0.45], [bottoms(j) bottoms(j)], 'w-', 'LineWidth', 1.0);
        end
    end
end

ax.YLim = [0 1];
ax.XLim = [0.3 5.7];
ax.XTick = [];
ax.YTick = 0:0.2:1;
ax.YTickLabel = compose('%d%%', round(ax.YTick * 100));
ax.FontSize = 11;
ylabel(ax, 'Percentage of Trials (%)', 'FontSize', 13, 'FontWeight', 'bold');
title(ax, 'Online Command-Delivery Accuracy', 'FontSize', 16, 'FontWeight', 'bold');

grid(ax, 'on');
ax.YGrid = 'on';
ax.XGrid = 'off';
ax.GridLineStyle = ':';
ax.GridAlpha = 0.55;
box(ax, 'on');

legend(ax, [bh(1) bh(2) bh(3)], {'Hit', 'Miss', 'Timeout'}, 'Location', 'southeast', 'FontSize', 7, 'Box', 'on');

add_bottom_label(fig, ax, 1, 'Session 2', 0.09, 11, 'normal');
add_bottom_label(fig, ax, 2, 'Session 3', 0.09, 11, 'normal');
add_bottom_label(fig, ax, 4, 'Session 2', 0.09, 11, 'normal');
add_bottom_label(fig, ax, 5, 'Session 3', 0.09, 11, 'normal');
add_bottom_label(fig, ax, 1.5, 'Onset', 0.01, 13, 'bold');
add_bottom_label(fig, ax, 4.5, 'Offset', 0.01, 13, 'bold');

ensure_parent(fullfile(outdir, 'dummy'));
exportgraphics(fig, fullfile(outdir, 'fig3_online_command_delivery_accuracy.png'), 'Resolution', 300);
exportgraphics(fig, fullfile(outdir, 'fig3_online_command_delivery_accuracy.pdf'), 'ContentType', 'vector');
close(fig);
end

function out = subject_session_props(T, label_col, outcomes, sessions)
subs = unique(T.subject);
rows = table();
for si = 1:numel(subs)
    for sj = 1:numel(sessions)
        mask = T.subject == subs(si) & T.session == sessions(sj);
        Ti = T(mask, :);
        if isempty(Ti)
            continue;
        end
        total = height(Ti);
        for oi = 1:numel(outcomes)
            count = sum(Ti.(label_col) == outcomes(oi));
            if count == 0
                continue;
            end
            rows = [rows; table(subs(si), sessions(sj), outcomes(oi), count / total, ...
                'VariableNames', {'subject', 'session', 'outcome', 'prop'})]; %#ok<AGROW>
        end
    end
end

G = groupsummary(rows, {'session', 'outcome'}, 'mean', 'prop');
out = table(G.session, G.outcome, G.mean_prop, 'VariableNames', {'session', 'outcome', 'mean_prop'});
end

function out = fill_missing_rows(summary, sessions, phases, outcomes)
session_col = [];
phase_col = strings(0,1);
outcome_col = strings(0,1);
mean_prop = [];
for si = 1:numel(sessions)
    for pi = 1:numel(phases)
        for oi = 1:numel(outcomes)
            mask = summary.session == sessions(si) & summary.phase == phases(pi) & summary.outcome == outcomes(oi);
            session_col(end + 1, 1) = sessions(si); %#ok<AGROW>
            phase_col(end + 1, 1) = phases(pi); %#ok<AGROW>
            outcome_col(end + 1, 1) = outcomes(oi); %#ok<AGROW>
            if any(mask)
                mean_prop(end + 1, 1) = summary.mean_prop(find(mask, 1)); %#ok<AGROW>
            else
                mean_prop(end + 1, 1) = 0.0; %#ok<AGROW>
            end
        end
    end
end
out = table(session_col, phase_col, outcome_col, mean_prop, ...
    'VariableNames', {'session', 'phase', 'outcome', 'mean_prop'});
end

function ensure_parent(pathstr)
folder = fileparts(pathstr);
if exist(folder, 'dir') ~= 7
    mkdir(folder);
end
end

function add_bottom_label(fig, ax, xval, label, y0, font_size, font_weight)
xnorm = ax.Position(1) + ax.Position(3) * ((xval - ax.XLim(1)) / diff(ax.XLim));
annotation(fig, 'textbox', [xnorm - 0.07, y0, 0.14, 0.05], 'String', label, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'LineStyle', 'none', 'FontSize', font_size, 'FontWeight', font_weight);
end
