function figure2_grand_avg_spectrogram(channel_name)
% Paper Figure 2: grand-average offline task spectrogram.

if nargin < 1 || isempty(channel_name)
    channel_name = 'C3';
end

clc;
script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(script_dir));
addpath(genpath(fullfile(repo_root, 'starter', 'functions')));

outdir = fullfile(repo_root, 'figures');
out_csv = fullfile(repo_root, 'fig_data', 'fig2_grand_avg_spectrogram_matrix.csv');
out_json = fullfile(repo_root, 'fig_data', 'fig2_grand_avg_spectrogram_meta.json');

dataset_root = fullfile(repo_root, 'BCI_Harmony_ExperimentalData');
if exist(dataset_root, 'dir') ~= 7
    dataset_root = fullfile(repo_root, 'BCI_course_EXP');
end
assert(exist(dataset_root, 'dir') == 7, 'Could not find the dataset folder.');
offline_root = fullfile(dataset_root, 'offline_data');
subs = dir(fullfile(offline_root, 'Sub_*'));
subs = subs([subs.isdir]);
assert(~isempty(subs), 'No offline subject folders found.');

subject_ids = [];
subject_mats = [];
kept_labels = {};
for i = 1:numel(subs)
    tok = regexp(subs(i).name, '^Sub_(\d+)$', 'tokens', 'once');
    if isempty(tok)
        continue;
    end
    sid = str2double(tok{1});
    [subj_mat, meta] = compute_subject_spectrogram(fullfile(offline_root, subs(i).name), upper(channel_name));
    subject_ids(end + 1, 1) = sid; %#ok<AGROW>
    if isempty(subject_mats)
        subject_mats = subj_mat;
    else
        subject_mats(:, :, end + 1) = subj_mat; %#ok<AGROW>
    end
    kept_labels = meta.kept_labels;
end
assert(~isempty(subject_ids), 'No valid offline spectrogram subjects were found.');

[subject_ids, order] = sort(subject_ids);
subject_mats = subject_mats(:, :, order);
grand_matrix = mean(subject_mats, 3);

meta = struct();
meta.subject_ids = subject_ids(:)';
meta.channel_name = upper(channel_name);
meta.num_windows = 16;
meta.freq_plot = 30:-1:8;
meta.kept_labels = kept_labels;
meta.n_subjects = numel(subject_ids);

ensure_parent(out_csv);
writematrix(grand_matrix, out_csv);
fid = fopen(out_json, 'w');
assert(fid ~= -1, 'Could not open %s for writing.', out_json);
fprintf(fid, '%s', jsonencode(meta, 'PrettyPrint', true));
fclose(fid);

fig = figure('Color', 'w', 'Position', [100 100 1000 750], 'Renderer', 'painters');
ax = axes(fig);
ax.Toolbar.Visible = 'off';
imagesc(ax, grand_matrix);
colormap(ax, 'jet');
cb = colorbar(ax);
cb.Label.String = 'Power (dB/Hz)';
cb.Label.FontSize = 12;

num_windows = double(meta.num_windows);
event_x = round([0.2 3.0 6.0 9.0 12.0 15.0 17.0] * num_windows);
event_labels = {'Countdown','Begin MI','Robot Moves','End MI','Robot Stops (rest)','Rest (moves)','Robot Returns'};
hold(ax, 'on');
for i = 1:numel(event_x)
    xline(ax, event_x(i), 'k-', 'LineWidth', 2.0);
    text(ax, event_x(i), -0.8, event_labels{i}, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 10, 'Clipping', 'off');
end

freq_plot = meta.freq_plot;
y_idx = 1:2:numel(freq_plot);
ax.XTick = event_x;
ax.XTickLabel = {'0','3','6','9','12','15','17'};
ax.YTick = y_idx;
ax.YTickLabel = string(freq_plot(y_idx));
ax.FontSize = 11;
xlabel(ax, 'Time (s)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel(ax, 'Frequency (Hz)', 'FontSize', 13, 'FontWeight', 'bold');
title(ax, sprintf('Grand Average Spectrogram, Electrode: %s', upper(channel_name)), ...
    'FontSize', 15, 'FontWeight', 'bold');
box(ax, 'off');
set(ax, 'TickDir', 'out');

ensure_parent(fullfile(outdir, 'dummy'));
exportgraphics(fig, fullfile(outdir, 'fig2_grand_avg_spectrogram.png'), 'Resolution', 300);
exportgraphics(fig, fullfile(outdir, 'fig2_grand_avg_spectrogram.pdf'), 'ContentType', 'vector');
close(fig);

cleanup_sopen(repo_root);
end

function [subject_mean, meta] = compute_subject_spectrogram(subject_dir, channel_name)
gdf_files = dir(fullfile(subject_dir, '*offline*', '*.gdf'));
assert(~isempty(gdf_files), 'No offline GDF files found in %s', subject_dir);
[~, order] = sort(fullfile(string({gdf_files.folder}), string({gdf_files.name})));
gdf_files = gdf_files(order);

keep_idx = [4:12 15:17 20:29];
subject_stack = [];
kept_labels = {};
for i = 1:numel(gdf_files)
    gdf_path = fullfile(gdf_files(i).folder, gdf_files(i).name);
    [signal, HDR] = sload(gdf_path);
    fs = double(HDR.SampleRate);
    stop_sample = double(HDR.EVENT.POS(end));
    data = double(signal(1:stop_sample, :));

    eeg = data(:, 1:64);
    eog = data(:, 65:min(67, size(data, 2)));
    [b, a] = butter(2, [0.1 45.0] / (fs / 2), 'bandpass');
    eeg = filtfilt(b, a, eeg);
    if ~isempty(eog)
        eog = filtfilt(b, a, eog);
        beta = eog \ eeg;
        eeg = eeg - eog * beta;
    end

    eeg = eeg(:, keep_idx);
    eeg = eeg - mean(eeg, 2);

    if isempty(kept_labels)
        labels = cellstr(strtrim(HDR.Label(keep_idx, :)));
        kept_labels = labels(:)';
    end
    spec_idx = find(strcmpi(kept_labels, channel_name), 1);
    assert(~isempty(spec_idx), 'Channel %s is not in the kept montage.', channel_name);

    countdown_pos = double(HDR.EVENT.POS(double(HDR.EVENT.TYP) == 300));
    countdown_pos = countdown_pos(countdown_pos - 2 * fs >= 1 & countdown_pos + 18 * fs <= size(eeg, 1));
    if isempty(countdown_pos)
        continue;
    end

    step_samples = round(0.0625 * fs);
    trial_len = round(18 / 0.0625);
    baseline_len = round((2 - 1) / 0.0625 + 1);
    freq_keep = 8:30;
    trial_psd = zeros(numel(freq_keep), trial_len, numel(countdown_pos));
    baseline_psd = zeros(numel(freq_keep), baseline_len, numel(countdown_pos));

    for ti = 1:numel(countdown_pos)
        curr = countdown_pos(ti);
        for wi = 1:trial_len
            [pxx, freq_axis] = pwelch(eeg(curr:(curr + fs - 1), spec_idx), hamming(round(0.5 * fs)), round(0.4 * fs), fs, fs);
            mask = freq_axis >= 8 & freq_axis <= 30;
            trial_psd(:, wi, ti) = pxx(mask);
            curr = curr + step_samples;
        end

        curr = countdown_pos(ti) - 2 * fs;
        for wi = 1:baseline_len
            [pxx, freq_axis] = pwelch(eeg(curr:(curr + fs - 1), spec_idx), hamming(round(0.5 * fs)), round(0.4 * fs), fs, fs);
            mask = freq_axis >= 8 & freq_axis <= 30;
            baseline_psd(:, wi, ti) = pxx(mask);
            curr = curr + step_samples;
        end
    end

    baseline_mean = mean(baseline_psd, 2);
    normalized = zeros(size(trial_psd));
    for ti = 1:size(trial_psd, 3)
        normalized(:, :, ti) = trial_psd(:, :, ti) ./ baseline_mean(:, 1, ti);
    end
    if isempty(subject_stack)
        subject_stack = flipud(mean(10 * log10(normalized), 3));
    else
        subject_stack(:, :, end + 1) = flipud(mean(10 * log10(normalized), 3)); %#ok<AGROW>
    end
end

assert(~isempty(subject_stack), 'No valid spectrogram trials found in %s', subject_dir);
subject_mean = mean(subject_stack, 3);
meta = struct('kept_labels', {kept_labels});
end

function ensure_parent(pathstr)
folder = fileparts(pathstr);
if exist(folder, 'dir') ~= 7
    mkdir(folder);
end
end

function cleanup_sopen(repo_root)
paths = {fullfile(repo_root, 'sopen.mat'), fullfile(repo_root, 'starter', 'sopen.mat')};
for i = 1:numel(paths)
    if exist(paths{i}, 'file') == 2
        delete(paths{i});
    end
end
end
