close all;
clear;
clc;

%% 1) Find the project and choose one offline file
% Run this script from either the repo root or from `starter/`.
% We first locate the project folders, add the shared MATLAB helpers,
% and pick the first offline GDF file as an example dataset.

[parent_dir, current_dir] = fileparts(pwd);
if strcmpi(current_dir, 'starter')
    repo_root = parent_dir;
else
    repo_root = pwd;
end

starter_dir = fullfile(repo_root, 'starter');
assert(exist(fullfile(starter_dir, 'functions'), 'dir') == 7, ...
    'Run this script from the project root or from starter/.');
addpath(genpath(fullfile(starter_dir, 'functions')));

dataset_root = fullfile(repo_root, 'BCI_Harmony_ExperimentalData');
if exist(dataset_root, 'dir') ~= 7
    dataset_root = fullfile(repo_root, 'BCI_course_EXP');
end
assert(exist(dataset_root, 'dir') == 7, 'Could not find the dataset folder.');

gdf_files = dir(fullfile(dataset_root, 'offline_data', 'Sub_*', '*offline*', '*.gdf'));
assert(~isempty(gdf_files), 'No offline .gdf files found.');
[~, order] = sort(fullfile(string({gdf_files.folder}), string({gdf_files.name})));
gdf_files = gdf_files(order);

target_file = fullfile(gdf_files(1).folder, gdf_files(1).name);
fprintf('Loading: %s\n', target_file);

%% 2) Open the GDF and read the basic metadata
% `sload` reads both the signal and the header. We keep the sample rate
% and channel labels because we need them for filtering and channel cleanup.
[signal, HDR] = sload(target_file);

if isfield(HDR, 'SampleRate')
    fs = double(HDR.SampleRate);
elseif isfield(HDR, 'SPR') && isfield(HDR, 'Dur') && HDR.Dur > 0
    fs = double(HDR.SPR / HDR.Dur);
else
    error('Could not infer sampling rate from the GDF header.');
end

if isfield(HDR, 'Label')
    if ischar(HDR.Label)
        chan_labels = strtrim(cellstr(HDR.Label));
    else
        chan_labels = cellstr(strtrim(string(HDR.Label(:))));
    end
else
    chan_labels = arrayfun(@(k) sprintf('Ch%d', k), 1:size(signal, 2), 'UniformOutput', false);
end

labels_upper = upper(string(chan_labels));
fprintf('Samples: %d | Channels: %d | Fs: %.2f Hz\n', size(signal, 1), size(signal, 2), fs);

%% 3) Keep EEG channels only
% We remove obvious non-EEG channels such as triggers and auxiliary sensors.
% `manual_remove` is the one place where you can quickly exclude channels
% that you do not want in the analysis.
manual_remove = ["Fp1", "Fp2"];

non_eeg_idx = find(any(contains(labels_upper, ["TRIG", "STATUS", "STI", "MARK", "AUX", "EMG", "ECG"]), 2));
eog_idx = find(any(contains(labels_upper, ["EOG", "HEOG", "VEOG"]), 2));
manual_idx = find(ismember(labels_upper, upper(manual_remove)));

eeg_idx = setdiff(1:size(signal, 2), unique([non_eeg_idx; eog_idx; manual_idx]));
assert(~isempty(eeg_idx), 'No EEG channels remain after channel removal.');

eeg_raw = signal(:, eeg_idx);
eeg_labels = chan_labels(eeg_idx);
fprintf('EEG channels kept: %d\n', numel(eeg_labels));

%% 4) Bandpass filter the EEG
% This is a simple 8-30 Hz motor-imagery style bandpass. We reuse your
% helper function so the filtering style matches the rest of the repo.
bp_hz = [8 30];
[b_bp, a_bp] = butter(2, bp_hz / (fs / 2), 'bandpass');
[eeg_bp, ~] = func_pseudo_online_filter(eeg_raw, struct('bpa', b_bp, 'bpb', a_bp), [], 'BPF');

%% 5) Regress out EOG activity when EOG channels exist
% `filterEOG` estimates regression weights from the EOG channels to the EEG.
% We then subtract that predicted eye-artifact contribution from the EEG.
if isempty(eog_idx)
    eeg_clean = eeg_bp;
    fprintf('No EOG channels detected; skipping EOG regression.\n');
else
    eog = signal(:, eog_idx);
    beta = filterEOG(eeg_bp, eog);
    eeg_clean = eeg_bp - eog * beta;
    fprintf('EOG regression applied using %d EOG channel(s).\n', numel(eog_idx));
end

%% 6) Make quick sanity-check plots
% Left: average power spectrum before/after preprocessing.
% Right: spectrogram of one cleaned EEG channel.
[pxx_raw, f] = pwelch(eeg_raw, round(2 * fs), [], [], fs);
[pxx_clean, ~] = pwelch(eeg_clean, round(2 * fs), [], [], fs);

figure('Color', 'w', 'Name', 'Offline preprocessing sanity check');

subplot(1, 2, 1);
plot(f, 10 * log10(mean(pxx_raw, 2) + eps), 'k', 'LineWidth', 1); hold on;
plot(f, 10 * log10(mean(pxx_clean, 2) + eps), 'b', 'LineWidth', 1.2);
xlim([0 60]);
grid on;
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('Raw EEG', 'Processed EEG', 'Location', 'northeast');
title('Average PSD');

subplot(1, 2, 2);
ch_plot = find(strcmpi(eeg_labels, 'C3'), 1);
if isempty(ch_plot)
    ch_plot = 1;
end
spectrogram(eeg_clean(:, ch_plot), round(2 * fs), round(1.5 * fs), 0:0.5:40, fs, 'yaxis');
ylim([0 40]);
title(sprintf('Spectrogram (%s)', eeg_labels{ch_plot}));

%% 7) Build the task spectrogram using your older pipeline style
% This cell mirrors the logic from your older spectrogram scripts:
% 1) use the first 64 channels as EEG and channels 65:67 as EOG regressors,
% 2) bandpass from 0.1 to 45 Hz,
% 3) keep the same reduced 22-channel montage,
% 4) apply CAR,
% 5) align trials to the countdown trigger,
% 6) compute a baseline-normalized PSD matrix, and
% 7) plot the result with the same event markers.
%
% Change `spectrogram_channel` if you want a different channel from the
% reduced montage. By default we plot C3.
spectrogram_channel = "C3";
spec_keep_idx = [4:12 15:17 20:29];
spec_labels = string(chan_labels(spec_keep_idx));
spec_ch_idx = find(strcmpi(spec_labels, spectrogram_channel), 1);
assert(~isempty(spec_ch_idx), 'Requested spectrogram channel was not found.');

spec_data = signal(1:double(HDR.EVENT.POS(end)), :);
spec_eeg = spec_data(:, 1:64);
spec_eog = spec_data(:, 65:min(67, size(spec_data, 2)));

[b_spec, a_spec] = butter(2, [0.1 45] / (fs / 2), 'bandpass');
spec_eeg = filtfilt(b_spec, a_spec, spec_eeg);
if ~isempty(spec_eog)
    spec_eog = filtfilt(b_spec, a_spec, spec_eog);
    spec_beta = filterEOG(spec_eeg, spec_eog);
    spec_eeg = spec_eeg - spec_eog * spec_beta;
end

spec_eeg = spec_eeg(:, spec_keep_idx);
spec_eeg = spec_eeg - mean(spec_eeg, 2);

countdown_pos = double(HDR.EVENT.POS(double(HDR.EVENT.TYP) == 300));
countdown_pos = countdown_pos(countdown_pos - 2 * fs >= 1 & countdown_pos + 18 * fs - 1 <= size(spec_eeg, 1));
assert(~isempty(countdown_pos), 'No valid countdown triggers were found for spectrogram trials.');

window_size = 1;
baseline_window_size = 2;
step_size = 0.0625;
num_windows = 1 / step_size;
num_trial_windows = 18 * num_windows;
num_baseline_windows = (baseline_window_size - window_size) / step_size + 1;
freq_keep = 8:30;

epochs_freq = zeros(numel(freq_keep), num_trial_windows, numel(countdown_pos));
baseline_epochs_freq = zeros(numel(freq_keep), num_baseline_windows, numel(countdown_pos));

for i_trial = 1:numel(countdown_pos)
    curr_pos = countdown_pos(i_trial);
    for a = 1:num_trial_windows
        [pxx, freq_axis] = pwelch(spec_eeg(curr_pos:curr_pos + fs - 1, spec_ch_idx), ...
            0.5 * fs, round(0.4 * fs), 4:40, fs);
        epochs_freq(:, a, i_trial) = pxx(freq_axis >= 8 & freq_axis <= 30);
        curr_pos = curr_pos + step_size * fs;
    end

    curr_pos = countdown_pos(i_trial) - 2 * fs;
    for a = 1:num_baseline_windows
        [pxx, freq_axis] = pwelch(spec_eeg(curr_pos:curr_pos + fs - 1, spec_ch_idx), ...
            0.5 * fs, round(0.4 * fs), 4:40, fs);
        baseline_epochs_freq(:, a, i_trial) = pxx(freq_axis >= 8 & freq_axis <= 30);
        curr_pos = curr_pos + step_size * fs;
    end
end

baseline_mean = reshape(mean(baseline_epochs_freq, 2), numel(freq_keep), []);
epochs_freq_bld = zeros(size(epochs_freq));
for i_trial = 1:size(epochs_freq, 3)
    epochs_freq_bld(:, :, i_trial) = epochs_freq(:, :, i_trial) ./ baseline_mean(:, i_trial);
end

PSD_matrix = flipud(mean(10 * log10(epochs_freq_bld), 3));
freq_plot = fliplr(freq_keep);

figure('Renderer', 'painters', 'Position', [10 10 1200 900], 'Color', 'w');
imagesc(PSD_matrix);
colormap('jet');
colorbar_handle = colorbar;
colorbar_handle.Label.String = 'Power (dB/Hz)';
clim([min(PSD_matrix(:)) max(PSD_matrix(:))]);

event_x = round([0.2 3.0 6.0 9.0 12.0 15.0 17.0] * num_windows);
event_labels = {'Countdown', 'Begin MI', 'Robot Moves', 'End MI', ...
    'Robot Stops (rest)', 'Rest (moves)', 'Robot Returns'};
marker_lines = xline(event_x, '-k', event_labels, 'FontSize', 12);
for i_line = 1:numel(marker_lines)
    marker_lines(i_line).LineWidth = 2;
end

xticks(event_x);
xticklabels({'0', '3', '6', '9', '12', '15', '17'});
yticks(1:2:numel(freq_plot));
yticklabels(freq_plot(1:2:end));
set(gca, 'FontSize', 14, 'Box', 'off');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title(sprintf('Task Spectrogram (%s)', spec_labels(spec_ch_idx)));

%% 8) Plot grand-average topoplots for each task condition
% Recent motor-imagery EEG papers usually show scalp maps of
% baseline-normalized mu/beta power (ERD/ERS). Here we use 8-30 Hz power
% in dB relative to the countdown part of each trial, average within each
% subject first, and then take the grand average across subjects.
% Negative values mean band-power suppression (ERD); positive values mean
% a power rebound (ERS). Because countdown is the reference, its map
% should stay close to 0 dB.
topo_band_hz = [8 30];
cond_names = {'Countdown', 'Start MI', 'Robot moves', 'Stop MI', 'Rest', 'Robot returns'};

all_labels_upper = upper(string(chan_labels(:)));
topo_drop = any(contains(all_labels_upper, ...
    ["TRIG", "STATUS", "STI", "MARK", "AUX", "EMG", "ECG", "EOG", "HEOG", "VEOG", "SENS"]), 2);
topo_labels = string(chan_labels(~topo_drop));

loc_hdr = struct();
loc_hdr.Label = cellstr(topo_labels(:));
loc_hdr.NS = numel(topo_labels);
loc_hdr = leadidcodexyz(loc_hdr);
loc_keep = all(~isnan(loc_hdr.ELEC.XYZ), 2);
topo_labels = topo_labels(loc_keep);
topo_xyz = loc_hdr.ELEC.XYZ(loc_keep, :);
topo_theta = -atan2d(topo_xyz(:, 2), topo_xyz(:, 1));
topo_radius = 0.5 - asind(topo_xyz(:, 3) ./ sqrt(sum(topo_xyz .^ 2, 2))) / 180;
topo_chanlocs = struct( ...
    'labels', cellstr(topo_labels(:)), ...
    'theta', num2cell(topo_theta), ...
    'radius', num2cell(topo_radius), ...
    'X', num2cell(topo_xyz(:, 1)), ...
    'Y', num2cell(topo_xyz(:, 2)), ...
    'Z', num2cell(topo_xyz(:, 3)));

[b_topo, a_topo] = butter(2, [0.1 45] / (fs / 2), 'bandpass');
[b_smr, a_smr] = butter(2, topo_band_hz / (fs / 2), 'bandpass');
subject_dirs = dir(fullfile(dataset_root, 'offline_data', 'Sub_*'));
subject_dirs = subject_dirs([subject_dirs.isdir]);
[~, order] = sort(string({subject_dirs.name}));
subject_dirs = subject_dirs(order);
subject_maps = nan(numel(topo_labels), numel(cond_names), numel(subject_dirs));

for i_subject = 1:numel(subject_dirs)
    run_files = dir(fullfile(subject_dirs(i_subject).folder, subject_dirs(i_subject).name, '*offline*', '*.gdf'));
    if isempty(run_files)
        continue;
    end

    [~, order] = sort(fullfile(string({run_files.folder}), string({run_files.name})));
    run_files = run_files(order);
    subject_sum = zeros(numel(topo_labels), numel(cond_names));
    subject_n = 0;

    for i_run = 1:numel(run_files)
        [run_signal, run_HDR] = sload(fullfile(run_files(i_run).folder, run_files(i_run).name));
        run_labels = strtrim(cellstr(run_HDR.Label));
        run_labels_upper = upper(string(run_labels(:)));
        [has_topo, topo_idx] = ismember(upper(topo_labels), run_labels_upper);
        if ~all(has_topo)
            continue;
        end

        run_eeg = double(run_signal(:, topo_idx));
        run_eeg = filtfilt(b_topo, a_topo, run_eeg);

        run_eog_idx = find(any(contains(run_labels_upper, ["EOG", "HEOG", "VEOG"]), 2));
        if ~isempty(run_eog_idx)
            run_eog = filtfilt(b_topo, a_topo, double(run_signal(:, run_eog_idx)));
            run_beta = filterEOG(run_eeg, run_eog);
            run_eeg = run_eeg - run_eog * run_beta;
        end

        run_eeg = run_eeg - mean(run_eeg, 2);
        run_smr = filtfilt(b_smr, a_smr, run_eeg);

        typ = double(run_HDR.EVENT.TYP(:));
        pos = double(run_HDR.EVENT.POS(:));
        p300 = pos(typ == 300);
        p100 = pos(typ == 100);
        p150 = pos(typ == 150);
        p500 = pos(typ == 500);
        p550 = pos(typ == 550);
        p900 = pos(typ == 900);
        p950 = pos(typ == 950);
        p2000 = pos(typ == 2000);
        n_trials = min([numel(p300), numel(p100), numel(p150), numel(p500), numel(p550), numel(p900), numel(p950)]);

        for i_trial = 1:n_trials
            trial_end = size(run_smr, 1);
            if i_trial < numel(p300)
                trial_end = min(trial_end, p300(i_trial + 1) - 1);
            elseif ~isempty(p2000)
                trial_end = min(trial_end, p2000(1) - 1);
            end

            stop_pos = min(p550(i_trial), p900(i_trial));
            segments = [
                p300(i_trial), p100(i_trial) - 1;
                p100(i_trial), p150(i_trial) - 1;
                p150(i_trial), p500(i_trial) - 1;
                p500(i_trial), stop_pos - 1;
                p900(i_trial), p950(i_trial) - 1;
                p950(i_trial), trial_end
            ];

            if any(segments(:, 1) < 1) || any(segments(:, 2) > size(run_smr, 1)) || any(diff(segments, 1, 2) < round(0.5 * fs))
                continue;
            end

            baseline_pow = mean(run_smr(segments(1, 1):segments(1, 2), :) .^ 2, 1) + eps;
            cond_map = zeros(numel(topo_labels), numel(cond_names));
            for i_cond = 1:numel(cond_names)
                seg = segments(i_cond, 1):segments(i_cond, 2);
                cond_pow = mean(run_smr(seg, :) .^ 2, 1) + eps;
                cond_map(:, i_cond) = 10 * log10(cond_pow ./ baseline_pow)';
            end

            subject_sum = subject_sum + cond_map;
            subject_n = subject_n + 1;
        end
    end

    if subject_n > 0
        subject_maps(:, :, i_subject) = subject_sum / subject_n;
    end
end

valid_subjects = squeeze(any(any(isfinite(subject_maps), 1), 2));
grand_maps = mean(subject_maps(:, :, valid_subjects), 3, 'omitnan');
clim_max = max(abs(grand_maps(:)));
blue_red = [linspace(0, 1, 128)', linspace(0, 1, 128)', ones(128, 1); ...
    ones(128, 1), linspace(1, 0, 128)', linspace(1, 0, 128)'];

figure('Color', 'w', 'Name', 'Grand-average SMR topoplots', 'Position', [50 50 1200 700]);
for i_cond = 1:numel(cond_names)
    subplot(2, 3, i_cond);
    topoplot(grand_maps(:, i_cond), topo_chanlocs, ...
        'maplimits', [-clim_max clim_max], ...
        'electrodes', 'off', ...
        'numcontour', 6);
    title(cond_names{i_cond}, 'FontSize', 14);
end
colormap(blue_red);
cb = colorbar('Position', [0.92 0.16 0.02 0.68]);
cb.Label.String = 'SMR power change (dB vs countdown)';
topo_title = annotation('textbox', [0.20 0.95 0.60 0.04], ...
    'String', sprintf('Grand-average topoplots (%d subjects, %d-%d Hz)', ...
    sum(valid_subjects), topo_band_hz(1), topo_band_hz(2)), ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'FontWeight', 'bold', ...
    'FontSize', 14);
