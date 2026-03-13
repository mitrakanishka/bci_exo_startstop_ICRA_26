function run_all_figures()
% Run all paper figure scripts from MATLAB.

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir, fullfile(script_dir, '_shared'));

figure_fns = {@figure2_grand_avg_spectrogram, @figure3_online_command_delivery, @figure4_online_decoding_time, @figure5_bias_shift_vs_identity, @figure6_auc_by_run_task_vs_fix};
figure_names = {'figure2_grand_avg_spectrogram', 'figure3_online_command_delivery', 'figure4_online_decoding_time', 'figure5_bias_shift_vs_identity', 'figure6_auc_by_run_task_vs_fix'};

for i = 1:numel(figure_fns)
    fprintf('\n=== Running %s ===\n', figure_names{i});
    figure_fns{i}();
end
end
