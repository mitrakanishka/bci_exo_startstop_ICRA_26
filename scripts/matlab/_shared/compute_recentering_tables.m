function [bias_tbl, auc_tbl] = compute_recentering_tables(of_dir, on_dir, max_runs, alpha, lamb, keep_frac)
% Recompute Figure 5/6 tables directly from the epoched offline/online data.

if nargin < 3 || isempty(max_runs), max_runs = 8; end
if nargin < 4 || isempty(alpha), alpha = 0.30; end
if nargin < 5 || isempty(lamb), lamb = 1e-3; end
if nargin < 6 || isempty(keep_frac), keep_frac = 0.90; end

phase_order = {"ONSET", "OFFSET"};
metric_order = {"Delta_pos", "Delta_neg", "Delta_sep"};
subjects = discover_overlapping_subjects(of_dir, on_dir);
Rrest_global = build_R_rest_global_from_dir(of_dir);

bias_subj = [];
bias_phase = strings(0,1);
bias_metric = strings(0,1);
bias_value = [];

auc_subj = [];
auc_run = [];
auc_phase = strings(0,1);
auc_scheme = strings(0,1);
auc_value = [];

for si = 1:numel(subjects)
    subject_id = subjects(si);
    offline_data = load(fullfile(of_dir, sprintf('subj_%03d_epoched_data_of.mat', subject_id)));
    online_data = load(fullfile(on_dir, sprintf('subj_%03d_epoched_data_on.mat', subject_id)));

    P0 = struct();
    P1 = struct();
    online_pos = struct();
    offline_neg = struct();
    for pi = 1:numel(phase_order)
        phase = phase_order{pi};
        if strcmp(phase, 'ONSET')
            P0.(phase) = logeu_mean(covs_from_cellarr(offline_data.data_rest_cell));
            P1.(phase) = logeu_mean(covs_from_cellarr(offline_data.data_bMI_cell));
            online_pos.(phase) = get_first_existing(online_data, {'on_startMI_cell','on_bMI_cell'});
            offline_neg.(phase) = pool_cell3d(offline_data, 'data_rest_cell');
        else
            P0.(phase) = logeu_mean(covs_from_cellarr(offline_data.data_dMI_cell));
            P1.(phase) = logeu_mean(covs_from_cellarr(offline_data.data_eMI_cell));
            online_pos.(phase) = get_first_existing(online_data, {'on_stopMI_cell','on_eMI_cell'});
            offline_neg.(phase) = pool_cell3d(offline_data, 'data_dMI_cell');
        end
    end

    acc = struct();
    for pi = 1:numel(phase_order)
        phase = phase_order{pi};
        for mi = 1:numel(metric_order)
            acc.(phase).(metric_order{mi}) = [];
        end
    end

    for run_index = 1:max_runs
        R_fix = build_R_fixation_robust(online_data, offline_data, run_index, Rrest_global, alpha, lamb, keep_frac);

        for pi = 1:numel(phase_order)
            phase = phase_order{pi};
            X_pos = pick3d(online_pos.(phase), run_index);
            X_neg = offline_neg.(phase);
            if isempty(X_pos) || isempty(X_neg)
                continue;
            end

            R_task = ref_from_3d(X_pos);
            I = eye(size(P0.(phase), 1));

            [~, ~, s_pos_id] = distances_and_margins(X_pos, P0.(phase), P1.(phase), I);
            [~, ~, s_neg_id] = distances_and_margins(X_neg, P0.(phase), P1.(phase), I);
            [~, ~, s_pos_task] = distances_and_margins(X_pos, P0.(phase), P1.(phase), R_task);
            [~, ~, s_neg_task] = distances_and_margins(X_neg, P0.(phase), P1.(phase), R_task);
            [~, ~, s_pos_fix] = distances_and_margins(X_pos, P0.(phase), P1.(phase), R_fix);
            [~, ~, s_neg_fix] = distances_and_margins(X_neg, P0.(phase), P1.(phase), R_fix);

            delta_pos = med_nan(s_pos_task) - med_nan(s_pos_id);
            delta_neg = med_nan(s_neg_task) - med_nan(s_neg_id);
            delta_sep = delta_pos - delta_neg;

            acc.(phase).Delta_pos(end + 1, 1) = delta_pos; %#ok<AGROW>
            acc.(phase).Delta_neg(end + 1, 1) = delta_neg; %#ok<AGROW>
            acc.(phase).Delta_sep(end + 1, 1) = delta_sep; %#ok<AGROW>

            auc_subj(end + 1, 1) = subject_id; %#ok<AGROW>
            auc_run(end + 1, 1) = run_index; %#ok<AGROW>
            auc_phase(end + 1, 1) = string(phase); %#ok<AGROW>
            auc_scheme(end + 1, 1) = "TASK"; %#ok<AGROW>
            auc_value(end + 1, 1) = auc_from_scores(s_pos_task, s_neg_task); %#ok<AGROW>

            auc_subj(end + 1, 1) = subject_id; %#ok<AGROW>
            auc_run(end + 1, 1) = run_index; %#ok<AGROW>
            auc_phase(end + 1, 1) = string(phase); %#ok<AGROW>
            auc_scheme(end + 1, 1) = "FIX"; %#ok<AGROW>
            auc_value(end + 1, 1) = auc_from_scores(s_pos_fix, s_neg_fix); %#ok<AGROW>
        end
    end

    for pi = 1:numel(phase_order)
        phase = phase_order{pi};
        for mi = 1:numel(metric_order)
            metric = metric_order{mi};
            vals = acc.(phase).(metric);
            bias_subj(end + 1, 1) = subject_id; %#ok<AGROW>
            bias_phase(end + 1, 1) = string(phase); %#ok<AGROW>
            bias_metric(end + 1, 1) = string(metric); %#ok<AGROW>
            if isempty(vals)
                bias_value(end + 1, 1) = NaN; %#ok<AGROW>
            else
                bias_value(end + 1, 1) = mean(vals, 'omitnan'); %#ok<AGROW>
            end
        end
    end
end

bias_tbl = table(bias_subj, bias_phase, bias_metric, bias_value, ...
    'VariableNames', {'subj', 'phase', 'metric', 'value'});
auc_tbl = table(auc_subj, auc_run, auc_phase, auc_scheme, auc_value, ...
    'VariableNames', {'subj', 'run', 'phase', 'scheme', 'auc'});

bias_tbl = sortrows(bias_tbl, {'subj', 'phase', 'metric'});
auc_tbl = sortrows(auc_tbl, {'subj', 'run', 'phase', 'scheme'});
end

function subjects = discover_overlapping_subjects(of_dir, on_dir)
of_files = dir(fullfile(of_dir, 'subj_*_epoched_data_of.mat'));
on_files = dir(fullfile(on_dir, 'subj_*_epoched_data_on.mat'));
subjects = intersect(parse_sub_ids(of_files), parse_sub_ids(on_files));
assert(~isempty(subjects), 'No overlapping epoched-data subjects found.');
end

function ids = parse_sub_ids(files)
ids = [];
for i = 1:numel(files)
    tok = regexp(files(i).name, 'subj_(\d+)_', 'tokens', 'once');
    if ~isempty(tok)
        ids(end + 1) = str2double(tok{1}); %#ok<AGROW>
    end
end
ids = unique(ids);
end

function value = get_first_existing(S, names)
value = [];
for i = 1:numel(names)
    if isfield(S, names{i})
        value = S.(names{i});
        return;
    end
end
end

function X = pick3d(cell_arr, run_index)
X = [];
if isempty(cell_arr) || run_index > size(cell_arr, 1)
    return;
end
X = to_epoch_tensor(cell_arr{run_index, 1});
end

function X = to_epoch_tensor(value)
if isempty(value)
    X = [];
    return;
end
X = double(value);
if ndims(X) == 2
    X = reshape(X, size(X, 1), size(X, 2), 1);
elseif ndims(X) ~= 3
    X = [];
end
end

function Xpool = pool_cell3d(S, field_name)
Xpool = [];
if ~isfield(S, field_name)
    return;
end
cell_arr = S.(field_name);
for i = 1:size(cell_arr, 1)
    Xi = to_epoch_tensor(cell_arr{i, 1});
    if isempty(Xi)
        continue;
    end
    if isempty(Xpool)
        Xpool = Xi;
    else
        Xpool = cat(3, Xpool, Xi);
    end
end
end

function Cstack = covs_from_cellarr(cell_arr)
Cstack = [];
for i = 1:size(cell_arr, 1)
    Xi = to_epoch_tensor(cell_arr{i, 1});
    if isempty(Xi)
        continue;
    end
    Cstack = cat_covs(Cstack, covs_from_3d(Xi));
end
end

function Cstack = covs_from_3d(X3)
if isempty(X3)
    Cstack = [];
    return;
end
[T, C, N] = size(X3);
Cstack = zeros(C, C, N);
denom = max(T - 1, 1);
for i = 1:N
    Xi = X3(:, :, i);
    Cstack(:, :, i) = ensure_spd((Xi' * Xi) / denom);
end
end

function C = cat_covs(A, B)
if isempty(A)
    C = B;
elseif isempty(B)
    C = A;
else
    C = cat(3, A, B);
end
end

function P = ref_from_3d(X3)
C = covs_from_3d(X3);
if isempty(C)
    P = eye(size(X3, 2));
else
    P = logeu_mean(C);
end
end

function [d0, d1, s] = distances_and_margins(X3, P0, P1, R)
if isempty(X3)
    d0 = [];
    d1 = [];
    s = [];
    return;
end
P0 = ensure_spd(P0);
P1 = ensure_spd(P1);
R = ensure_spd(R);
[V, D] = eig((R + R') / 2);
lam = max(real(diag(D)), 1e-12);
Rminv = V * diag(1 ./ sqrt(lam)) * V';

Cstack = covs_from_3d(X3);
n = size(Cstack, 3);
d0 = zeros(n, 1);
d1 = zeros(n, 1);
for i = 1:n
    Ci = ensure_spd(Cstack(:, :, i));
    Cw = ensure_spd(Rminv * Ci * Rminv');
    d0(i) = spd_dist(Cw, P0);
    d1(i) = spd_dist(Cw, P1);
end
s = d0 - d1;
end

function d = spd_dist(A, B)
A = ensure_spd(A);
B = ensure_spd(B);
[V, D] = eig(B);
lam = max(real(diag(D)), 1e-12);
Binv2 = V * diag(1 ./ sqrt(lam)) * V';
M = ensure_spd(Binv2 * A * Binv2');
ev = eig(M);
ev = max(real(ev), 1e-12);
d = sqrt(sum(log(ev).^2));
end

function L = spd_log(A)
A = ensure_spd(A);
[V, D] = eig(A);
lam = max(real(diag(D)), 1e-12);
L = V * diag(log(lam)) * V';
L = (L + L') / 2;
end

function A = spd_exp(L)
L = (L + L') / 2;
[V, D] = eig(L);
A = V * diag(exp(real(diag(D)))) * V';
A = ensure_spd(A);
end

function P = logeu_mean(Cstack)
assert(~isempty(Cstack), 'Cannot compute a prototype from an empty covariance stack.');
cdim = size(Cstack, 1);
Lsum = zeros(cdim, cdim);
for i = 1:size(Cstack, 3)
    Lsum = Lsum + spd_log(Cstack(:, :, i));
end
P = spd_exp(Lsum / size(Cstack, 3));
end

function A = ensure_spd(A)
A = double(A);
A = (A + A') / 2;
n = size(A, 1);
epsw = 1e-9 * trace(A) / max(n, 1);
if ~isfinite(epsw) || epsw <= 0
    epsw = 1e-9;
end
A = A + epsw * eye(n);
end

function value = med_nan(x)
x = double(x(:));
x = x(isfinite(x));
if isempty(x)
    value = NaN;
else
    value = median(x);
end
end

function auc = auc_from_scores(s_pos, s_neg)
s_pos = double(s_pos(:));
s_neg = double(s_neg(:));
s_pos = s_pos(isfinite(s_pos));
s_neg = s_neg(isfinite(s_neg));
n_pos = numel(s_pos);
n_neg = numel(s_neg);
if n_pos == 0 || n_neg == 0
    auc = NaN;
    return;
end
ranks = rankdata_average([s_pos; s_neg]);
ranks_pos = ranks(1:n_pos);
auc = (sum(ranks_pos) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg);
end

function ranks = rankdata_average(x)
[sorted, order] = sort(x);
ranks = zeros(size(x));
i = 1;
while i <= numel(sorted)
    j = i;
    while j < numel(sorted) && sorted(j + 1) == sorted(i)
        j = j + 1;
    end
    avg_rank = 0.5 * (i + j);
    ranks(order(i:j)) = avg_rank;
    i = j + 1;
end
end

function Rrest_global = build_R_rest_global_from_dir(of_dir)
Cstack = [];
cguess = [];
files = dir(fullfile(of_dir, 'subj_*_epoched_data_of.mat'));
for i = 1:numel(files)
    data = load(fullfile(files(i).folder, files(i).name));
    if ~isfield(data, 'data_rest_cell')
        continue;
    end
    Cstack = cat_covs(Cstack, covs_from_cellarr(data.data_rest_cell));
    if isempty(cguess)
        X0 = to_epoch_tensor(data.data_rest_cell{1, 1});
        if ~isempty(X0)
            cguess = size(X0, 2);
        end
    end
end
if isempty(Cstack)
    Rrest_global = eye(max(1, cguess));
else
    Rrest_global = logeu_mean(Cstack);
end
end

function X_fix = get_fix_for_run(online_data, offline_data, run_index)
fix_names = {'on_fix_cell','on_fixation_cell','on_preCue_cell','on_pre_fix_cell','on_rest_cell','on_countdown_cell'};
X_fix = [];
for i = 1:numel(fix_names)
    if isfield(online_data, fix_names{i})
        X_fix = pick3d(online_data.(fix_names{i}), run_index);
        if ~isempty(X_fix)
            return;
        end
    end
end
if isfield(offline_data, 'data_rest_cell')
    X_fix = pick3d(offline_data.data_rest_cell, run_index);
    if ~isempty(X_fix)
        return;
    end
    X_fix = pool_cell3d(offline_data, 'data_rest_cell');
end
end

function Rfix = build_R_fixation_robust(online_data, offline_data, run_index, Rrest_global, alpha, lamb, keep_frac)
X_fix = get_fix_for_run(online_data, offline_data, run_index);
Cfix = covs_from_3d(X_fix);
if isempty(Cfix)
    Rfix = ensure_spd(Rrest_global);
    return;
end

Mu_fix = logeu_mean(Cfix);
d = zeros(size(Cfix, 3), 1);
for i = 1:size(Cfix, 3)
    d(i) = spd_dist(Cfix(:, :, i), Mu_fix);
end
keep_n = max(1, round(keep_frac * size(Cfix, 3)));
[~, order] = sort(d);
Cfix_kept = Cfix(:, :, order(1:keep_n));
Rfix_run = logeu_mean(Cfix_kept);
L_mix = (1 - alpha) * spd_log(Rrest_global) + alpha * spd_log(Rfix_run);
Rfix = spd_exp(L_mix);
cdim = size(Rfix, 1);
Rfix = (1 - lamb) * Rfix + lamb * (trace(Rfix) / cdim) * eye(cdim);
Rfix = ensure_spd(Rfix);
end
