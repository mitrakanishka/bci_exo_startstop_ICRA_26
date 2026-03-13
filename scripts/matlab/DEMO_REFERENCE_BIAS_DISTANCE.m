function DEMO_REFERENCE_BIAS_DISTANCE()
% GROUP_REFERENCE_BIAS_DISTANCE_ALL()
% Group-level diagnostic of unilateral, task-based recentering bias.
%
% OFFLINE negatives:
%   ONSET : pos = ONLINE Start-MI  (per run),  neg = OFFLINE REST (pooled)
%   OFFSET: pos = ONLINE Stop-MI   (per run),  neg = OFFLINE Doing-MI (pooled)
%
% For each subject and each run (1..8), we compute median s = d0 - d1 under:
%   - Identity (no recenter, R = I)
%   - Biased   (R = mean covariance of that run's pos windows)
% Then we compute Delta s = s_biased - s_identity for:
%   - pos (prompted class)
%   - neg (offline negative pool)
% and Delta separation: Delta_sep = Delta_pos - Delta_neg
%
% Outputs:
%   1) Console printouts:
%        - Subject list
%        - Per-run group mean +/- SD for Delta_pos, Delta_neg, Delta_sep
%        - OVERALL (per-subject mean across runs) group mean +/- SD
%        - Wilcoxon signed-rank tests (per run and overall) vs 0 for Delta_sep
%          and also for Delta_pos / Delta_neg
%
% Folder layout assumed:
%   BCI_Harmony_ExperimentalData/epoched_data_of/subj_###_epoched_data_of.mat
%   BCI_Harmony_ExperimentalData/epoched_data_on/subj_###_epoched_data_on.mat

clc;

repo_root = resolve_repo_root(fileparts(mfilename('fullpath')));
root     = fullfile(repo_root,'BCI_Harmony_ExperimentalData');
of_dir   = fullfile(root,'epoched_data_of');
on_dir   = fullfile(root,'epoched_data_on');
assert(exist(of_dir,'dir')==7 && exist(on_dir,'dir')==7, ...
    'Could not locate epoched data folders under %s', root);

% Discover overlapping subject IDs
of_files = dir(fullfile(of_dir, 'subj_*_epoched_data_of.mat'));
on_files = dir(fullfile(on_dir, 'subj_*_epoched_data_on.mat'));
sub_of   = parse_sub_ids(of_files);
sub_on   = parse_sub_ids(on_files);
subs     = intersect(sub_of, sub_on);
assert(~isempty(subs), 'No overlapping subjects found.');

fprintf('Found %d subjects: %s \n', numel(subs), sprintf('%03d ', subs));

S = numel(subs);
RMAX = 8;  % use runs 1..8 only

% Allocate: rows = subjects, cols = runs
Dpos_on  = nan(S, RMAX);  Dneg_on  = nan(S, RMAX);  Dsep_on  = nan(S, RMAX);
Dpos_off = nan(S, RMAX);  Dneg_off = nan(S, RMAX);  Dsep_off = nan(S, RMAX);

for si = 1:S
    sid = subs(si);
    Sof = load(fullfile(of_dir, sprintf('subj_%03d_epoched_data_of.mat', sid)));
    Son = load(fullfile(on_dir, sprintf('subj_%03d_epoched_data_on.mat', sid)));

    % Resolve online field names
    onStart = get_first_existing(Son, {'on_startMI_cell','on_bMI_cell'});
    onStop  = get_first_existing(Son, {'on_stopMI_cell','on_eMI_cell'});
    assert(~isempty(onStart) && ~isempty(onStop), 'Missing online epochs for subject %03d.', sid);

    % OFFLINE negative pools
    X_on_neg_pool  = pool_cell3d(Sof, 'data_rest_cell'); % REST
    X_off_neg_pool = pool_cell3d(Sof, 'data_dMI_cell');  % Doing-MI

    % OFFLINE prototypes (match training)
    % ONSET: P0 = REST, P1 = Begin-MI
    P0_on  = proto_from_cell(Sof.data_rest_cell);
    P1_on  = proto_from_cell(Sof.data_bMI_cell);
    % OFFSET: P0 = Doing-MI, P1 = End-MI
    P0_off = proto_from_cell(Sof.data_dMI_cell);
    P1_off = proto_from_cell(Sof.data_eMI_cell);

    nRunsAvail = max([size(onStart,1), size(onStop,1)]);
    nRuns = min(RMAX, nRunsAvail);

    for r = 1:nRuns
        % ONSET block (pos = online Start-MI of run r, neg = offline REST pooled)
        X_on_pos = pick3d(onStart, r);
        if ~isempty(X_on_pos) && ~isempty(X_on_neg_pool)
            % Identity
            [~, ~, s_pos_id] = distances_and_margins(X_on_pos,     P0_on, P1_on, eye(size(P0_on)));
            [~, ~, s_neg_id] = distances_and_margins(X_on_neg_pool, P0_on, P1_on, eye(size(P0_on)));
            % Biased (R = mean cov of run's Start-MI)
            R_on = ref_from_3d(X_on_pos);
            [~, ~, s_pos_bi] = distances_and_margins(X_on_pos,     P0_on, P1_on, R_on);
            [~, ~, s_neg_bi] = distances_and_margins(X_on_neg_pool, P0_on, P1_on, R_on);

            dpos = med_nan(s_pos_bi) - med_nan(s_pos_id);
            dneg = med_nan(s_neg_bi) - med_nan(s_neg_id);
            Dpos_on(si,r) = dpos;
            Dneg_on(si,r) = dneg;
            Dsep_on(si,r) = dpos - dneg;
        end

        % OFFSET block (pos = online Stop-MI of run r, neg = offline Doing-MI pooled)
        X_off_pos = pick3d(onStop, r);
        if ~isempty(X_off_pos) && ~isempty(X_off_neg_pool)
            % Identity
            [~, ~, s_pos_id] = distances_and_margins(X_off_pos,     P0_off, P1_off, eye(size(P0_off)));
            [~, ~, s_neg_id] = distances_and_margins(X_off_neg_pool, P0_off, P1_off, eye(size(P0_off)));
            % Biased (R = mean cov of run's Stop-MI)
            R_off = ref_from_3d(X_off_pos);
            [~, ~, s_pos_bi] = distances_and_margins(X_off_pos,     P0_off, P1_off, R_off);
            [~, ~, s_neg_bi] = distances_and_margins(X_off_neg_pool, P0_off, P1_off, R_off);

            dpos = med_nan(s_pos_bi) - med_nan(s_pos_id);
            dneg = med_nan(s_neg_bi) - med_nan(s_neg_id);
            Dpos_off(si,r) = dpos;
            Dneg_off(si,r) = dneg;
            Dsep_off(si,r) = dpos - dneg;
        end
    end
end

%% --------- Group summaries (per run) ----------
fprintf('\n=== Group Delta s (Biased - Identity) -- ONSET -- per run (mean +/- SD) ===\n');
print_run_means(Dpos_on,  'pos');
print_run_means(Dneg_on,  'neg');
fprintf('--- Delta separation (ONSET) ---\n');
print_run_means(Dsep_on,  'Delta_sep');

fprintf('\n=== Group Delta s (Biased - Identity) -- OFFSET -- per run (mean +/- SD) ===\n');
print_run_means(Dpos_off, 'pos');
print_run_means(Dneg_off, 'neg');
fprintf('--- Delta separation (OFFSET) ---\n');
print_run_means(Dsep_off, 'Delta_sep');

%% --------- OVERALL (per-subject mean across runs) ----------
sub_pos_on  = row_mean(Dpos_on);
sub_neg_on  = row_mean(Dneg_on);
sub_sep_on  = row_mean(Dsep_on);
sub_pos_off = row_mean(Dpos_off);
sub_neg_off = row_mean(Dneg_off);
sub_sep_off = row_mean(Dsep_off);

fprintf('\n=== OVERALL (per-subject mean across runs) -- group mean +/- SD ===\n');
fprintf('ONSET : Delta_pos = %+0.3f +/- %0.3f,  Delta_neg = %+0.3f +/- %0.3f,  Delta_sep = %+0.3f +/- %0.3f  (n=%d)\n', ...
    mean(sub_pos_on,'omitnan'),  std(sub_pos_on,'omitnan'), ...
    mean(sub_neg_on,'omitnan'),  std(sub_neg_on,'omitnan'), ...
    mean(sub_sep_on,'omitnan'),  std(sub_sep_on,'omitnan'),  sum(isfinite(sub_sep_on)));
fprintf('OFFSET: Delta_pos = %+0.3f +/- %0.3f,  Delta_neg = %+0.3f +/- %0.3f,  Delta_sep = %+0.3f +/- %0.3f  (n=%d)\n', ...
    mean(sub_pos_off,'omitnan'), std(sub_pos_off,'omitnan'), ...
    mean(sub_neg_off,'omitnan'), std(sub_neg_off,'omitnan'), ...
    mean(sub_sep_off,'omitnan'), std(sub_sep_off,'omitnan'), sum(isfinite(sub_sep_off)));

%% --------- Wilcoxon signed-rank tests (vs 0) ----------
% Per run tests for Delta_sep (primary effect), and also for components
fprintf('\n=== Wilcoxon signed-rank vs 0 (per run) ===\n');
fprintf('Delta_sep (ONSET):\n');
run_signrank(Dsep_on);
fprintf('Delta_sep (OFFSET):\n');
run_signrank(Dsep_off);

fprintf('Delta_pos (ONSET):\n');
run_signrank(Dpos_on);
fprintf('Delta_pos (OFFSET):\n');
run_signrank(Dpos_off);

fprintf('Delta_neg (ONSET):\n');
run_signrank(Dneg_on);
fprintf('Delta_neg (OFFSET):\n');
run_signrank(Dneg_off);

% Overall (per-subject mean across runs)
fprintf('\n=== Wilcoxon signed-rank vs 0 (overall; per-subject mean across runs) ===\n');
print_signrank('ONSET  Delta_sep',  sub_sep_on);
print_signrank('OFFSET Delta_sep',  sub_sep_off);
print_signrank('ONSET  Delta_pos',  sub_pos_on);
print_signrank('OFFSET Delta_pos',  sub_pos_off);
print_signrank('ONSET  Delta_neg',  sub_neg_on);
print_signrank('OFFSET Delta_neg',  sub_neg_off);

end

%% ======================= Helpers =======================

function ids = parse_sub_ids(dd)
ids = [];
for i = 1:numel(dd)
    m = regexp(dd(i).name,'subj_(\d+)_','tokens','once');
    if ~isempty(m), ids(end+1) = str2double(m{1}); end %#ok<AGROW>
end
ids = unique(ids);
end

function val = get_first_existing(S, names)
val = [];
for i = 1:numel(names)
    if isfield(S, names{i})
        val = S.(names{i});
        return;
    end
end
end

function X = pick3d(cellArr, r)
X = [];
if ~isempty(cellArr) && r <= size(cellArr,1)
    X = cellArr{r,1};
end
end

function Xpool = pool_cell3d(S, fieldname)
Xpool = [];
assert(isfield(S, fieldname), 'Missing field %s in struct.', fieldname);
cellArr = S.(fieldname);
for r = 1:size(cellArr,1)
    Xi = cellArr{r,1};
    if isempty(Xi), continue; end
    if isempty(Xpool), Xpool = Xi; else, Xpool = cat(3, Xpool, Xi); end
end
end

function P = proto_from_cell(cellArr)
% Log-Euclidean mean covariance across all windows in a cell array {run}[T x C x N]
Cstack = [];
for r = 1:size(cellArr,1)
    X = cellArr{r,1};
    if isempty(X), continue; end
    C = covs_from_3d(X);
    Cstack = cat_covs(Cstack, C);
end
if isempty(Cstack), error('Empty pool when building prototype.'); end
P = logeu_mean(Cstack);
end

function [d0, d1, s] = distances_and_margins(X3, P0, P1, R)
% For each window Xi in X3 (T x C x N):
% - Compute C = cov(Xi)
% - Whiten by R^{-1/2}
% - Affine-invariant distances to P0 and P1
% - s = d0 - d1
if isempty(X3)
    d0 = []; d1 = []; s = [];
    return;
end
P0 = ensure_spd(P0); P1 = ensure_spd(P1); R = ensure_spd(R);
[V,D] = eig((R+R')/2); lam = max(real(diag(D)), 1e-12);
Rminv = V * diag(1./sqrt(lam)) * V';

Cstack = covs_from_3d(X3);
n = size(Cstack,3);
d0 = zeros(n,1); d1 = zeros(n,1);
for i = 1:n
    Ci = ensure_spd(Cstack(:,:,i));
    Cw = (Rminv * Ci * Rminv'); Cw = (Cw + Cw')/2;
    d0(i) = spd_dist(Cw, P0);
    d1(i) = spd_dist(Cw, P1);
end
s = d0 - d1;
end

function C = covs_from_3d(X3)
% X3: [T x C x N] -> covariances (C x C x N)
if isempty(X3), C = []; return; end
[T, Cc, N] = size(X3);
C = zeros(Cc, Cc, N);
for i = 1:N
    Xi = X3(:,:,i);
    Ci = (Xi.'*Xi) ./ max(T-1,1);
    Ci = ensure_spd(Ci);
    C(:,:,i) = Ci;
end
end

function R = ref_from_3d(X3)
% Reference = log-Euclidean mean of covariances of X3
C = covs_from_3d(X3);
if isempty(C)
    R = eye(size(X3,2));
else
    R = logeu_mean(C);
end
end

function d = spd_dist(A, B)
% Affine-invariant distance via log-eigs of B^{-1/2} A B^{-1/2}
A = ensure_spd(A); B = ensure_spd(B);
[V,D] = eig((B+B')/2); lam = max(real(diag(D)), 1e-12);
Binv2 = V * diag(1./sqrt(lam)) * V';
M = (Binv2 * A * Binv2'); M = (M + M')/2;
ev = eig(M); ev = max(real(ev),1e-12);
d = sqrt(sum(log(ev).^2));
end

function P = logeu_mean(Cstack)
% Log-Euclidean mean of SPD stack
n = size(Cstack,3); Cdim = size(Cstack,1);
Lsum = zeros(Cdim);
for i = 1:n
    Ci = ensure_spd(Cstack(:,:,i));
    [V,D] = eig((Ci+Ci')/2); lam = max(real(diag(D)), 1e-12);
    Li = V * diag(log(lam)) * V';
    Lsum = Lsum + Li;
end
Lbar = Lsum / n;
[V,D] = eig((Lbar+Lbar')/2); lam = real(diag(D));
P = V * diag(exp(lam)) * V';
P = ensure_spd(P);
end

function C = cat_covs(A,B)
if isempty(A), C = B;
elseif isempty(B), C = A;
else, C = cat(3,A,B);
end
end

function A = ensure_spd(A)
A = (A + A')/2;
n = size(A,1);
epsw = 1e-9 * trace(A)/max(n,1);
if ~isfinite(epsw) || epsw<=0, epsw = 1e-9; end
A = A + epsw*eye(n);
end

function m = med_nan(x)
x = x(:); x = x(isfinite(x));
if isempty(x), m = NaN; else, m = median(x); end
end

function v = row_mean(M)
% mean across runs per subject
v = mean(M, 2, 'omitnan');
end

function print_run_means(M, label)
% M: S x R
mu = mean(M, 1, 'omitnan');
sd = std(M, 0, 1, 'omitnan');
fprintf('%-6s:', label);
for r = 1:min(8, numel(mu))
    if isfinite(mu(r))
        fprintf(' %+0.3f+/-%0.3f ', mu(r), sd(r));
    else
        fprintf('   NaN+/-NaN ');
    end
end
fprintf('\n');
end

function run_signrank(M)
% For each run: Wilcoxon signed-rank vs 0 across subjects
R = size(M,2);
for r = 1:min(8,R)
    x = M(:,r); x = x(isfinite(x));
    if numel(x) >= 2
        [p,~,stats] = signrank(x, 0, 'method','approximate');
        z = field_or(stats,'zval', NaN);
        fprintf('  Run %d: z=%+0.3f, p=%0.4f, n=%d\n', r, z, p, numel(x));
    else
        fprintf('  Run %d: insufficient data\n', r);
    end
end
end

function print_signrank(name, x)
x = x(isfinite(x));
if numel(x) >= 2
    [p,~,stats] = signrank(x, 0, 'method','approximate');
    z = field_or(stats,'zval', NaN);
    fprintf('  %s: z=%+0.3f, p=%0.4f, n=%d\n', name, z, p, numel(x));
else
    fprintf('  %s: insufficient data\n', name);
end
end

function v = field_or(S, name, defaultVal)
if isstruct(S) && isfield(S, name)
    v = S.(name);
else
    v = defaultVal;
end
end

function repo_root = resolve_repo_root(start_dir)
% Resolve repository root by walking upward until BCI_Harmony_ExperimentalData is found.
repo_root = start_dir;
for k = 1:8
    has_data = (exist(fullfile(repo_root,'BCI_Harmony_ExperimentalData'),'dir') == 7);
    if has_data
        return;
    end
    parent = fileparts(repo_root);
    if strcmp(parent, repo_root)
        break;
    end
    repo_root = parent;
end
error(['Could not locate repository root containing BCI_Harmony_ExperimentalData. ', ...
       'Run from the project or place data folder under repo root.']);
end
