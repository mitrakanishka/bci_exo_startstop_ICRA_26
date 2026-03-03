function DEMO_TASK_vs_FIX_PERF()
% Compare Task-based vs Robust Fixation recentering under a frozen offline threshold (Identity-trained).
% Metrics: per-subject mean Balanced Accuracy (BA) across runs, operating-point drift, and ROC AUC.
% Paired Wilcoxon tests: TASK - FIX, reported for ONSET and OFFSET.

clc;
repo_root = resolve_repo_root(fileparts(mfilename('fullpath')));
root   = fullfile(repo_root,'BCI_course_EXP');
of_dir = fullfile(root,'epoched_data_of');
on_dir = fullfile(root,'epoched_data_on');
assert(exist(of_dir,'dir')==7 && exist(on_dir,'dir')==7, ...
    'Could not locate epoched data folders under %s', root);

fprintf('=== FIXED-THRESHOLD EVAL: TASK vs FIX (no params) ===\n');

% --- Discover overlapping subjects
of_files = dir(fullfile(of_dir,'subj_*_epoched_data_of.mat'));
on_files = dir(fullfile(on_dir,'subj_*_epoched_data_on.mat'));
subs = intersect(parse_sub_ids(of_files), parse_sub_ids(on_files));
assert(~isempty(subs), 'No overlapping subjects found.');
fprintf('Subjects: %s\n', sprintf('%03d ', subs));

% --- Build global REST reference for robust fixation
OPTS = struct('alpha',0.30, 'lambda',1e-3, 'keep_frac',0.90);
RREST_GLOBAL = build_R_rest_global_from_dir(of_dir);

% --- Output containers (rows: subjects)   [:,1]=TASK, [:,2]=FIX
nS = numel(subs);
BA_on  = nan(nS,2);  DR_on  = nan(nS,2);  AUC_on  = nan(nS,2);
BA_off = nan(nS,2);  DR_off = nan(nS,2);  AUC_off = nan(nS,2);

for si = 1:nS
    sid = subs(si);
    Sof = load(fullfile(of_dir, sprintf('subj_%03d_epoched_data_of.mat', sid)));
    Son = load(fullfile(on_dir, sprintf('subj_%03d_epoched_data_on.mat', sid)));

    % --- Prototypes
    P0_on  = proto_from_cell(Sof.data_rest_cell);  P1_on  = proto_from_cell(Sof.data_bMI_cell);
    P0_off = proto_from_cell(Sof.data_dMI_cell);   P1_off = proto_from_cell(Sof.data_eMI_cell);

    % --- Frozen thresholds τ from OFFLINE Identity (fair to both schemes)
    tau_on  = pick_threshold( ...
        margins_from_cell(Sof.data_rest_cell, P0_on,  P1_on), ...
        margins_from_cell(Sof.data_bMI_cell,  P0_on,  P1_on) );
    tau_off = pick_threshold( ...
        margins_from_cell(Sof.data_dMI_cell,  P0_off, P1_off), ...
        margins_from_cell(Sof.data_eMI_cell,  P0_off, P1_off) );

    % --- Online positives & negative pools
    onStart = get_first_existing(Son, {'on_startMI_cell','on_bMI_cell'});   % ONSET positives (per run)
    onStop  = get_first_existing(Son, {'on_stopMI_cell','on_eMI_cell'});    % OFFSET positives (per run)
    X_on_neg_pool  = pool_cell3d(Sof, 'data_rest_cell'); % ONSET negatives
    X_off_neg_pool = pool_cell3d(Sof, 'data_dMI_cell');  % OFFSET negatives

    RMAX = 8; nRuns = min(RMAX, max([size(onStart,1), size(onStop,1)]));

    % Per-run accumulators
    BA_on_task=[]; BA_on_fix=[]; DR_on_task=[]; DR_on_fix=[]; AUC_on_task=[]; AUC_on_fix=[];
    BA_off_task=[];BA_off_fix=[];DR_off_task=[];DR_off_fix=[];AUC_off_task=[];AUC_off_fix=[];

    for r = 1:nRuns
        % ----- Build references
        X_on_pos  = pick3d(onStart, r);
        X_off_pos = pick3d(onStop,  r);

        % Task-based references (from prompted class of THIS run)
        R_task_on  = ref_from_3d(X_on_pos);
        R_task_off = ref_from_3d(X_off_pos);
        % Robust Fixation-based
        R_fix = build_R_fixation_robust(Son, Sof, r, RREST_GLOBAL, OPTS);

        % ----- ONSET: pos = online Start-MI (run r), neg = offline REST pool
        if ~isempty(X_on_pos) && ~isempty(X_on_neg_pool)
            % Task margins
            [~,~, sp_t] = distances_and_margins(X_on_pos,      P0_on,P1_on, R_task_on);
            [~,~, sn_t] = distances_and_margins(X_on_neg_pool,  P0_on,P1_on, R_task_on);
            % Fix margins
            [~,~, sp_f] = distances_and_margins(X_on_pos,      P0_on,P1_on, R_fix);
            [~,~, sn_f] = distances_and_margins(X_on_neg_pool,  P0_on,P1_on, R_fix);

            % BA @ fixed τ (predict pos if s<τ)
            BA_on_task(end+1) = balanced_acc(sp_t, sn_t, tau_on);
            BA_on_fix(end+1)  = balanced_acc(sp_f, sn_f, tau_on);
            % Drift relative to τ (sum of |median shifts|)
            DR_on_task(end+1) = abs(med_nan(sp_t)-tau_on) + abs(med_nan(sn_t)-tau_on);
            DR_on_fix(end+1)  = abs(med_nan(sp_f)-tau_on) + abs(med_nan(sn_f)-tau_on);
            % Threshold-free AUC (lower s = more positive)
            AUC_on_task(end+1)= auc_from_scores(sp_t, sn_t);
            AUC_on_fix(end+1) = auc_from_scores(sp_f, sn_f);
        end

        % ----- OFFSET: pos = online Stop-MI (run r), neg = offline Doing-MI pool
        if ~isempty(X_off_pos) && ~isempty(X_off_neg_pool)
            % Task margins
            [~,~, sp_t] = distances_and_margins(X_off_pos,      P0_off,P1_off, R_task_off);
            [~,~, sn_t] = distances_and_margins(X_off_neg_pool,  P0_off,P1_off, R_task_off);
            % Fix margins
            [~,~, sp_f] = distances_and_margins(X_off_pos,      P0_off,P1_off, R_fix);
            [~,~, sn_f] = distances_and_margins(X_off_neg_pool,  P0_off,P1_off, R_fix);

            BA_off_task(end+1) = balanced_acc(sp_t, sn_t, tau_off);
            BA_off_fix(end+1)  = balanced_acc(sp_f, sn_f, tau_off);

            DR_off_task(end+1) = abs(med_nan(sp_t)-tau_off) + abs(med_nan(sn_t)-tau_off);
            DR_off_fix(end+1)  = abs(med_nan(sp_f)-tau_off) + abs(med_nan(sn_f)-tau_off);

            AUC_off_task(end+1)= auc_from_scores(sp_t, sn_t);
            AUC_off_fix(end+1) = auc_from_scores(sp_f, sn_f);
        end
    end

    % Per-subject means across runs
    BA_on(si,:)  = [mean(BA_on_task,'omitnan'),  mean(BA_on_fix,'omitnan')];
    DR_on(si,:)  = [mean(DR_on_task,'omitnan'),  mean(DR_on_fix,'omitnan')];
    AUC_on(si,:) = [mean(AUC_on_task,'omitnan'), mean(AUC_on_fix,'omitnan')];

    BA_off(si,:)  = [mean(BA_off_task,'omitnan'),  mean(BA_off_fix,'omitnan')];
    DR_off(si,:)  = [mean(DR_off_task,'omitnan'),  mean(DR_off_fix,'omitnan')];
    AUC_off(si,:) = [mean(AUC_off_task,'omitnan'), mean(AUC_off_fix,'omitnan')];
end

% --- Group summaries
prt = @(name,x) fprintf('%-28s %0.3f +/- %0.3f\n', name, mean(x,'omitnan'), std(x,0,'omitnan'));

fprintf('\nONSET Balanced Accuracy (TASK vs FIX)\n');  prt('TASK', BA_on(:,1));  prt('FIX', BA_on(:,2));
fprintf('ONSET ROC AUC (TASK vs FIX)\n');            prt('TASK', AUC_on(:,1)); prt('FIX', AUC_on(:,2));
fprintf('ONSET Operating-point drift\n');            prt('TASK', DR_on(:,1));  prt('FIX', DR_on(:,2));
paired_signrank('ONSET  (TASK - FIX) BA',    BA_on(:,1)  - BA_on(:,2));
paired_signrank('ONSET  (TASK - FIX) AUC',   AUC_on(:,1) - AUC_on(:,2));
paired_signrank('ONSET  (TASK - FIX) Drift', DR_on(:,1)  - DR_on(:,2));

fprintf('\nOFFSET Balanced Accuracy (TASK vs FIX)\n'); prt('TASK', BA_off(:,1)); prt('FIX', BA_off(:,2));
fprintf('OFFSET ROC AUC (TASK vs FIX)\n');           prt('TASK', AUC_off(:,1)); prt('FIX', AUC_off(:,2));
fprintf('OFFSET Operating-point drift\n');           prt('TASK', DR_off(:,1));  prt('FIX', DR_off(:,2));
paired_signrank('OFFSET (TASK - FIX) BA',    BA_off(:,1) - BA_off(:,2));
paired_signrank('OFFSET (TASK - FIX) AUC',   AUC_off(:,1)- AUC_off(:,2));
paired_signrank('OFFSET (TASK - FIX) Drift', DR_off(:,1) - DR_off(:,2));

fprintf('\nDone.\n');
end

%% ======================= Helpers (same conventions as your pipeline) =======================
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
C = covs_from_3d(X3);
if isempty(C)
    R = eye(size(X3,2));
else
    R = logeu_mean(C);
end
end

function d = spd_dist(A, B)
A = ensure_spd(A); B = ensure_spd(B);
[V,D] = eig((B+B')/2); lam = max(real(diag(D)), 1e-12);
Binv2 = V * diag(1./sqrt(lam)) * V';
M = (Binv2 * A * Binv2'); M = (M + M')/2;
ev = eig(M); ev = max(real(ev),1e-12);
d = sqrt(sum(log(ev).^2));
end

function P = logeu_mean(Cstack)
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

function s = margins_from_cell(cellArr, P0, P1)
C = covs_from_cellarr(cellArr); if isempty(C), s = []; return; end
d0 = zeros(size(C,3),1); d1 = d0;
for i=1:numel(d0), d0(i)=spd_dist(C(:,:,i),P0); d1(i)=spd_dist(C(:,:,i),P1); end
s = d0 - d1;
end

function Cstack = covs_from_cellarr(cellArr)
Cstack = [];
for r = 1:size(cellArr,1)
    X = cellArr{r,1}; if isempty(X), continue; end
    Cstack = cat_covs(Cstack, covs_from_3d(X));
end
end

function tau = pick_threshold(s0, s1)
% Grid on quantiles to maximize balanced accuracy
x = [s0(:); s1(:)];
if isempty(x) || all(~isfinite(x)), tau = 0; return; end
th = quantile(x(isfinite(x)), linspace(0.01,0.99,101));
BA = arrayfun(@(t) balanced_acc(s1,s0,t), th);
[~,ix] = max(BA); tau = th(ix);
end

function ba = balanced_acc(s_pos, s_neg, tau)
% predict pos if s<tau
TPR = mean(s_pos < tau, 'omitnan');
TNR = mean(s_neg >= tau, 'omitnan');
ba = 0.5*(TPR + TNR);
end

function auc = auc_from_scores(s_pos, s_neg)
% AUC = P(s_pos < s_neg). Uses rank-based Mann–Whitney formulation.
sp = s_pos(isfinite(s_pos)); sn = s_neg(isfinite(s_neg));
np = numel(sp); nn = numel(sn);
if np==0 || nn==0, auc = NaN; return; end
all_scores = [sp(:); sn(:)];
[~,~,rk] = unique(all_scores); rk = double(rk);  % average ranks not needed with unique()
rp = rk(1:np);
auc = (sum(rp) - np*(np+1)/2) / (np*nn);
end

function v = field_or(S, name, defaultVal)
if isstruct(S) && isfield(S, name)
    v = S.(name);
else
    v = defaultVal;
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

function paired_signrank(name, diff_vec)
x = diff_vec(isfinite(diff_vec));
if numel(x) >= 2
    [p,~,stats] = signrank(x, 0, 'method','approximate');
    z = field_or(stats,'zval', NaN);
    fprintf('  %s: z=%+0.3f, p=%0.4f, n=%d\n', name, z, p, numel(x));
else
    fprintf('  %s: insufficient data\n', name);
end
end

function X_fix = get_fix_for_run(Son, Sof, r)
% Try to find neutral/fixation windows for this run from the online struct.
% Fallbacks: offline REST of the same run; else the pooled offline REST.
X_fix = [];
candidates = {'on_fix_cell','on_fixation_cell','on_preCue_cell','on_pre_fix_cell','on_rest_cell','on_countdown_cell'};
for k = 1:numel(candidates)
    if isfield(Son, candidates{k})
        X_fix = pick3d(Son.(candidates{k}), r);
        if ~isempty(X_fix), break; end
    end
end
if isempty(X_fix) && isfield(Sof,'data_rest_cell')
    X_fix = pick3d(Sof.data_rest_cell, r);  % per-run offline REST if exists
end
if isempty(X_fix)
    X_fix = pool_cell3d(Sof, 'data_rest_cell');  % global offline REST pool
end
end

function Rrest_global = build_R_rest_global_from_dir(of_dir)
% Build a single global REST reference across all offline subjects/runs.
files = dir(fullfile(of_dir, 'subj_*_epoched_data_of.mat'));
Cstack = []; Cguess = [];
for k = 1:numel(files)
    S = load(fullfile(of_dir, files(k).name));
    if isfield(S,'data_rest_cell') && ~isempty(S.data_rest_cell)
        Cstack = cat_covs(Cstack, covs_from_cellarr(S.data_rest_cell));
        if isempty(Cguess)
            X0 = S.data_rest_cell{1,1};
            if ~isempty(X0), Cguess = size(X0,2); end
        end
    end
end
if isempty(Cstack)
    if isempty(Cguess), Cguess = 32; end
    warning('Global REST pool empty; using identity(%d).', Cguess);
    Rrest_global = eye(Cguess);
else
    Rrest_global = logeu_mean(Cstack);
end
end

function Rfix = build_R_fixation_robust(Son, Sof, r, Rrest_global, opts)
% Robust, label-agnostic fixation reference for run r.
% - neutral windows via get_fix_for_run(...)
% - outlier trim by SPD distance (keep_frac)
% - shrinkage to global REST in log-Euclidean space (alpha)
% - diagonal loading for stability (lambda)
if ~isfield(opts,'alpha'),     opts.alpha     = 0.30; end
if ~isfield(opts,'lambda'),    opts.lambda    = 1e-3; end
if ~isfield(opts,'keep_frac'), opts.keep_frac = 0.90; end

X_fix = get_fix_for_run(Son, Sof, r);
Cfix  = covs_from_3d(X_fix);
if isempty(Cfix)
    Rfix = ensure_spd(Rrest_global);
    return;
end

% Outlier trimming
Mu_fix = logeu_mean(Cfix);
nC = size(Cfix,3); d = zeros(nC,1);
for i = 1:nC
    d(i) = spd_dist(Cfix(:,:,i), Mu_fix);
end
[~,ord] = sort(d,'ascend');
keep_n  = max(1, round(opts.keep_frac * nC));
Cfix_kept = Cfix(:,:,ord(1:keep_n));

% Shrinkage in log-Euclidean space
Rfix_run = logeu_mean(Cfix_kept);
L_rest   = spd_log(Rrest_global);
L_fix    = spd_log(Rfix_run);
L_mix    = (1-opts.alpha)*L_rest + opts.alpha*L_fix;
Rfix     = spd_exp(L_mix);

% Diagonal loading & SPD cleanup
Cdim = size(Rfix,1);
Rfix = (1-opts.lambda)*Rfix + opts.lambda*(trace(Rfix)/Cdim)*eye(Cdim);
Rfix = ensure_spd(Rfix);
end

function L = spd_log(A)
A = ensure_spd(A); A = (A + A')/2;
[V,D] = eig(A); lam = max(real(diag(D)), 1e-12);
L = V * diag(log(lam)) * V';
L = (L + L')/2;
end

function A = spd_exp(L)
L = (L + L')/2;
[V,D] = eig(L); lam = real(diag(D));
A = V * diag(exp(lam)) * V';
A = ensure_spd(A);
end

function repo_root = resolve_repo_root(start_dir)
% Resolve repository root by walking upward until BCI_course_EXP is found.
repo_root = start_dir;
for k = 1:8
    has_data = (exist(fullfile(repo_root,'BCI_course_EXP'),'dir') == 7);
    if has_data
        return;
    end
    parent = fileparts(repo_root);
    if strcmp(parent, repo_root)
        break;
    end
    repo_root = parent;
end
error(['Could not locate repository root containing BCI_course_EXP. ', ...
       'Run from the project or place data folder under repo root.']);
end
