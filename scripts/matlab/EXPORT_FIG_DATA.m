function EXPORT_FIG_DATA()
% EXPORT_FIG_DATA
% Extracts all values needed for:
%  - Figure 1 (bias bars: Δpos, Δneg, Δsep for TASK vs Identity)
%  - Figure 2 (ROC traces for TASK vs FIX; AUCs; Δsep inset for TASK & FIX)
% Outputs CSVs under ./fig_data/
%
% Assumed layout:
%   BCI_course_EXP/
%     epoched_data_of/subj_###_epoched_data_of.mat
%     epoched_data_on/subj_###_epoched_data_on.mat

clc; fprintf('Exporting fig data...\n');

repo_root = resolve_repo_root(fileparts(mfilename('fullpath')));
root   = fullfile(repo_root,'BCI_course_EXP');
of_dir = fullfile(root,'epoched_data_of');
on_dir = fullfile(root,'epoched_data_on');
assert(exist(of_dir,'dir')==7 && exist(on_dir,'dir')==7, ...
    'Could not locate epoched data folders under %s', root);
outdir = fullfile(repo_root,'fig_data'); if ~exist(outdir,'dir'), mkdir(outdir); end

% --- discover overlapping subjects
of_files = dir(fullfile(of_dir,'subj_*_epoched_data_of.mat'));
on_files = dir(fullfile(on_dir,'subj_*_epoched_data_on.mat'));
subs = intersect(parse_sub_ids(of_files), parse_sub_ids(on_files));
assert(~isempty(subs),'No overlapping subjects found.');
fprintf('Subjects: %s\n', sprintf('%03d ', subs));

% --- robust fixation opts + global REST
OPTS = struct('alpha',0.30,'lambda',1e-3,'keep_frac',0.90);
RREST_GLOBAL = build_R_rest_global_from_dir(of_dir);

% --- outputs
T_bias = table();  % Fig1: per-subject Δpos/Δneg/Δsep (TASK vs Identity), means across runs
T_roc  = table();  % Fig2: per (subject,run,phase,scheme,fpr)->tpr
T_auc  = table();  % Fig2: per (subject,run,phase,scheme) AUC
T_dsep = table();  % Fig2 inset: per-subject Δsep for TASK & FIX vs Identity (means across runs)

FPR_GRID = linspace(0,1,201);

for si = 1:numel(subs)
    sid = subs(si);
    Sof = load(fullfile(of_dir, sprintf('subj_%03d_epoched_data_of.mat', sid)));
    Son = load(fullfile(on_dir, sprintf('subj_%03d_epoched_data_on.mat', sid)));

    % offline prototypes
    P0_on  = proto_from_cell(Sof.data_rest_cell);
    P1_on  = proto_from_cell(Sof.data_bMI_cell);
    P0_off = proto_from_cell(Sof.data_dMI_cell);
    P1_off = proto_from_cell(Sof.data_eMI_cell);

    % online positives & offline negative pools
    onStart = get_first_existing(Son, {'on_startMI_cell','on_bMI_cell'});
    onStop  = get_first_existing(Son, {'on_stopMI_cell','on_eMI_cell'});
    X_on_neg_pool  = pool_cell3d(Sof, 'data_rest_cell'); % ONSET negatives
    X_off_neg_pool = pool_cell3d(Sof, 'data_dMI_cell');  % OFFSET negatives

    RMAX = 8; nRuns = min(RMAX, max([size(onStart,1), size(onStop,1)]));

    % accumulators for subject-level means across runs
    acc.on.delta_pos = []; acc.on.delta_neg = []; acc.on.delta_sep = [];
    acc.off.delta_pos= []; acc.off.delta_neg= []; acc.off.delta_sep= [];
    acc.on.dsep_task = []; acc.on.dsep_fix  = [];
    acc.off.dsep_task= []; acc.off.dsep_fix = [];

    for r = 1:nRuns
        % references
        X_on_pos  = pick3d(onStart, r);
        X_off_pos = pick3d(onStop,  r);

        R_task_on  = ref_from_3d(X_on_pos);
        R_task_off = ref_from_3d(X_off_pos);
        R_fix      = build_R_fixation_robust(Son, Sof, r, RREST_GLOBAL, OPTS);

        % ===== ONSET =====
        if ~isempty(X_on_pos) && ~isempty(X_on_neg_pool)
            % Identity
            [~,~, sp_id] = distances_and_margins(X_on_pos,      P0_on,P1_on, eye(size(P0_on)));
            [~,~, sn_id] = distances_and_margins(X_on_neg_pool,  P0_on,P1_on, eye(size(P0_on)));
            % Task
            [~,~, sp_t]  = distances_and_margins(X_on_pos,      P0_on,P1_on, R_task_on);
            [~,~, sn_t]  = distances_and_margins(X_on_neg_pool,  P0_on,P1_on, R_task_on);
            % Fix
            [~,~, sp_f]  = distances_and_margins(X_on_pos,      P0_on,P1_on, R_fix);
            [~,~, sn_f]  = distances_and_margins(X_on_neg_pool,  P0_on,P1_on, R_fix);

            % --- deltas (Task vs Identity) → Fig1
            dpos = med_nan(sp_t) - med_nan(sp_id);
            dneg = med_nan(sn_t) - med_nan(sn_id);
            dsep = dpos - dneg;
            acc.on.delta_pos(end+1) = dpos;
            acc.on.delta_neg(end+1) = dneg;
            acc.on.delta_sep(end+1) = dsep;

            % --- Δsep (Task & Fix vs Identity) → Fig2 inset
            dsep_task = (med_nan(sp_t)-med_nan(sp_id)) - (med_nan(sn_t)-med_nan(sn_id));
            dsep_fix  = (med_nan(sp_f)-med_nan(sp_id)) - (med_nan(sn_f)-med_nan(sn_id));
            acc.on.dsep_task(end+1) = dsep_task;
            acc.on.dsep_fix(end+1)  = dsep_fix;

            % --- ROC traces (per run) → Fig2 main
            [fpr_t, tpr_t] = roc_on_grid(sp_t, sn_t, FPR_GRID);
            [fpr_f, tpr_f] = roc_on_grid(sp_f, sn_f, FPR_GRID);
            T_roc = add_rows(T_roc, sid, r, 'ONSET','TASK', fpr_t, tpr_t);
            T_roc = add_rows(T_roc, sid, r, 'ONSET','FIX',  fpr_f, tpr_f);

            % --- AUC per subject-run → Figure 2 stats
            T_auc = add_auc_row(T_auc, sid, r, 'ONSET','TASK', auc_from_scores(sp_t, sn_t));
            T_auc = add_auc_row(T_auc, sid, r, 'ONSET','FIX',  auc_from_scores(sp_f, sn_f));
        end

        % ===== OFFSET =====
        if ~isempty(X_off_pos) && ~isempty(X_off_neg_pool)
            % Identity
            [~,~, sp_id] = distances_and_margins(X_off_pos,      P0_off,P1_off, eye(size(P0_off)));
            [~,~, sn_id] = distances_and_margins(X_off_neg_pool,  P0_off,P1_off, eye(size(P0_off)));
            % Task
            [~,~, sp_t]  = distances_and_margins(X_off_pos,      P0_off,P1_off, R_task_off);
            [~,~, sn_t]  = distances_and_margins(X_off_neg_pool,  P0_off,P1_off, R_task_off);
            % Fix
            [~,~, sp_f]  = distances_and_margins(X_off_pos,      P0_off,P1_off, R_fix);
            [~,~, sn_f]  = distances_and_margins(X_off_neg_pool,  P0_off,P1_off, R_fix);

            % --- deltas (Task vs Identity) → Fig1
            dpos = med_nan(sp_t) - med_nan(sp_id);
            dneg = med_nan(sn_t) - med_nan(sn_id);
            dsep = dpos - dneg;
            acc.off.delta_pos(end+1) = dpos;
            acc.off.delta_neg(end+1) = dneg;
            acc.off.delta_sep(end+1) = dsep;

            % --- Δsep (Task & Fix) → Fig2 inset
            dsep_task = (med_nan(sp_t)-med_nan(sp_id)) - (med_nan(sn_t)-med_nan(sn_id));
            dsep_fix  = (med_nan(sp_f)-med_nan(sp_id)) - (med_nan(sn_f)-med_nan(sn_id));
            acc.off.dsep_task(end+1) = dsep_task;
            acc.off.dsep_fix(end+1)  = dsep_fix;

            % --- ROC traces → Fig2
            [fpr_t, tpr_t] = roc_on_grid(sp_t, sn_t, FPR_GRID);
            [fpr_f, tpr_f] = roc_on_grid(sp_f, sn_f, FPR_GRID);
            T_roc = add_rows(T_roc, sid, r, 'OFFSET','TASK', fpr_t, tpr_t);
            T_roc = add_rows(T_roc, sid, r, 'OFFSET','FIX',  fpr_f, tpr_f);

            % --- AUC per subject-run
            T_auc = add_auc_row(T_auc, sid, r, 'OFFSET','TASK', auc_from_scores(sp_t, sn_t));
            T_auc = add_auc_row(T_auc, sid, r, 'OFFSET','FIX',  auc_from_scores(sp_f, sn_f));
        end
    end

    % --- pack subject-level Δ means for Fig1
    T_bias = [T_bias;
        mk_bias_row(sid,'ONSET','Delta_pos',  mean(acc.on.delta_pos, 'omitnan'));
        mk_bias_row(sid,'ONSET','Delta_neg',  mean(acc.on.delta_neg, 'omitnan'));
        mk_bias_row(sid,'ONSET','Delta_sep',  mean(acc.on.delta_sep, 'omitnan'));
        mk_bias_row(sid,'OFFSET','Delta_pos', mean(acc.off.delta_pos,'omitnan'));
        mk_bias_row(sid,'OFFSET','Delta_neg', mean(acc.off.delta_neg,'omitnan'));
        mk_bias_row(sid,'OFFSET','Delta_sep', mean(acc.off.delta_sep,'omitnan')) ];

    % --- Δsep inset (per-subject means across runs) for TASK & FIX
    T_dsep = [T_dsep;
        table(sid, categorical({'ONSET'}) , categorical({'TASK'}),  mean(acc.on.dsep_task,'omitnan'), ...
              'VariableNames',{'subj','phase','scheme','dsep'});
        table(sid, categorical({'ONSET'}) , categorical({'FIX' }),  mean(acc.on.dsep_fix, 'omitnan'), ...
              'VariableNames',{'subj','phase','scheme','dsep'});
        table(sid, categorical({'OFFSET'}), categorical({'TASK'}),  mean(acc.off.dsep_task,'omitnan'), ...
              'VariableNames',{'subj','phase','scheme','dsep'});
        table(sid, categorical({'OFFSET'}), categorical({'FIX' }),  mean(acc.off.dsep_fix, 'omitnan'), ...
              'VariableNames',{'subj','phase','scheme','dsep'}) ];
end

% ----- write CSVs
writetable(T_bias, fullfile(outdir,'fig1_bias_task_vs_identity.csv'));
writetable(T_roc,  fullfile(outdir,'fig2_roc_traces.csv'));
writetable(T_auc,  fullfile(outdir,'fig2_auc_subject_run.csv'));
writetable(T_dsep, fullfile(outdir,'fig2_dsep_inset.csv'));

fprintf('Wrote:\n  %s\n  %s\n  %s\n  %s\n', ...
    fullfile(outdir,'fig1_bias_task_vs_identity.csv'), ...
    fullfile(outdir,'fig2_roc_traces.csv'), ...
    fullfile(outdir,'fig2_auc_subject_run.csv'), ...
    fullfile(outdir,'fig2_dsep_inset.csv'));
end

%% ======================= Local helpers =======================

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
% Log-Euclidean mean covariance across all windows in a {runs}[T x C x N] cell array
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

function row = mk_bias_row(sid, phase, metric, val)
row = table( sid, categorical({phase}), categorical({metric}), val, ...
             'VariableNames', {'subj','phase','metric','value'} );
end

function T = add_rows(T, sid, run, phase, scheme, fpr, tpr)
n = numel(fpr);
TT = table( repmat(sid,n,1), repmat(run,n,1), ...
            repmat(categorical({phase}),n,1), repmat(categorical({scheme}),n,1), ...
            fpr(:), tpr(:), ...
            'VariableNames', {'subj','run','phase','scheme','fpr','tpr'} );
if isempty(T), T = TT; else, T = [T; TT]; end
end

function T = add_auc_row(T, sid, run, phase, scheme, auc)
TT = table(sid, run, categorical({phase}), categorical({scheme}), auc, ...
           'VariableNames', {'subj','run','phase','scheme','auc'});
if isempty(T), T = TT; else, T = [T; TT]; end
end

function [fpr_grid, tpr_grid] = roc_on_grid(sp, sn, fpr_grid)
% Rank-based ROC on a uniform FPR grid (lower s = more positive)
sp = sp(:); sn = sn(:);
sp = sp(isfinite(sp)); sn = sn(isfinite(sn));
if isempty(sp) || isempty(sn)
    tpr_grid = nan(size(fpr_grid));
    return;
end
scores = [sp; sn];
labels = [ones(numel(sp),1); zeros(numel(sn),1)];
[sc,ord] = sort(scores,'ascend'); labs = labels(ord);
tp = cumsum(labs==1); fp = cumsum(labs==0);
P  = sum(labs==1);    N  = sum(labs==0);
TPR = tp / max(P,1);
FPR = fp / max(N,1);
% unique FPR for interpolation
[Funiq, iu] = unique(FPR);
Tuniq = TPR(iu);
if numel(Funiq)==1
    tpr_grid = repmat(Tuniq, size(fpr_grid));
else
    tpr_grid = interp1(Funiq, Tuniq, fpr_grid, 'linear','extrap');
    tpr_grid = max(0,min(1,tpr_grid));
end
end

function auc = auc_from_scores(sp, sn)
% Mann–Whitney U formulation with tie handling
sp = sp(isfinite(sp)); sn = sn(isfinite(sn));
np = numel(sp); nn = numel(sn);
if np==0 || nn==0, auc = NaN; return; end
all_scores = [sp(:); sn(:)];
rk = tiedrank(all_scores);      % ascending ranks (lower s = more positive)
rp = rk(1:np);
auc = (sum(rp) - np*(np+1)/2) / (np*nn);
end

% ======================= Robust fixation builder =======================
function X_fix = get_fix_for_run(Son, Sof, r)
% Find neutral/fixation windows for run r; fallbacks to offline REST.
X_fix = [];
candidates = {'on_fix_cell','on_fixation_cell','on_preCue_cell','on_pre_fix_cell','on_rest_cell','on_countdown_cell'};
for k = 1:numel(candidates)
    if isfield(Son, candidates{k})
        X_fix = pick3d(Son.(candidates{k}), r);
        if ~isempty(X_fix), break; end
    end
end
if isempty(X_fix) && isfield(Sof,'data_rest_cell')
    X_fix = pick3d(Sof.data_rest_cell, r);  % per-run offline REST if available
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

function Cstack = covs_from_cellarr(cellArr)
Cstack = [];
for r = 1:size(cellArr,1)
    X = cellArr{r,1};
    if isempty(X), continue; end
    Cstack = cat_covs(Cstack, covs_from_3d(X));
end
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

function Rfix = build_R_fixation_robust(Son, Sof, r, Rrest_global, opts)
% Robust, label-agnostic fixation reference for run r.
% - neutral windows
% - outlier trim by SPD distance (keep_frac)
% - shrinkage to global REST in log-Euclidean space (alpha)
% - diagonal loading (lambda)
if ~isfield(opts,'alpha'),     opts.alpha     = 0.30; end
if ~isfield(opts,'lambda'),    opts.lambda    = 1e-3; end
if ~isfield(opts,'keep_frac'), opts.keep_frac = 0.90; end

X_fix = get_fix_for_run(Son, Sof, r);
Cfix  = covs_from_3d(X_fix);
if isempty(Cfix)
    Rfix = ensure_spd(Rrest_global);
    return;
end

% Outlier trimming in SPD space
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
