function GROUP_BIAS_AUC_ALL()
% GROUP_BIAS_AUC_ALL  (8-run version)
% Group-level plot of AUC (Biased vs No-recenter) with shaded 95% CI bands.
% - ONSET:  pos = on_startMI, neg = online fixation
% - OFFSET: pos = on_stopMI,  neg = online fixation
% "No-recenter": identity reference.
% "Biased": run-level reference = mean covariance of that run's positive windows.
%
% Hard limit to 8 runs across printing, figures, averaging, and stats.

clc;

repo_root = resolve_repo_root(fileparts(mfilename('fullpath')));
root     = fullfile(repo_root,'BCI_course_EXP');
of_dir   = fullfile(root,'epoched_data_of');
on_dir   = fullfile(root,'epoched_data_on');
assert(exist(of_dir,'dir')==7 && exist(on_dir,'dir')==7, ...
    'Could not locate epoched data folders under %s', root);

% ---- Discover subjects present in BOTH folders ----
of_files = dir(fullfile(of_dir, 'subj_*_epoched_data_of.mat'));
on_files = dir(fullfile(on_dir, 'subj_*_epoched_data_on.mat'));
sub_of   = parse_sub_ids(of_files);
sub_on   = parse_sub_ids(on_files);
subs     = intersect(sub_of, sub_on);

if isempty(subs)
    error('No overlapping subjects found in epoched_data_of/ and epoched_data_on/.');
end
fprintf('Found %d subjects: %s\n', numel(subs), sprintf('%03d ', subs));

% ---- Store per-subject per-run AUC vectors (fixed length = 8) ----
NRUNS = 8;   % <<--- hard cap to 8
all_on_norc  = cell(numel(subs),1);
all_on_bias  = cell(numel(subs),1);
all_off_norc = cell(numel(subs),1);
all_off_bias = cell(numel(subs),1);

for si = 1:numel(subs)
    sid = subs(si);

    Sof = load(fullfile(of_dir,  sprintf('subj_%03d_epoched_data_of.mat', sid)));
    Son = load(fullfile(on_dir,  sprintf('subj_%03d_epoched_data_on.mat', sid)));

    % Online field names (robust to earlier variants)
    onStart = get_first_existing(Son, {'on_startMI_cell','on_bMI_cell'});
    onStop  = get_first_existing(Son, {'on_stopMI_cell','on_eMI_cell'});
    onFix   = get_first_existing(Son, {'on_fix3s_cell','on_fix_cell'});
    assert(~isempty(onStart) && ~isempty(onStop) && ~isempty(onFix), ...
        'Missing online epoch fields for subject %03d.', sid);

    % Offline prototypes (No-recenter training)
    % ONSET: Rest(0) vs Begin-MI(1)
    P0_on = proto_from_cell(Sof.data_rest_cell);
    P1_on = proto_from_cell(Sof.data_bMI_cell);

    % OFFSET: Doing-MI(0) vs End-MI(1)
    P0_off = proto_from_cell(Sof.data_dMI_cell);
    P1_off = proto_from_cell(Sof.data_eMI_cell);

    % Available runs, then clip to NRUNS
    nRunsAvail = max([size(onStart,1), size(onStop,1), size(onFix,1)]);
    nRuns = min(nRunsAvail, NRUNS);  % <<--- clip to 8

    on_auc_norc  = nan(1, NRUNS);  % fixed length arrays
    on_auc_bias  = nan(1, NRUNS);
    off_auc_norc = nan(1, NRUNS);
    off_auc_bias = nan(1, NRUNS);

    for r = 1:nRuns
        Xfix  = pick3d(onFix,   r);   % [T x C x Nf]
        Xon   = pick3d(onStart, r);   % [T x C x Np_on]
        Xoff  = pick3d(onStop,  r);   % [T x C x Np_off]

        % ---- ONSET: pos=startMI, neg=fixation ----
        if ~isempty(Xon) && ~isempty(Xfix)
            % No-recenter
            [s_pos, s_neg] = score_windows_vs_protos(Xon, Xfix, P0_on, P1_on, eye(size(P0_on)));
            on_auc_norc(r) = auc_from_scores([ones(numel(s_pos),1); zeros(numel(s_neg),1)], ...
                                             [s_pos; s_neg]);
            % Biased (run reference = mean cov of Xon)
            R_on = ref_from_3d(Xon);
            [s_posb, s_negb] = score_windows_vs_protos(Xon, Xfix, P0_on, P1_on, R_on);
            on_auc_bias(r) = auc_from_scores([ones(numel(s_posb),1); zeros(numel(s_negb),1)], ...
                                             [s_posb; s_negb]);
        end

        % ---- OFFSET: pos=stopMI, neg=fixation ----
        if ~isempty(Xoff) && ~isempty(Xfix)
            % No-recenter
            [s_pos, s_neg] = score_windows_vs_protos(Xoff, Xfix, P0_off, P1_off, eye(size(P0_off)));
            off_auc_norc(r) = auc_from_scores([ones(numel(s_pos),1); zeros(numel(s_neg),1)], ...
                                              [s_pos; s_neg]);
            % Biased (run reference = mean cov of Xoff)
            R_off = ref_from_3d(Xoff);
            [s_posb, s_negb] = score_windows_vs_protos(Xoff, Xfix, P0_off, P1_off, R_off);
            off_auc_bias(r) = auc_from_scores([ones(numel(s_posb),1); zeros(numel(s_negb),1)], ...
                                              [s_posb; s_negb]);
        end
    end

    all_on_norc{si}  = on_auc_norc;
    all_on_bias{si}  = on_auc_bias;
    all_off_norc{si} = off_auc_norc;
    all_off_bias{si} = off_auc_bias;
end

% ---- Group per-run means + SD + 95% bootstrap CI across subjects ----
B = 10000; alpha = 0.95;

grp_on_norc  = perrun_boot_stats(all_on_norc,  B, alpha, NRUNS); % <<--- enforce 8 runs
grp_on_bias  = perrun_boot_stats(all_on_bias,  B, alpha, NRUNS);
grp_off_norc = perrun_boot_stats(all_off_norc, B, alpha, NRUNS);
grp_off_bias = perrun_boot_stats(all_off_bias, B, alpha, NRUNS);

% ---- Print per-run group means ± SD (8 runs only) ----
fprintf('\n=== Per-run group AUC: mean ± SD (No-recenter / Biased) ===\n');

fprintf('Onset  (runs 1..%d):\n', NRUNS);
print_run_means('  No-recenter', grp_on_norc.mean, grp_on_norc.std);
print_run_means('  Biased     ', grp_on_bias.mean, grp_on_bias.std);

fprintf('Offset (runs 1..%d):\n', NRUNS);
print_run_means('  No-recenter', grp_off_norc.mean, grp_off_norc.std);
print_run_means('  Biased     ', grp_off_bias.mean, grp_off_bias.std);

% ---- Figure: group mean curves + shaded CI (8 runs) ----
figure('Color','w','Name','Group AUC per run — No-recenter vs Biased (8 runs)');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

colN = [0.50 0.50 0.50];   % grey (No-recenter)
colB = [0.22 0.49 0.72];   % blue (Biased)
lw   = 2.0; fa = 0.18;

% ========== ONSET ==========
nexttile; hold on; title('Onset'); xlabel('Run #'); ylabel('AUC');
hChance = yline(0.5,':','Chance','Color',[0 0 0 0.5]); set(hChance,'HandleVisibility','off');

hN_CI = plot_ribbon(grp_on_norc, colN, fa);
hN    = plot(grp_on_norc.x, grp_on_norc.mean, 'Color', colN, 'LineWidth', lw);

hB_CI = plot_ribbon(grp_on_bias, colB, fa);
hB    = plot(grp_on_bias.x,  grp_on_bias.mean,  'Color', colB, 'LineWidth', lw);

legend([hN_CI, hN, hB_CI, hB], ...
       {'No-recenter (CI)', 'No-recenter mean', 'Biased (CI)', 'Biased mean'}, ...
       'Location','southoutside','Orientation','horizontal','Box','off','FontSize',8);

xlim([1 NRUNS]); ylim([0 1]); box off;

% ========== OFFSET ==========
nexttile; hold on; title('Offset'); xlabel('Run #'); ylabel('AUC');
hChance = yline(0.5,':','Chance','Color',[0 0 0 0.5]); set(hChance,'HandleVisibility','off');

hN_CI2 = plot_ribbon(grp_off_norc, colN, fa);
hN2    = plot(grp_off_norc.x, grp_off_norc.mean, 'Color', colN, 'LineWidth', lw);

hB_CI2 = plot_ribbon(grp_off_bias, colB, fa);
hB2    = plot(grp_off_bias.x,  grp_off_bias.mean,  'Color', colB, 'LineWidth', lw);

xlim([1 NRUNS]); ylim([0 1]); box off;

% ---- Group summary & stats (subject-level means across first 8 runs) ----
% Per subject, mean AUC across first NRUNS runs for each condition
sub_on_norc  = cellfun(@(v) nanmean(v(1:NRUNS)), all_on_norc);
sub_on_bias  = cellfun(@(v) nanmean(v(1:NRUNS)), all_on_bias);
sub_off_norc = cellfun(@(v) nanmean(v(1:NRUNS)), all_off_norc);
sub_off_bias = cellfun(@(v) nanmean(v(1:NRUNS)), all_off_bias);

% Print group means ± SD
fprintf('\n=== Group means ± SD (AUC across runs 1..%d) ===\n', NRUNS);
fprintf('Onset  No-recenter: %0.3f ± %0.3f\n', mean(sub_on_norc,'omitnan'),  std(sub_on_norc,'omitnan'));
fprintf('Onset  Biased     : %0.3f ± %0.3f\n', mean(sub_on_bias,'omitnan'),  std(sub_on_bias,'omitnan'));
fprintf('Offset No-recenter: %0.3f ± %0.3f\n', mean(sub_off_norc,'omitnan'), std(sub_off_norc,'omitnan'));
fprintf('Offset Biased     : %0.3f ± %0.3f\n', mean(sub_off_bias,'omitnan'), std(sub_off_bias,'omitnan'));

% ΔAUC per subject (biased − no-recenter), bootstrap CI and Wilcoxon
sub_d_on  = sub_on_bias  - sub_on_norc;
sub_d_off = sub_off_bias - sub_off_norc;

B = 10000; alpha = 0.95;
[mu_on,  ci_on]  = mean_ci_bootstrap(sub_d_on,  B, alpha);
[mu_off, ci_off] = mean_ci_bootstrap(sub_d_off, B, alpha);

[p_on,  ~, stats_on]  = signrank(sub_d_on,  0, 'method','approximate');
[p_off, ~, stats_off] = signrank(sub_d_off, 0, 'method','approximate');
z_on  = field_or(stats_on,  'zval', NaN);
z_off = field_or(stats_off, 'zval', NaN);

fprintf('\n=== Group ΔAUC (Biased − No-recenter) over runs 1..%d ===\n', NRUNS);
fprintf('Onset : mean = %+0.3f, 95%% CI [%+0.3f, %+0.3f], Wilcoxon z=%.3f, p=%.4f, n=%d\n', ...
    mu_on, ci_on(1), ci_on(2), z_on,  p_on,  numel(sub_d_on));
fprintf('Offset: mean = %+0.3f, 95%% CI [%+0.3f, %+0.3f], Wilcoxon z=%.3f, p=%.4f, n=%d\n', ...
    mu_off, ci_off(1), ci_off(2), z_off, p_off, numel(sub_d_off));

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
if isempty(cellArr), return; end
if r <= size(cellArr,1) && ~isempty(cellArr{r,1})
    X = cellArr{r,1};
end
end

function P = proto_from_cell(cellArr)
% Log-Euclidean mean covariance over all windows from a cell array {run}[T x C x N]
Cstack = [];
for r = 1:size(cellArr,1)
    X = cellArr{r,1};
    if isempty(X), continue; end
    C = covs_from_3d(X);
    Cstack = cat_covs(Cstack, C);
end
if isempty(Cstack), error('Empty training pool when building prototype.'); end
P = logeu_mean(Cstack);
end

function C = covs_from_3d(X3)
% X3: [T x C x N] → covariances (C x C x N)
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
% Run-level reference: log-Euclidean mean of covariances from X3
C = covs_from_3d(X3);
if isempty(C)
    R = eye(size(X3,2));
else
    R = logeu_mean(C);
end
end

function [s_pos, s_neg] = score_windows_vs_protos(Xpos, Xneg, P0, P1, R)
% Score = (d0 - d1), using affine-invariant distance after whitening by R
Cpos = covs_from_3d(Xpos); Cneg = covs_from_3d(Xneg);
s_pos = score_covs(Cpos, P0, P1, R);
s_neg = score_covs(Cneg, P0, P1, R);
end

function s = score_covs(Cstack, P0, P1, R)
if isempty(Cstack), s = []; return; end
P0 = ensure_spd(P0); P1 = ensure_spd(P1); R = ensure_spd(R);
% R^{-1/2}
[V,D] = eig((R+R')/2); lam = max(real(diag(D)), 1e-12);
Rminv = V * diag(1./sqrt(lam)) * V';
n = size(Cstack,3);
s = zeros(n,1);
for i = 1:n
    C = ensure_spd(Cstack(:,:,i));
    Cw = (Rminv * C * Rminv'); Cw = (Cw + Cw')/2;
    d0 = spd_dist(Cw, P0);
    d1 = spd_dist(Cw, P1);
    s(i) = (d0 - d1);  % larger → closer to class 1
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

function C = cat_covs(A,B)
if isempty(A), C = B; elseif isempty(B), C = A; else, C = cat(3,A,B); end
end

function P = logeu_mean(Cstack)
% Log-Euclidean mean of SPD stack
n = size(Cstack,3); C = size(Cstack,1);
Lsum = zeros(C);
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

function A = ensure_spd(A)
A = (A + A')/2;
n = size(A,1);
epsw = 1e-9 * trace(A)/max(n,1);
if ~isfinite(epsw) || epsw<=0, epsw = 1e-9; end
A = A + epsw*eye(n);
end

function auc = auc_from_scores(y, s)
% y in {0,1}, s: larger means class 1 more likely
y = y(:); s = s(:);
mask = isfinite(y) & isfinite(s);
y = y(mask); s = s(mask);
if numel(unique(y)) < 2, auc = NaN; return; end
% Mann-Whitney U / rank-AUC
r = tiedrank(s);
pos = (y==1); neg = (y==0);
n1 = sum(pos); n0 = sum(neg);
sumRanksPos = sum(r(pos));
auc = (sumRanksPos - n1*(n1+1)/2) / (n1*n0);
end

function [mu, ci] = mean_ci_bootstrap(x, B, alpha)
x = x(:); x = x(isfinite(x));
n = numel(x);
if n==0, mu = NaN; ci = [NaN NaN]; return; end
if nargin < 2, B = 10000; end
if nargin < 3, alpha = 0.95; end
mu = mean(x);
boots = zeros(B,1);
for b = 1:B
    idx = randi(n, n, 1);
    boots(b) = mean(x(idx));
end
lo = (1-alpha)/2; hi = 1-lo;
ci = quantile(boots, [lo hi]);
end

function out = perrun_boot_stats(cellVecs, B, alpha, NRUNS)
% cellVecs: {S x 1}, each row a 1xR vector (NaNs allowed); restrict to 1..NRUNS
S = numel(cellVecs);
X = nan(S, NRUNS);
for s = 1:S
    v = cellVecs{s};
    v = v(:).';                % row
    v = v(1:min(numel(v), NRUNS));     % clip to NRUNS
    X(s,1:numel(v)) = v;
end
out.x    = 1:NRUNS;
out.mean = nanmean(X,1);
out.std  = nanstd(X,0,1);
out.n    = sum(isfinite(X),1);
[out.lo, out.hi] = deal(nan(1,NRUNS));
for r = 1:NRUNS
    xr = X(:,r); xr = xr(isfinite(xr));
    if numel(xr) >= 2
        [m, ci] = mean_ci_bootstrap(xr, B, alpha);
        out.mean(r) = m;
        out.lo(r)   = ci(1);
        out.hi(r)   = ci(2);
    else
        out.lo(r) = out.mean(r);
        out.hi(r) = out.mean(r);
    end
end
end

function h = plot_ribbon(st, color, faceAlpha)
% st has fields: x, mean, lo, hi
x  = st.x;
lo = st.lo;
hi = st.hi;
h = fill([x, fliplr(x)], [lo, fliplr(hi)], color, ...
         'FaceAlpha', faceAlpha, 'EdgeColor', 'none', 'HandleVisibility', 'on');
uistack(h,'bottom');  % keep the ribbon behind the mean line
end

function v = field_or(S, name, defaultVal)
if isstruct(S) && isfield(S, name)
    v = S.(name);
else
    v = defaultVal;
end
end

function print_run_means(label, m, s)
% Pretty print as one row: label:  m1±s1  m2±s2  ... (assumes 8 runs)
fprintf('%s : ', label);
for i = 1:numel(m)
    if i == numel(m)
        fprintf('%0.3f±%0.3f\n', m(i), s(i));
    else
        fprintf('%0.3f±%0.3f  ', m(i), s(i));
    end
end
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
