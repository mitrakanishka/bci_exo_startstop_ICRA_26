function res = paired_wilcoxon_exact(x, y)
% Two-sided Wilcoxon signed-rank test with exact p-values for small n.

x = x(:);
y = y(:);
assert(numel(x) == numel(y), 'x and y must have the same number of elements.');

d = y - x;
d = d(isfinite(d));
d = d(d ~= 0);
res = struct('n', numel(d), 'statistic', NaN, 'pvalue', NaN);
if numel(d) < 2
    return;
end

abs_d = abs(d);
ranks = rankdata_average(abs_d);
w_plus = sum(ranks(d > 0));
w_minus = sum(ranks(d < 0));
w = min(w_plus, w_minus);
res.statistic = w;

if numel(d) <= 25
    res.pvalue = exact_p_from_signed_ranks(ranks, w);
    return;
end

n = numel(d);
mean_w = n * (n + 1) / 4;
var_w = n * (n + 1) * (2 * n + 1) / 24;
[uvals, ~, ic] = unique(abs_d);
counts = accumarray(ic, 1, [numel(uvals), 1]);
if ~isempty(counts)
    var_w = var_w - sum(counts .* (counts + 1) .* (2 * counts + 1)) / 48;
end
if var_w <= 0
    return;
end
z = (w - mean_w) / sqrt(var_w);
res.pvalue = max(0, min(1, 2 * (1 - normcdf_local(abs(z)))));
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

function p = exact_p_from_signed_ranks(ranks, w_obs)
scaled = round(ranks * 2);
total = sum(scaled);
obs = round(w_obs * 2);
counts = zeros(total + 1, 1);
counts(1) = 1;
for i = 1:numel(scaled)
    r = scaled(i);
    counts((r + 1):end) = counts((r + 1):end) + counts(1:(end - r));
end
n = numel(scaled);
denom = 2 ^ n;
p = 0;
for s = 0:total
    c = counts(s + 1);
    if c == 0
        continue;
    end
    if min(s, total - s) <= obs
        p = p + c / denom;
    end
end
p = max(0, min(1, p));
end

function p = normcdf_local(z)
p = 0.5 * (1 + erf(z / sqrt(2)));
end
