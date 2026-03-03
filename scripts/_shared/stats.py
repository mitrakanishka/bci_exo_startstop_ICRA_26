"""Small statistics helpers with no SciPy dependency."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


@dataclass
class WilcoxonResult:
    n: int
    statistic: float
    pvalue: float


def _rankdata_average(x: np.ndarray) -> np.ndarray:
    order = np.argsort(x)
    ranks = np.empty(len(x), dtype=float)
    i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and x[order[j + 1]] == x[order[i]]:
            j += 1
        avg_rank = 0.5 * (i + j) + 1.0
        ranks[order[i : j + 1]] = avg_rank
        i = j + 1
    return ranks


def _normal_cdf(z: float) -> float:
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def _exact_p_from_signed_ranks(ranks: np.ndarray, w_obs: float) -> float:
    """
    Exact two-sided p-value for Wilcoxon statistic W=min(W+,W-) via subset-sum DP.

    Ranks may include 0.5 increments from average-rank ties, so we scale by 2.
    """
    scaled = np.rint(ranks * 2.0).astype(int)
    total = int(np.sum(scaled))
    obs = int(round(w_obs * 2.0))

    counts = np.zeros(total + 1, dtype=np.int64)
    counts[0] = 1
    for r in scaled:
        counts[r:] += counts[:-r]

    n = len(scaled)
    denom = float(1 << n)
    p = 0.0
    for s, c in enumerate(counts):
        if c == 0:
            continue
        if min(s, total - s) <= obs:
            p += c / denom
    return max(0.0, min(1.0, float(p)))


def paired_wilcoxon(x: np.ndarray, y: np.ndarray) -> WilcoxonResult:
    """Two-sided Wilcoxon signed-rank p-value (exact for n<=25, else asymptotic)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.shape != y.shape:
        raise ValueError("x and y must have the same shape")

    d = y - x
    d = d[np.isfinite(d)]
    d = d[d != 0]
    n = d.size
    if n < 2:
        return WilcoxonResult(n=n, statistic=float("nan"), pvalue=float("nan"))

    abs_d = np.abs(d)
    ranks = _rankdata_average(abs_d)
    w_plus = float(np.sum(ranks[d > 0]))
    w_minus = float(np.sum(ranks[d < 0]))
    w = min(w_plus, w_minus)

    # Exact p-values are cheap and stable for our paper sample sizes (n<=8).
    if n <= 25:
        return WilcoxonResult(n=n, statistic=w, pvalue=_exact_p_from_signed_ranks(ranks, w))

    mean_w = n * (n + 1) / 4.0
    var_w = n * (n + 1) * (2 * n + 1) / 24.0
    tie_counts = np.unique(abs_d, return_counts=True)[1]
    if tie_counts.size:
        tie_term = np.sum(tie_counts * (tie_counts + 1) * (2 * tie_counts + 1))
        var_w -= tie_term / 48.0
    if var_w <= 0:
        return WilcoxonResult(n=n, statistic=w, pvalue=float("nan"))

    z = (w - mean_w) / math.sqrt(var_w)
    p = 2.0 * (1.0 - _normal_cdf(abs(z)))
    p = max(0.0, min(1.0, p))
    return WilcoxonResult(n=n, statistic=w, pvalue=p)
