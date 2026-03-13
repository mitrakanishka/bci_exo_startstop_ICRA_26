#!/usr/bin/env python3
"""Paper Figure 5: Task-based recentering bias shift vs identity."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parent
sys.path.insert(0, str(THIS_DIR))

from _shared.recentering_analysis import compute_recentering_tables
from _shared.stats import paired_wilcoxon


PHASE_ORDER = ["ONSET", "OFFSET"]
METRIC_ORDER = ["Delta_pos", "Delta_neg", "Delta_sep"]
FONT_FAMILY = "DejaVu Sans"

COLORS = {
    "ONSET": {"Delta_pos": "#2E7D32", "Delta_neg": "#1E88E5", "Delta_sep": "#C9C9C9"},
    "OFFSET": {"Delta_pos": "#CD5C5C", "Delta_neg": "#7E57C2", "Delta_sep": "#C9C9C9"},
}


def _bootstrap_ci(values: np.ndarray, n_boot: int = 10000, seed: int = 0) -> tuple[float, float]:
    values = values[np.isfinite(values)]
    if values.size == 0:
        return (np.nan, np.nan)
    rng = np.random.RandomState(seed)
    means = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        sample = values[rng.randint(0, values.size, size=values.size)]
        means[i] = np.mean(sample)
    return (float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5)))


def make_figure(df: pd.DataFrame, out_png: Path, out_pdf: Path, out_csv: Path) -> None:
    df = df.copy()
    df["phase"] = pd.Categorical(df["phase"], PHASE_ORDER)
    df["metric"] = pd.Categorical(df["metric"], METRIC_ORDER)

    rows = []
    for ph in PHASE_ORDER:
        for met in METRIC_ORDER:
            vals = df[(df["phase"] == ph) & (df["metric"] == met)]["value"].to_numpy(dtype=float)
            mean = float(np.nanmean(vals))
            lo, hi = _bootstrap_ci(vals, seed=(7 if ph == "ONSET" else 17) + METRIC_ORDER.index(met))
            rows.append({"phase": ph, "metric": met, "mean": mean, "ci_low": lo, "ci_high": hi, "n": int(np.sum(np.isfinite(vals)))})

    summary = pd.DataFrame(rows)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out_csv, index=False)

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 5.4), dpi=300, sharey=False)
    fig.suptitle(
        r"Task Recentering vs Identity ($\Delta$ Median Margin)",
        fontsize=16,
        fontweight="bold",
        fontfamily=FONT_FAMILY,
        y=0.95,
    )

    labels = {
        "ONSET": [r"$\Delta_{pos}$" + "\n(Start-MI)", r"$\Delta_{neg}$" + "\n(REST)", r"$\Delta_{sep}$"],
        "OFFSET": [r"$\Delta_{pos}$" + "\n(Stop-MI)", r"$\Delta_{neg}$" + "\n(Maintain-MI)", r"$\Delta_{sep}$"],
    }

    for i, ph in enumerate(PHASE_ORDER):
        ax = axes[i]
        vals_phase = summary[summary["phase"] == ph].set_index("metric").loc[METRIC_ORDER]
        x = np.arange(len(METRIC_ORDER))
        means = vals_phase["mean"].to_numpy(dtype=float)
        lo = vals_phase["ci_low"].to_numpy(dtype=float)
        hi = vals_phase["ci_high"].to_numpy(dtype=float)
        yerr = np.vstack((means - lo, hi - means))

        colors = [COLORS[ph][m] for m in METRIC_ORDER]
        ax.bar(x, means, yerr=yerr, color=colors, edgecolor="black", linewidth=0.5, capsize=4)

        ax.axhline(0.0, color="0.55", lw=1, ls="--")
        ax.set_xticks(x)
        ax.set_xticklabels(labels[ph], fontsize=10, fontweight="bold", fontfamily=FONT_FAMILY)
        ax.set_title(ph.title(), fontsize=16, fontweight="bold", fontfamily=FONT_FAMILY)
        ax.set_xlabel("")
        if i == 0:
            ax.set_ylabel(
                r"$\Delta$ vs $S_I$ (median margin)",
                fontsize=18,
                fontweight="bold",
                fontfamily=FONT_FAMILY,
            )

        dsep = df[(df["phase"] == ph) & (df["metric"] == "Delta_sep")]["value"].to_numpy(dtype=float)
        res = paired_wilcoxon(np.zeros_like(dsep), dsep)
        ptxt = (
            rf"$\Delta_{{sep}}$ vs $S_I$: p={res.pvalue:.4f}"
            if np.isfinite(res.pvalue)
            else rf"$\Delta_{{sep}}$ vs $S_I$: n={res.n}"
        )
        ax.text(
            0.5,
            0.04,
            ptxt,
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=15,
            fontfamily=FONT_FAMILY,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.85"),
        )

        ax.grid(axis="y", linestyle=":", linewidth=0.8, alpha=0.45)
        ax.spines[["top", "right"]].set_visible(False)
        for side in ["left", "bottom"]:
            ax.spines[side].set_visible(True)
            ax.spines[side].set_color("black")
            ax.spines[side].set_linewidth(0.8)
        ax.tick_params(axis="both", labelsize=14)
        for tick in ax.get_xticklabels():
            tick.set_fontweight("bold")
            tick.set_fontfamily(FONT_FAMILY)
        for tick in ax.get_yticklabels():
            tick.set_fontfamily(FONT_FAMILY)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_png)
    fig.savefig(out_pdf)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--of-dir",
        type=Path,
        default=REPO_ROOT / "BCI_Harmony_ExperimentalData" / "epoched_data_of",
        help="Path to offline epoched-data directory.",
    )
    parser.add_argument(
        "--on-dir",
        type=Path,
        default=REPO_ROOT / "BCI_Harmony_ExperimentalData" / "epoched_data_on",
        help="Path to online epoched-data directory.",
    )
    parser.add_argument(
        "--subject-csv",
        type=Path,
        default=REPO_ROOT / "data" / "processed" / "fig5_bias_task_vs_identity.csv",
        help="Path to save per-subject bias metrics.",
    )
    parser.add_argument(
        "--max-runs",
        type=int,
        default=8,
        help="Maximum number of online runs per subject to analyze.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=REPO_ROOT / "figures" / "paper",
        help="Output figure directory.",
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=REPO_ROOT / "data" / "processed" / "fig5_summary_stats.csv",
        help="Path to save mean/CI summary table.",
    )
    args = parser.parse_args()

    df, _ = compute_recentering_tables(args.of_dir, args.on_dir, max_runs=args.max_runs)
    args.subject_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.subject_csv, index=False)

    make_figure(
        df,
        out_png=args.outdir / "fig5_bias_shift_vs_identity.png",
        out_pdf=args.outdir / "fig5_bias_shift_vs_identity.pdf",
        out_csv=args.summary_csv,
    )
    print("Saved Figure 5 outputs to", args.outdir)
    print("Saved subject-level bias metrics to", args.subject_csv)


if __name__ == "__main__":
    main()
