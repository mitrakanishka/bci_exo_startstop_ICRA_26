#!/usr/bin/env python3
"""Paper Figure 6: AUC by run for TASK vs fixation-based recentering."""

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
SCHEME_ORDER = ["TASK", "FIX"]
SCHEME_COLORS = {"TASK": "#E6862A", "FIX": "#2979FF"}
FONT_FAMILY = "DejaVu Sans"
LINE_WIDTH = 2.6
RUNS = np.arange(1, 9)


def run_to_session_label(r: int) -> str:
    return f"S2,R{r}" if r <= 4 else f"S3,R{r - 4}"


def _csv_summary(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, float]]:
    df = df.copy()
    df["phase"] = pd.Categorical(df["phase"], PHASE_ORDER)
    df["scheme"] = pd.Categorical(df["scheme"], SCHEME_ORDER)

    rows = []
    pvals: dict[str, float] = {}
    for phase in PHASE_ORDER:
        for scheme in SCHEME_ORDER:
            for run in RUNS:
                vals = df[(df["phase"] == phase) & (df["scheme"] == scheme) & (df["run"] == run)][
                    "auc"
                ].to_numpy(dtype=float)
                vals = vals[np.isfinite(vals)]
                if vals.size == 0:
                    rows.append(
                        {
                            "phase": phase,
                            "scheme": scheme,
                            "run": int(run),
                            "mean": np.nan,
                            "sd": np.nan,
                            "ci_low": np.nan,
                            "ci_high": np.nan,
                            "n": 0,
                        }
                    )
                    continue
                m = float(np.mean(vals))
                sd = float(np.std(vals, ddof=0))
                half = 1.96 * sd / np.sqrt(vals.size)
                rows.append(
                    {
                        "phase": phase,
                        "scheme": scheme,
                        "run": int(run),
                        "mean": m,
                        "sd": sd,
                        "ci_low": m - half,
                        "ci_high": m + half,
                        "n": int(vals.size),
                    }
                )

        per_subj = (
            df[df["phase"] == phase]
            .groupby(["subj", "scheme"])["auc"]
            .mean()
            .unstack("scheme")
            .reindex(columns=SCHEME_ORDER)
            .dropna()
        )
        if per_subj.empty:
            pvals[phase] = float("nan")
        else:
            res = paired_wilcoxon(per_subj["TASK"].to_numpy(), per_subj["FIX"].to_numpy())
            pvals[phase] = float(res.pvalue)

    return pd.DataFrame(rows), pvals


def make_figure(summary: pd.DataFrame, pvals: dict[str, float], out_png: Path, out_pdf: Path, out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out_csv, index=False)

    fig, axes = plt.subplots(2, 1, figsize=(7.2, 8.0), dpi=300, sharey=True)
    fig.suptitle(
        "AUC: Task vs Fixation-Based Recentering",
        fontsize=16,
        fontweight="bold",
        fontfamily=FONT_FAMILY,
        y=0.95,
    )

    for i, phase in enumerate(PHASE_ORDER):
        ax = axes[i]
        phase_df = summary[summary["phase"] == phase]

        for scheme in SCHEME_ORDER:
            sch_df = phase_df[phase_df["scheme"] == scheme].sort_values("run")
            x = sch_df["run"].to_numpy(dtype=float)
            y = sch_df["mean"].to_numpy(dtype=float)
            lo = sch_df["ci_low"].to_numpy(dtype=float)
            hi = sch_df["ci_high"].to_numpy(dtype=float)

            ax.fill_between(x, lo, hi, color=SCHEME_COLORS[scheme], alpha=0.14, linewidth=0)
            ax.plot(
                x,
                y,
                color=SCHEME_COLORS[scheme],
                lw=LINE_WIDTH,
                marker="o",
                ms=3.2,
                label=scheme,
            )

        ax.axhline(0.5, color="black", lw=LINE_WIDTH, ls="--", zorder=0)
        ax.set_title(phase.title(), fontsize=13, fontweight="bold", fontfamily=FONT_FAMILY, pad=2)
        ax.set_ylabel("AUC", fontsize=13, fontweight="bold", fontfamily=FONT_FAMILY)
        ax.set_xlim(1, 8)
        ax.set_ylim(0.30, 1.01)
        ax.set_xticks(RUNS)
        ax.set_xticklabels(
            [run_to_session_label(int(r)) for r in RUNS],
            fontsize=11,
            fontfamily=FONT_FAMILY,
        )

        p = pvals.get(phase, np.nan)
        txt = f"p={p:.4f}" if np.isfinite(p) else "p=NA"
        ax.text(
            0.98,
            0.04,
            txt,
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=11,
            fontfamily=FONT_FAMILY,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.85"),
        )

        if i == 0:
            leg_handles = [
                plt.Line2D([0], [0], color=SCHEME_COLORS["FIX"], lw=LINE_WIDTH, marker="o", label="FIX"),
                plt.Line2D([0], [0], color=SCHEME_COLORS["TASK"], lw=LINE_WIDTH, marker="o", label="TASK"),
                plt.Line2D([0], [0], color="black", lw=LINE_WIDTH, ls="--", label="Chance"),
            ]
            ax.legend(
                handles=leg_handles,
                loc="lower left",
                frameon=True,
                prop={"family": FONT_FAMILY, "size": 10},
            )

        ax.grid(axis="both", linestyle="-", linewidth=0.6, alpha=0.2)
        ax.spines[["top", "right"]].set_visible(False)
        for side in ["left", "bottom"]:
            ax.spines[side].set_visible(True)
            ax.spines[side].set_color("black")
            ax.spines[side].set_linewidth(0.8)
        for tick in ax.get_yticklabels():
            tick.set_fontfamily(FONT_FAMILY)

    axes[-1].set_xlabel(
        "Session Number, Run Number",
        fontsize=13,
        fontweight="bold",
        fontfamily=FONT_FAMILY,
        labelpad=10,
    )

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
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
        "--subject-run-csv",
        type=Path,
        default=REPO_ROOT / "data" / "fig6_auc_subject_run.csv",
        help="Path to save per-subject per-run AUC values.",
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
        default=REPO_ROOT / "figures",
        help="Output figure directory.",
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=REPO_ROOT / "data" / "fig6_run_level_summary.csv",
        help="Path to save run-level summary table.",
    )
    args = parser.parse_args()

    _, df = compute_recentering_tables(args.of_dir, args.on_dir, max_runs=args.max_runs)
    args.subject_run_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.subject_run_csv, index=False)
    summary, pvals = _csv_summary(df)

    make_figure(
        summary,
        pvals=pvals,
        out_png=args.outdir / "fig6_auc_by_run_task_vs_fix.png",
        out_pdf=args.outdir / "fig6_auc_by_run_task_vs_fix.pdf",
        out_csv=args.summary_csv,
    )
    print("Saved Figure 6 outputs to", args.outdir)
    print("Saved subject-run AUC table to", args.subject_run_csv)


if __name__ == "__main__":
    main()
