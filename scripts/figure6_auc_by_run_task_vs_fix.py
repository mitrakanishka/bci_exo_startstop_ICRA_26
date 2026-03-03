#!/usr/bin/env python3
"""Paper Figure 6: AUC by run for TASK vs fixation-based recentering."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR))

from _shared.stats import paired_wilcoxon


PHASE_ORDER = ["ONSET", "OFFSET"]
SCHEME_ORDER = ["TASK", "FIX"]
SCHEME_COLORS = {"TASK": "#E6862A", "FIX": "#2979FF"}
RUNS = np.arange(1, 9)
N_SUBJ = 8
Z95 = 1.96 / np.sqrt(N_SUBJ)

# Paper-faithful run-level means and SDs from the final notebook figure.
ON_TASK_MEAN = np.array([0.716, 0.700, 0.673, 0.681, 0.401, 0.398, 0.427, 0.439])
ON_TASK_SD = np.array([0.333, 0.339, 0.340, 0.293, 0.300, 0.303, 0.345, 0.369])
ON_FIX_MEAN = np.array([0.898, 0.869, 0.838, 0.833, 0.936, 0.898, 0.840, 0.818])
ON_FIX_SD = np.array([0.075, 0.232, 0.149, 0.166, 0.073, 0.101, 0.179, 0.164])

OFF_TASK_MEAN = np.array([0.698, 0.576, 0.648, 0.688, 0.666, 0.575, 0.585, 0.517])
OFF_TASK_SD = np.array([0.349, 0.399, 0.342, 0.285, 0.305, 0.354, 0.408, 0.355])
OFF_FIX_MEAN = np.array([0.786, 0.777, 0.744, 0.786, 0.941, 0.881, 0.909, 0.833])
OFF_FIX_SD = np.array([0.287, 0.286, 0.281, 0.233, 0.124, 0.195, 0.079, 0.254])

PAPER_PVALS = {"ONSET": 0.0117, "OFFSET": 0.0251}


def run_to_session_label(r: int) -> str:
    return f"S2,R{r}" if r <= 4 else f"S3,R{r - 4}"


def _paper_summary() -> pd.DataFrame:
    rows = []
    for phase, t_mean, t_sd, f_mean, f_sd in [
        ("ONSET", ON_TASK_MEAN, ON_TASK_SD, ON_FIX_MEAN, ON_FIX_SD),
        ("OFFSET", OFF_TASK_MEAN, OFF_TASK_SD, OFF_FIX_MEAN, OFF_FIX_SD),
    ]:
        for run_idx, run in enumerate(RUNS):
            for scheme, means, sds in [("TASK", t_mean, t_sd), ("FIX", f_mean, f_sd)]:
                m = float(means[run_idx])
                sd = float(sds[run_idx])
                half = float(Z95 * sd)
                rows.append(
                    {
                        "phase": phase,
                        "scheme": scheme,
                        "run": int(run),
                        "mean": m,
                        "sd": sd,
                        "ci_low": m - half,
                        "ci_high": m + half,
                        "n": N_SUBJ,
                    }
                )
    return pd.DataFrame(rows)


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
    fig.suptitle("AUC by run: Task vs Fixation-based recentering", fontsize=20, fontweight="bold", y=0.99)

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
            ax.plot(x, y, color=SCHEME_COLORS[scheme], lw=2.0, marker="o", ms=3.2, label=scheme)

        ax.axhline(0.5, color="0.55", lw=1, ls="--", zorder=0)
        ax.set_title(phase.title(), fontsize=13, fontweight="bold", pad=2)
        ax.set_ylabel("AUC", fontsize=12, fontweight="bold")
        ax.set_xlim(1, 8)
        ax.set_ylim(0.30, 1.01)
        ax.set_xticks(RUNS)
        ax.set_xticklabels([run_to_session_label(int(r)) for r in RUNS], fontsize=11)

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
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.85"),
        )

        if i == 0:
            leg_handles = [
                plt.Line2D([0], [0], color=SCHEME_COLORS["TASK"], lw=2.0, marker="o", label="TASK"),
                plt.Line2D([0], [0], color=SCHEME_COLORS["FIX"], lw=2.0, marker="o", label="FIX"),
                plt.Line2D([0], [0], color="0.55", lw=1, ls="--", label="Chance"),
            ]
            ax.legend(handles=leg_handles, loc="lower left", frameon=True, fontsize=11)

        ax.grid(axis="both", linestyle="-", linewidth=0.6, alpha=0.2)
        ax.spines[["top", "right"]].set_visible(False)

    axes[-1].set_xlabel("Session number, Run Number", fontsize=12, fontweight="bold")

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_png)
    fig.savefig(out_pdf)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source",
        choices=["paper", "csv"],
        default="paper",
        help="Use paper-frozen run summaries (default) or recompute from input CSV.",
    )
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=Path("data") / "processed" / "fig6_auc_subject_run.csv",
        help="Input CSV with columns: subj, run, phase, scheme, auc (used when --source csv).",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("figures") / "paper",
        help="Output figure directory.",
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=Path("data") / "processed" / "fig6_run_level_summary.csv",
        help="Path to save run-level summary table.",
    )
    args = parser.parse_args()

    if args.source == "paper":
        summary = _paper_summary()
        pvals = PAPER_PVALS.copy()
    else:
        df = pd.read_csv(args.input_csv)
        summary, pvals = _csv_summary(df)

    make_figure(
        summary,
        pvals=pvals,
        out_png=args.outdir / "fig6_auc_by_run_task_vs_fix.png",
        out_pdf=args.outdir / "fig6_auc_by_run_task_vs_fix.pdf",
        out_csv=args.summary_csv,
    )
    print("Saved Figure 6 outputs to", args.outdir)
    print("Source mode:", args.source)


if __name__ == "__main__":
    main()
