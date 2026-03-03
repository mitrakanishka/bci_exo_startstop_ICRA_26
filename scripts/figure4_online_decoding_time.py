#!/usr/bin/env python3
"""Paper Figure 4: Online decoding-time distributions (onset/offset)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.patches as mpatches

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR))

from _shared.online_log_parser import load_all_online_logs
from _shared.stats import paired_wilcoxon


PHASE_ORDER = ["Onset", "Offset"]
SESSION_ORDER = [2, 3]
COLORS = {2: "#4C72B0", 3: "#55A868"}


def build_rt_table(df_trials: pd.DataFrame) -> pd.DataFrame:
    df = df_trials[df_trials["session"].isin([2, 3])].copy()

    onset = df[df["onset_label"] == "hit"][["subject", "session", "onset_rt_s"]].copy()
    onset = onset.rename(columns={"onset_rt_s": "rt_s"})
    onset["phase"] = "Onset"

    offset = df[(df["offset_attempted"]) & (df["offset_label"] == "hit")][
        ["subject", "session", "offset_rt_s"]
    ].copy()
    offset = offset.rename(columns={"offset_rt_s": "rt_s"})
    offset["rt_s"] = offset["rt_s"] + 2.0
    offset["phase"] = "Offset"

    out = pd.concat([onset, offset], ignore_index=True)
    out = out[np.isfinite(out["rt_s"])].copy()
    out = out[out["rt_s"] <= 6.0].copy()
    out["session"] = out["session"].astype(int)
    return out


def _paired_session_stats(df_rt: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for phase in PHASE_ORDER:
        med = (
            df_rt[df_rt["phase"] == phase]
            .groupby(["subject", "session"])["rt_s"]
            .median()
            .unstack("session")
            .reindex(columns=[2, 3])
            .dropna()
        )
        if med.empty:
            rows.append({"phase": phase, "n": 0, "w_stat": np.nan, "p_value": np.nan, "mean_delta_s": np.nan, "sd_delta_s": np.nan})
            continue

        d = med[3].to_numpy() - med[2].to_numpy()
        res = paired_wilcoxon(med[2].to_numpy(), med[3].to_numpy())
        sd = float(np.std(d, ddof=1)) if len(d) > 1 else np.nan
        rows.append(
            {
                "phase": phase,
                "n": int(res.n),
                "w_stat": float(res.statistic),
                "p_value": float(res.pvalue),
                "mean_delta_s": float(np.mean(d)),
                "sd_delta_s": sd,
            }
        )
    return pd.DataFrame(rows)


def make_figure(df_rt: pd.DataFrame, out_png: Path, out_pdf: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 3.4), dpi=300)

    centers = {"Onset": 1.0, "Offset": 2.0}
    offsets = {2: -0.13, 3: 0.13}
    width = 0.34

    for phase in PHASE_ORDER:
        for sess in SESSION_ORDER:
            vals_real = df_rt[(df_rt["phase"] == phase) & (df_rt["session"] == sess)]["rt_s"].to_numpy()
            if vals_real.size == 0:
                continue
            pos = centers[phase] + offsets[sess]
            vals_violin = vals_real
            if phase == "Offset":
                # Anchor offset violins at 0 s to match paper panel geometry.
                vals_violin = np.append(vals_real, 0.0)

            v = ax.violinplot(
                vals_violin,
                positions=[pos],
                widths=width,
                showmeans=False,
                showmedians=False,
                showextrema=False,
            )
            for body in v["bodies"]:
                body.set_facecolor(COLORS[sess])
                body.set_edgecolor("none")
                body.set_alpha(0.35)

            jitter = (np.random.RandomState(42).rand(vals_real.size) - 0.5) * 0.10
            ax.scatter(
                np.full(vals_real.size, pos) + jitter,
                vals_real,
                s=11,
                alpha=0.40,
                color=COLORS[sess],
                edgecolors="none",
            )

            m = float(np.mean(vals_real))
            ax.hlines(m, pos - width * 0.24, pos + width * 0.24, colors="black", linewidth=1.0)

    ax.set_xlim(0.55, 2.45)
    ax.set_ylim(0, 6)
    ax.set_yticks(np.arange(0, 7, 1))
    ax.set_xticks([1.0, 2.0])
    ax.set_xticklabels(PHASE_ORDER, fontsize=16, fontweight="bold")
    ax.set_ylabel("Decoding time (s)", fontsize=16, fontweight="bold")
    ax.set_title("Online Decoding Time", fontsize=22, fontweight="bold", pad=10)

    handles = [
        mpatches.Patch(color=COLORS[2], alpha=0.5, label="Session 2"),
        mpatches.Patch(color=COLORS[3], alpha=0.5, label="Session 3"),
    ]
    ax.legend(handles=handles, loc="upper right", frameon=True, fontsize=9)

    ax.grid(axis="y", linestyle="-", linewidth=1.0, alpha=0.32)
    ax.spines[["top", "right"]].set_visible(False)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_png)
    fig.savefig(out_pdf)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--log-root",
        type=Path,
        default=Path("BCI_course_EXP") / "online_python_log",
        help="Path to online_python_log directory.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("figures") / "paper",
        help="Output figure directory.",
    )
    parser.add_argument(
        "--stats-out",
        type=Path,
        default=Path("data") / "processed" / "fig4_session2_vs3_stats.csv",
        help="Path to save paired-session stats.",
    )
    args = parser.parse_args()

    df_trials = load_all_online_logs(args.log_root)
    df_rt = build_rt_table(df_trials)

    make_figure(
        df_rt=df_rt,
        out_png=args.outdir / "fig4_online_decoding_time.png",
        out_pdf=args.outdir / "fig4_online_decoding_time.pdf",
    )

    stats = _paired_session_stats(df_rt)
    args.stats_out.parent.mkdir(parents=True, exist_ok=True)
    stats.to_csv(args.stats_out, index=False)

    print("Saved Figure 4 outputs to", args.outdir)
    print(stats.to_string(index=False))


if __name__ == "__main__":
    main()
