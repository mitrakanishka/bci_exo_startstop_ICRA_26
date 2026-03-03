#!/usr/bin/env python3
"""Paper Figure 3: Online command-delivery accuracy (onset/offset)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import PercentFormatter

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR))

from _shared.online_log_parser import load_all_online_logs


OUTCOMES = ["hit", "miss", "timeout"]
PHASES = ["Onset", "Offset"]
SESSIONS = [2, 3]
COLOR_MAP = {"hit": "#55A868", "miss": "#E17C05", "timeout": "#4C78A8"}


def _proportions(df: pd.DataFrame, label_col: str) -> pd.DataFrame:
    counts = (
        df.groupby(["subject", "session", label_col]).size().rename("n").reset_index()
    )
    totals = counts.groupby(["subject", "session"])["n"].sum().rename("total").reset_index()
    out = counts.merge(totals, on=["subject", "session"], how="left")
    out["prop"] = out["n"] / out["total"]
    return out


def _phase_session_means(df_trials: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    df_s23 = df_trials[df_trials["session"].isin(SESSIONS)].copy()

    on = _proportions(df_s23[df_s23["onset_label"].isin(OUTCOMES)], "onset_label")
    on = (
        on.groupby(["session", "onset_label"])["prop"]
        .mean()
        .rename("mean_prop")
        .reset_index()
        .rename(columns={"onset_label": "outcome"})
    )
    on["phase"] = "Onset"

    off_input = df_s23[df_s23["offset_attempted"]].copy()
    off_input = off_input[off_input["offset_label"].isin(OUTCOMES)]
    off = _proportions(off_input, "offset_label")
    off = (
        off.groupby(["session", "offset_label"])["prop"]
        .mean()
        .rename("mean_prop")
        .reset_index()
        .rename(columns={"offset_label": "outcome"})
    )
    off["phase"] = "Offset"

    return on, off


def make_figure(df_trials: pd.DataFrame, out_png: Path, out_pdf: Path, out_csv: Path) -> None:
    on, off = _phase_session_means(df_trials)
    mean_stack = pd.concat([on, off], ignore_index=True)

    grid = pd.MultiIndex.from_product(
        [SESSIONS, PHASES, OUTCOMES], names=["session", "phase", "outcome"]
    )
    mean_stack = (
        mean_stack.set_index(["session", "phase", "outcome"])["mean_prop"]
        .reindex(grid)
        .fillna(0.0)
        .rename("mean_prop")
        .reset_index()
    )
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    mean_stack.to_csv(out_csv, index=False)

    fig, ax = plt.subplots(figsize=(5.8, 3.4), dpi=300)

    w_bar = 0.42
    gap = 1.25
    x_on = 0.0
    x_off = x_on + gap

    xpos = {
        ("Onset", 2): x_on,
        ("Onset", 3): x_on + w_bar,
        ("Offset", 2): x_off,
        ("Offset", 3): x_off + w_bar,
    }

    def get_mean(session: int, phase: str, outcome: str) -> float:
        row = mean_stack[
            (mean_stack["session"] == session)
            & (mean_stack["phase"] == phase)
            & (mean_stack["outcome"] == outcome)
        ]
        if row.empty:
            return 0.0
        return float(row["mean_prop"].iloc[0])

    for phase in PHASES:
        for session in SESSIONS:
            x0 = xpos[(phase, session)]
            bottom = 0.0
            for outcome in OUTCOMES:
                h = get_mean(session, phase, outcome)
                ax.bar(
                    x0,
                    h,
                    width=w_bar * 0.9,
                    bottom=bottom,
                    color=COLOR_MAP[outcome],
                    edgecolor="black",
                    linewidth=0.5,
                )
                if bottom > 0:
                    ax.plot(
                        [x0 - (w_bar * 0.9) / 2.0, x0 + (w_bar * 0.9) / 2.0],
                        [bottom, bottom],
                        color="white",
                        lw=1.0,
                        zorder=3,
                    )
                if outcome == "hit" and h > 0.06:
                    ax.text(
                        x0,
                        bottom + 0.5 * h,
                        f"{h*100:.0f}%",
                        ha="center",
                        va="center",
                        fontsize=13,
                        color="white",
                        fontweight="bold",
                    )
                bottom += h

    ax.set_ylim(0.0, 1.0)
    ax.yaxis.set_major_formatter(PercentFormatter(1.0))
    ax.set_xticks([])
    ax.set_ylabel("")
    ax.set_title("Online command-delivery accuracy", fontsize=20, fontweight="bold", pad=12)

    xf = ax.get_xaxis_transform()
    for phase in PHASES:
        for session in SESSIONS:
            ax.text(
                xpos[(phase, session)],
                -0.01,
                f"Session {session}",
                ha="center",
                va="top",
                transform=xf,
                fontsize=11,
            )

    ax.text(
        (x_on + (x_on + w_bar)) / 2,
        -0.12,
        "Onset",
        ha="center",
        va="top",
        transform=xf,
        fontsize=13,
        fontweight="bold",
    )
    ax.text(
        (x_off + (x_off + w_bar)) / 2,
        -0.12,
        "Offset",
        ha="center",
        va="top",
        transform=xf,
        fontsize=13,
        fontweight="bold",
    )

    handles = [
        plt.Line2D([0], [0], color=COLOR_MAP[k], lw=8, label=k.capitalize()) for k in OUTCOMES
    ]
    ax.legend(handles=handles, loc="lower right", frameon=True, fontsize=9)

    ax.grid(axis="y", linestyle=":", linewidth=1.0, alpha=0.55)
    ax.spines[["top", "right"]].set_visible(False)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.subplots_adjust(bottom=0.22)
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
        "--csv-out",
        type=Path,
        default=Path("data") / "processed" / "fig3_group_composition.csv",
        help="Path to save aggregated composition table.",
    )
    args = parser.parse_args()

    df = load_all_online_logs(args.log_root)
    make_figure(
        df_trials=df,
        out_png=args.outdir / "fig3_online_command_delivery_accuracy.png",
        out_pdf=args.outdir / "fig3_online_command_delivery_accuracy.pdf",
        out_csv=args.csv_out,
    )
    print("Saved Figure 3 outputs to", args.outdir)


if __name__ == "__main__":
    main()
