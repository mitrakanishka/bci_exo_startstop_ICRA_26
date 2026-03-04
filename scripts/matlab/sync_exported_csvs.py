#!/usr/bin/env python3
"""Map MATLAB EXPORT_FIG_DATA outputs into this repo's processed-data paths.

Expected input files (created by EXPORT_FIG_DATA.m):
  fig_data/fig1_bias_task_vs_identity.csv
  fig_data/fig2_auc_subject_run.csv

Mapped outputs:
  data/processed/fig1_bias_task_vs_identity.csv
  data/processed/fig2_auc_subject_run.csv
"""

from __future__ import annotations

from pathlib import Path
import shutil


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    src_dir = root / "fig_data"
    dst_dir = root / "data" / "processed"
    dst_dir.mkdir(parents=True, exist_ok=True)

    mapping = {
        src_dir / "fig1_bias_task_vs_identity.csv": dst_dir / "fig1_bias_task_vs_identity.csv",
        src_dir / "fig2_auc_subject_run.csv": dst_dir / "fig2_auc_subject_run.csv",
    }

    missing = [str(p) for p in mapping if not p.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing expected MATLAB export files:\n  " + "\n  ".join(missing)
        )

    for src, dst in mapping.items():
        shutil.copy2(src, dst)
        print(f"Copied {src} -> {dst}")


if __name__ == "__main__":
    main()
