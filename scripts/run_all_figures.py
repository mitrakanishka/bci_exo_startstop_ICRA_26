#!/usr/bin/env python3
"""Run all paper figure scripts (Fig. 1-4)."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


SCRIPTS = [
    "figure1_bias_shift_vs_identity.py",
    "figure2_auc_by_run_task_vs_fix.py",
    "figure3_online_command_delivery.py",
    "figure4_online_decoding_time.py",
]


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent
    for script in SCRIPTS:
        path = script_dir / script
        print(f"\n=== Running {script} ===")
        subprocess.run([sys.executable, str(path)], check=True, cwd=repo_root)


if __name__ == "__main__":
    main()
