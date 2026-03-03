#!/usr/bin/env python3
"""Run all paper figure scripts (Fig. 3-6)."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


SCRIPTS = [
    "figure3_online_command_delivery.py",
    "figure4_online_decoding_time.py",
    "figure5_bias_shift_vs_identity.py",
    "figure6_auc_by_run_task_vs_fix.py",
]


def main() -> None:
    root = Path(__file__).resolve().parent
    for script in SCRIPTS:
        path = root / script
        print(f"\n=== Running {script} ===")
        subprocess.run([sys.executable, str(path)], check=True)


if __name__ == "__main__":
    main()
