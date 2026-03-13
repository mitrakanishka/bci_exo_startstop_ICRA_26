"""Utilities for parsing online start/stop log files."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

import pandas as pd


LOG_NAME_RE = re.compile(
    r"sub_(?P<sub>\d+)_session_(?P<sess>\d+)_run_(?P<run>\d+)_online_python_log\.txt$",
    re.IGNORECASE,
)

HEADER_RE = re.compile(
    r"Subject\s+(?P<sub>\d+)\s*,\s*Session\s+(?P<sess>\d+)\s*,\s*Run\s+(?P<run>\d+)",
    re.IGNORECASE,
)

TRIAL_RE = re.compile(
    r"Trial:\s*(?P<trial>-?\d+).*?"
    r"Onset\s+Success:\s*(?P<on_succ>-?\d+).*?"
    r"Onset\s+Decode\s*time:\s*(?P<on_rt>-?\d*\.?\d+)s.*?"
    r"Offset\s+Success:\s*(?P<off_succ>-?\d+).*?"
    r"Offset\s+Decode\s*time:\s*(?P<off_rt>-?\d*\.?\d+)s",
    re.IGNORECASE,
)


def _decode_label(code: int, phase: str) -> str:
    if code == 1:
        return "hit"
    if code == 0:
        return "miss"
    if code == -2:
        return "timeout"
    if phase == "offset" and code == -1:
        return "not_attempted"
    return "unknown"


def _clean_rt(rt: str) -> float:
    v = float(rt)
    return float("nan") if v < 0 else v


def parse_online_log_file(path: Path) -> pd.DataFrame:
    """Parse one `online_python_log` file to a tidy trial table."""
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()

    m = LOG_NAME_RE.search(path.name)
    sub_file = int(m.group("sub")) if m else None
    ses_file = int(m.group("sess")) if m else None
    run_file = int(m.group("run")) if m else None

    sub_head = ses_head = run_head = None
    for line in lines[:3]:
        mh = HEADER_RE.search(line)
        if mh:
            sub_head = int(mh.group("sub"))
            ses_head = int(mh.group("sess"))
            run_head = int(mh.group("run"))
            break

    subject = sub_file if sub_file is not None else sub_head
    session = ses_file if ses_file is not None else ses_head
    run = run_file if run_file is not None else run_head
    if subject is None or session is None or run is None:
        raise ValueError(f"Could not parse subject/session/run from {path}")

    rows = []
    for line in lines:
        mt = TRIAL_RE.search(line)
        if not mt:
            continue

        onset_code = int(mt.group("on_succ"))
        offset_code = int(mt.group("off_succ"))

        rows.append(
            {
                "subject": subject,
                "session": session,
                "run": run,
                "trial": int(mt.group("trial")),
                "onset_code": onset_code,
                "onset_label": _decode_label(onset_code, "onset"),
                "onset_rt_s": _clean_rt(mt.group("on_rt")),
                "offset_code": offset_code,
                "offset_label": _decode_label(offset_code, "offset"),
                "offset_rt_s": _clean_rt(mt.group("off_rt")),
                # Match original analysis notebook logic used for paper figures.
                "offset_attempted": bool(offset_code != -1),
                "file": str(path),
            }
        )

    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(["subject", "session", "run", "trial"]).reset_index(drop=True)


def load_all_online_logs(log_root: Path) -> pd.DataFrame:
    """Load all online logs under `BCI_Harmony_ExperimentalData/online_python_log`."""
    files: Iterable[Path] = sorted(log_root.glob("Sub_*/*.txt"))
    if not files:
        raise FileNotFoundError(f"No online python log files found under: {log_root}")

    frames = []
    for fp in files:
        parsed = parse_online_log_file(fp)
        if not parsed.empty:
            frames.append(parsed)

    if not frames:
        raise RuntimeError("No trials parsed from online logs.")

    return pd.concat(frames, ignore_index=True)
