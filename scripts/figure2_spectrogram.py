#!/usr/bin/env python3
"""Paper Figure 2: Grand-average offline task spectrogram."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent.parent / ".matplotlib"))

import matplotlib as mpl
import matplotlib.pyplot as plt
import mne
import numpy as np
from scipy.signal import butter, filtfilt, welch

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parent
FONT_FAMILY = "DejaVu Sans"
KEEP_IDX = np.array([4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 16, 17, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29]) - 1
EVENT_SECONDS = np.array([0.2, 3.0, 6.0, 9.0, 12.0, 15.0, 17.0], dtype=float)
EVENT_LABELS = [
    "Countdown",
    "Begin MI",
    "Robot Moves",
    "End MI",
    "Robot Stops (rest)",
    "Rest (moves)",
    "Robot Returns",
]


def resolve_dataset_root() -> Path:
    for name in ("BCI_Harmony_ExperimentalData", "BCI_course_EXP"):
        path = REPO_ROOT / name
        if path.exists():
            return path
    raise FileNotFoundError("Could not find the dataset folder.")


def discover_subjects(offline_root: Path, requested: list[int] | None) -> list[int]:
    if requested:
        return requested
    subject_ids = []
    for path in sorted(offline_root.glob("Sub_*")):
        try:
            subject_ids.append(int(path.name.split("_")[1]))
        except (IndexError, ValueError):
            continue
    if not subject_ids:
        raise FileNotFoundError(f"No subject folders found under {offline_root}")
    return subject_ids


def extract_events(raw: mne.io.BaseRaw) -> tuple[np.ndarray, dict[str, int]]:
    events, event_id = mne.events_from_annotations(raw, verbose="ERROR")
    code_map = {desc: int(desc) for desc in event_id}
    event_codes = np.array([code_map[desc] for desc in raw.annotations.description], dtype=int)
    event_samples = np.rint(raw.annotations.onset * raw.info["sfreq"]).astype(int)
    event_zero = np.zeros_like(event_codes)
    return np.column_stack([event_samples, event_zero, event_codes]), code_map


def compute_subject_spectrogram(subject_id: int, offline_root: Path, channel: str) -> tuple[np.ndarray, dict]:
    gdf_files = sorted((offline_root / f"Sub_{subject_id}").glob("*offline*/*.gdf"))
    if not gdf_files:
        raise FileNotFoundError(f"No offline GDF files found for subject {subject_id}")

    subject_matrices = []
    kept_labels: list[str] | None = None
    fs_ref: float | None = None

    for gdf_path in gdf_files:
        raw = mne.io.read_raw_gdf(gdf_path, preload=True, verbose="ERROR")
        events, _ = extract_events(raw)

        fs = float(raw.info["sfreq"])
        if fs_ref is None:
            fs_ref = fs
        elif fs != fs_ref:
            raise ValueError(f"Inconsistent sampling rate for subject {subject_id}: {fs} vs {fs_ref}")

        data = raw.get_data().T
        stop_sample = int(events[-1, 0])
        data = data[:stop_sample, :]

        eeg = data[:, :64]
        eog = data[:, 64:67]

        b, a = butter(2, np.array([0.1, 45.0]) / (fs / 2.0), btype="bandpass")
        eeg = filtfilt(b, a, eeg, axis=0)
        if eog.size:
            eog = filtfilt(b, a, eog, axis=0)
            beta, *_ = np.linalg.lstsq(eog, eeg, rcond=None)
            eeg = eeg - eog @ beta

        eeg = eeg[:, KEEP_IDX]
        eeg = eeg - np.mean(eeg, axis=1, keepdims=True)

        if kept_labels is None:
            kept_labels = [raw.ch_names[idx] for idx in KEEP_IDX]

        try:
            channel_idx = next(i for i, name in enumerate(kept_labels) if name.upper() == channel.upper())
        except StopIteration as exc:
            raise ValueError(f"Channel {channel} is not in the kept montage: {kept_labels}") from exc

        countdown_pos = events[events[:, 2] == 300, 0]
        countdown_pos = countdown_pos[(countdown_pos - 2 * fs >= 0) & (countdown_pos + 18 * fs <= eeg.shape[0])]
        if countdown_pos.size == 0:
            continue

        step_samples = int(round(0.0625 * fs))
        trial_len = int(18 * (1 / 0.0625))
        baseline_len = int((2 - 1) / 0.0625 + 1)
        freq_keep = np.arange(8, 31)

        trial_psd = np.zeros((freq_keep.size, trial_len, countdown_pos.size), dtype=float)
        baseline_psd = np.zeros((freq_keep.size, baseline_len, countdown_pos.size), dtype=float)

        for trial_idx, start in enumerate(countdown_pos.astype(int)):
            curr = start
            for t in range(trial_len):
                freq_axis, pxx = welch(
                    eeg[curr : curr + int(fs), channel_idx],
                    fs=fs,
                    window="hamming",
                    nperseg=int(0.5 * fs),
                    noverlap=int(round(0.4 * fs)),
                    nfft=int(fs),
                )
                trial_psd[:, t, trial_idx] = pxx[(freq_axis >= 8) & (freq_axis <= 30)]
                curr += step_samples

            curr = start - int(2 * fs)
            for t in range(baseline_len):
                freq_axis, pxx = welch(
                    eeg[curr : curr + int(fs), channel_idx],
                    fs=fs,
                    window="hamming",
                    nperseg=int(0.5 * fs),
                    noverlap=int(round(0.4 * fs)),
                    nfft=int(fs),
                )
                baseline_psd[:, t, trial_idx] = pxx[(freq_axis >= 8) & (freq_axis <= 30)]
                curr += step_samples

        baseline_mean = baseline_psd.mean(axis=1)
        normalized = np.zeros_like(trial_psd)
        for trial_idx in range(trial_psd.shape[2]):
            normalized[:, :, trial_idx] = trial_psd[:, :, trial_idx] / baseline_mean[:, trial_idx][:, None]

        subject_matrices.append(np.flipud(np.mean(10 * np.log10(normalized), axis=2)))

    if not subject_matrices or kept_labels is None or fs_ref is None:
        raise RuntimeError(f"No valid spectrogram trials found for subject {subject_id}")

    subject_mean = np.mean(np.stack(subject_matrices, axis=2), axis=2)
    meta = {
        "subject_id": subject_id,
        "channel_name": channel.upper(),
        "num_windows": int(1 / 0.0625),
        "freq_plot": list(range(30, 7, -1)),
        "kept_labels": kept_labels,
        "n_runs": len(gdf_files),
    }
    return subject_mean, meta


def export_processed_data(matrix: np.ndarray, meta: dict, matrix_csv: Path, meta_json: Path) -> None:
    matrix_csv.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(matrix_csv, matrix, delimiter=",")
    with meta_json.open("w", encoding="utf-8") as handle:
        json.dump(meta, handle, indent=2)


def make_figure(matrix: np.ndarray, meta: dict, out_png: Path, out_pdf: Path) -> None:
    fig, ax = plt.subplots(figsize=(10.0, 7.5), dpi=300)
    image = ax.imshow(matrix, aspect="auto", cmap="jet", origin="upper")

    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label("Power (dB/Hz)", fontsize=12, fontfamily=FONT_FAMILY)

    event_x = np.rint(EVENT_SECONDS * float(meta["num_windows"])).astype(int)
    for x, label in zip(event_x, EVENT_LABELS):
        ax.axvline(x, color="black", linewidth=2.0)
        ax.text(
            x,
            -0.8,
            label,
            ha="center",
            va="bottom",
            fontsize=10,
            fontfamily=FONT_FAMILY,
            clip_on=False,
        )

    freq_plot = np.asarray(meta["freq_plot"], dtype=int)
    y_idx = np.arange(0, len(freq_plot), 2)
    ax.set_xticks(event_x)
    ax.set_xticklabels(["0", "3", "6", "9", "12", "15", "17"], fontfamily=FONT_FAMILY, fontsize=11)
    ax.set_yticks(y_idx)
    ax.set_yticklabels(freq_plot[y_idx], fontfamily=FONT_FAMILY, fontsize=11)
    ax.set_xlabel("Time (s)", fontsize=13, fontweight="bold", fontfamily=FONT_FAMILY)
    ax.set_ylabel("Frequency (Hz)", fontsize=13, fontweight="bold", fontfamily=FONT_FAMILY)
    ax.set_title(
        f"Grand Average Spectrogram, Electrode: {meta['channel_name']}",
        fontsize=15,
        fontweight="bold",
        fontfamily=FONT_FAMILY,
        pad=22,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_png)
    fig.savefig(out_pdf)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--channel", type=str, default="C3", help="Channel to plot from the reduced montage.")
    parser.add_argument(
        "--subjects",
        type=int,
        nargs="*",
        default=None,
        help="Optional subject IDs. Default is all offline subjects.",
    )
    parser.add_argument(
        "--matrix-csv",
        type=Path,
        default=REPO_ROOT / "data" / "processed" / "fig2_spectrogram_matrix.csv",
        help="Path to save the grand-average PSD matrix.",
    )
    parser.add_argument(
        "--meta-json",
        type=Path,
        default=REPO_ROOT / "data" / "processed" / "fig2_spectrogram_meta.json",
        help="Path to save figure metadata.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=REPO_ROOT / "figures" / "paper",
        help="Output figure directory.",
    )
    args = parser.parse_args()

    dataset_root = resolve_dataset_root()
    offline_root = dataset_root / "offline_data"
    subject_ids = discover_subjects(offline_root, args.subjects)

    subject_matrices = []
    kept_meta = None
    for subject_id in subject_ids:
        subject_matrix, subject_meta = compute_subject_spectrogram(subject_id, offline_root, args.channel)
        subject_matrices.append(subject_matrix)
        if kept_meta is None:
            kept_meta = subject_meta

    grand_average = np.mean(np.stack(subject_matrices, axis=2), axis=2)
    meta = {
        "channel_name": args.channel.upper(),
        "num_windows": kept_meta["num_windows"],
        "freq_plot": kept_meta["freq_plot"],
        "kept_labels": kept_meta["kept_labels"],
        "subject_ids": subject_ids,
        "n_subjects": len(subject_ids),
    }

    export_processed_data(grand_average, meta, args.matrix_csv, args.meta_json)
    make_figure(
        matrix=grand_average,
        meta=meta,
        out_png=args.outdir / "fig2_spectrogram.png",
        out_pdf=args.outdir / "fig2_spectrogram.pdf",
    )
    print("Saved Figure 2 outputs to", args.outdir)


if __name__ == "__main__":
    main()
