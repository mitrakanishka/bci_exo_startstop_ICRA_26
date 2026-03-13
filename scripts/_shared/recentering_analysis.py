"""Shared helpers for recomputing recentering analyses from epoched .mat data."""

from __future__ import annotations

from pathlib import Path
import re

import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.stats import rankdata

FIXATION_CANDIDATES = (
    "on_fix_cell",
    "on_fixation_cell",
    "on_preCue_cell",
    "on_pre_fix_cell",
    "on_rest_cell",
    "on_countdown_cell",
)

ONLINE_POSITIVE_FIELDS = {
    "ONSET": ("on_startMI_cell", "on_bMI_cell"),
    "OFFSET": ("on_stopMI_cell", "on_eMI_cell"),
}

OFFLINE_NEGATIVE_FIELDS = {"ONSET": "data_rest_cell", "OFFSET": "data_dMI_cell"}
OFFLINE_POSITIVE_FIELDS = {"ONSET": "data_bMI_cell", "OFFSET": "data_eMI_cell"}
OFFLINE_PROTO0_FIELDS = {"ONSET": "data_rest_cell", "OFFSET": "data_dMI_cell"}
OFFLINE_PROTO1_FIELDS = {"ONSET": "data_bMI_cell", "OFFSET": "data_eMI_cell"}
METRIC_ORDER = ("Delta_pos", "Delta_neg", "Delta_sep")
PHASE_ORDER = ("ONSET", "OFFSET")


def _parse_subject_ids(paths: list[Path]) -> list[int]:
    ids: set[int] = set()
    for path in paths:
        match = re.search(r"subj_(\d+)_", path.name)
        if match:
            ids.add(int(match.group(1)))
    return sorted(ids)


def discover_overlapping_subjects(of_dir: Path, on_dir: Path) -> list[int]:
    of_ids = _parse_subject_ids(sorted(of_dir.glob("subj_*_epoched_data_of.mat")))
    on_ids = _parse_subject_ids(sorted(on_dir.glob("subj_*_epoched_data_on.mat")))
    subs = sorted(set(of_ids).intersection(on_ids))
    if not subs:
        raise FileNotFoundError("No overlapping epoched-data subjects found.")
    return subs


def load_subject_pair(of_dir: Path, on_dir: Path, subject_id: int) -> tuple[dict, dict]:
    of_path = of_dir / f"subj_{subject_id:03d}_epoched_data_of.mat"
    on_path = on_dir / f"subj_{subject_id:03d}_epoched_data_on.mat"
    if not of_path.exists() or not on_path.exists():
        raise FileNotFoundError(f"Missing epoched data for subject {subject_id:03d}")
    return (
        loadmat(of_path, squeeze_me=False, struct_as_record=False),
        loadmat(on_path, squeeze_me=False, struct_as_record=False),
    )


def _to_epoch_tensor(value: object) -> np.ndarray | None:
    if not isinstance(value, np.ndarray) or value.size == 0:
        return None
    X = np.asarray(value, dtype=float)
    if X.ndim == 2:
        X = X[:, :, None]
    if X.ndim != 3:
        return None
    return X


def get_first_existing(data: dict, names: tuple[str, ...]) -> np.ndarray | None:
    for name in names:
        if name in data:
            return data[name]
    return None


def pick3d(cell_arr: np.ndarray | None, run_index: int) -> np.ndarray | None:
    if cell_arr is None or run_index - 1 >= cell_arr.shape[0]:
        return None
    return _to_epoch_tensor(cell_arr[run_index - 1, 0])


def pool_cell3d(data: dict, field_name: str) -> np.ndarray | None:
    if field_name not in data:
        return None
    chunks: list[np.ndarray] = []
    for row_idx in range(data[field_name].shape[0]):
        X = _to_epoch_tensor(data[field_name][row_idx, 0])
        if X is not None:
            chunks.append(X)
    return np.concatenate(chunks, axis=2) if chunks else None


def ensure_spd(A: np.ndarray) -> np.ndarray:
    A = np.asarray(A, dtype=float)
    A = 0.5 * (A + A.T)
    n = A.shape[0]
    epsw = 1e-9 * np.trace(A) / max(n, 1)
    if not np.isfinite(epsw) or epsw <= 0:
        epsw = 1e-9
    return A + epsw * np.eye(n)


def covs_from_3d(X3: np.ndarray | None) -> np.ndarray | None:
    if X3 is None or X3.size == 0:
        return None
    T, C, N = X3.shape
    denom = max(T - 1, 1)
    Cstack = np.empty((C, C, N), dtype=float)
    for idx in range(N):
        Xi = X3[:, :, idx]
        Cstack[:, :, idx] = ensure_spd((Xi.T @ Xi) / denom)
    return Cstack


def cat_covs(A: np.ndarray | None, B: np.ndarray | None) -> np.ndarray | None:
    if A is None:
        return B
    if B is None:
        return A
    return np.concatenate([A, B], axis=2)


def covs_from_cellarr(cell_arr: np.ndarray) -> np.ndarray | None:
    Cstack = None
    for row_idx in range(cell_arr.shape[0]):
        X = _to_epoch_tensor(cell_arr[row_idx, 0])
        if X is not None:
            Cstack = cat_covs(Cstack, covs_from_3d(X))
    return Cstack


def logeu_mean(Cstack: np.ndarray) -> np.ndarray:
    n_cov = Cstack.shape[2]
    cdim = Cstack.shape[0]
    Lsum = np.zeros((cdim, cdim), dtype=float)
    for idx in range(n_cov):
        Ci = ensure_spd(Cstack[:, :, idx])
        lam, V = np.linalg.eigh(Ci)
        lam = np.maximum(np.real(lam), 1e-12)
        Lsum += V @ np.diag(np.log(lam)) @ V.T
    Lbar = Lsum / n_cov
    lam, V = np.linalg.eigh(0.5 * (Lbar + Lbar.T))
    return ensure_spd(V @ np.diag(np.exp(np.real(lam))) @ V.T)


def ref_from_3d(X3: np.ndarray) -> np.ndarray:
    C = covs_from_3d(X3)
    if C is None:
        return np.eye(X3.shape[1])
    return logeu_mean(C)


def spd_log(A: np.ndarray) -> np.ndarray:
    A = ensure_spd(A)
    lam, V = np.linalg.eigh(A)
    lam = np.maximum(np.real(lam), 1e-12)
    L = V @ np.diag(np.log(lam)) @ V.T
    return 0.5 * (L + L.T)


def spd_exp(L: np.ndarray) -> np.ndarray:
    L = 0.5 * (L + L.T)
    lam, V = np.linalg.eigh(L)
    return ensure_spd(V @ np.diag(np.exp(np.real(lam))) @ V.T)


def spd_dist(A: np.ndarray, B: np.ndarray) -> float:
    A = ensure_spd(A)
    B = ensure_spd(B)
    lam, V = np.linalg.eigh(B)
    lam = np.maximum(np.real(lam), 1e-12)
    Binv2 = V @ np.diag(1.0 / np.sqrt(lam)) @ V.T
    BA = Binv2 @ A @ Binv2.T
    M = 0.5 * (BA + BA.T)
    ev = np.linalg.eigvalsh(M)
    ev = np.maximum(np.real(ev), 1e-12)
    return float(np.sqrt(np.sum(np.log(ev) ** 2)))


def distances_and_margins(
    X3: np.ndarray | None,
    P0: np.ndarray,
    P1: np.ndarray,
    R: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if X3 is None or X3.size == 0:
        return np.array([]), np.array([]), np.array([])

    P0 = ensure_spd(P0)
    P1 = ensure_spd(P1)
    R = ensure_spd(R)
    lam, V = np.linalg.eigh(R)
    lam = np.maximum(np.real(lam), 1e-12)
    Rminv = V @ np.diag(1.0 / np.sqrt(lam)) @ V.T

    Cstack = covs_from_3d(X3)
    n_cov = Cstack.shape[2]
    d0 = np.empty(n_cov, dtype=float)
    d1 = np.empty(n_cov, dtype=float)
    for idx in range(n_cov):
        Ci = ensure_spd(Cstack[:, :, idx])
        Cw = ensure_spd(Rminv @ Ci @ Rminv.T)
        d0[idx] = spd_dist(Cw, P0)
        d1[idx] = spd_dist(Cw, P1)
    return d0, d1, d0 - d1


def med_nan(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return float(np.median(x)) if x.size else np.nan


def auc_from_scores(s_pos: np.ndarray, s_neg: np.ndarray) -> float:
    s_pos = np.asarray(s_pos, dtype=float)
    s_neg = np.asarray(s_neg, dtype=float)
    s_pos = s_pos[np.isfinite(s_pos)]
    s_neg = s_neg[np.isfinite(s_neg)]
    n_pos = s_pos.size
    n_neg = s_neg.size
    if n_pos == 0 or n_neg == 0:
        return np.nan

    ranks = rankdata(np.concatenate([s_pos, s_neg]), method="average")
    ranks_pos = ranks[:n_pos]
    return float((np.sum(ranks_pos) - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def build_R_rest_global_from_dir(of_dir: Path) -> np.ndarray:
    Cstack = None
    cguess = None
    for path in sorted(of_dir.glob("subj_*_epoched_data_of.mat")):
        data = loadmat(path, squeeze_me=False, struct_as_record=False)
        if "data_rest_cell" not in data:
            continue
        Cstack = cat_covs(Cstack, covs_from_cellarr(data["data_rest_cell"]))
        if cguess is None:
            X0 = _to_epoch_tensor(data["data_rest_cell"][0, 0])
            if X0 is not None:
                cguess = X0.shape[1]
    if Cstack is None:
        return np.eye(cguess or 32)
    return logeu_mean(Cstack)


def get_fix_for_run(online_data: dict, offline_data: dict, run_index: int) -> np.ndarray | None:
    for name in FIXATION_CANDIDATES:
        if name in online_data:
            X = pick3d(online_data[name], run_index)
            if X is not None:
                return X
    if "data_rest_cell" in offline_data:
        X = pick3d(offline_data["data_rest_cell"], run_index)
        if X is not None:
            return X
    return pool_cell3d(offline_data, "data_rest_cell")


def build_R_fixation_robust(
    online_data: dict,
    offline_data: dict,
    run_index: int,
    Rrest_global: np.ndarray,
    alpha: float = 0.30,
    lamb: float = 1e-3,
    keep_frac: float = 0.90,
) -> np.ndarray:
    X_fix = get_fix_for_run(online_data, offline_data, run_index)
    Cfix = covs_from_3d(X_fix)
    if Cfix is None:
        return ensure_spd(Rrest_global)

    Mu_fix = logeu_mean(Cfix)
    d = np.array([spd_dist(Cfix[:, :, idx], Mu_fix) for idx in range(Cfix.shape[2])], dtype=float)
    keep_n = max(1, round(keep_frac * Cfix.shape[2]))
    Cfix_kept = Cfix[:, :, np.argsort(d)[:keep_n]]

    Rfix_run = logeu_mean(Cfix_kept)
    L_mix = (1.0 - alpha) * spd_log(Rrest_global) + alpha * spd_log(Rfix_run)
    Rfix = spd_exp(L_mix)
    cdim = Rfix.shape[0]
    Rfix = (1.0 - lamb) * Rfix + lamb * (np.trace(Rfix) / cdim) * np.eye(cdim)
    return ensure_spd(Rfix)


def compute_recentering_tables(
    of_dir: Path,
    on_dir: Path,
    max_runs: int = 8,
    alpha: float = 0.30,
    lamb: float = 1e-3,
    keep_frac: float = 0.90,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    subjects = discover_overlapping_subjects(of_dir, on_dir)
    Rrest_global = build_R_rest_global_from_dir(of_dir)

    bias_rows: list[dict[str, float | int | str]] = []
    auc_rows: list[dict[str, float | int | str]] = []

    for subject_id in subjects:
        offline_data, online_data = load_subject_pair(of_dir, on_dir, subject_id)

        P0 = {
            phase: logeu_mean(covs_from_cellarr(offline_data[OFFLINE_PROTO0_FIELDS[phase]]))
            for phase in PHASE_ORDER
        }
        P1 = {
            phase: logeu_mean(covs_from_cellarr(offline_data[OFFLINE_PROTO1_FIELDS[phase]]))
            for phase in PHASE_ORDER
        }

        online_pos = {
            phase: get_first_existing(online_data, ONLINE_POSITIVE_FIELDS[phase])
            for phase in PHASE_ORDER
        }
        offline_neg = {
            phase: pool_cell3d(offline_data, OFFLINE_NEGATIVE_FIELDS[phase])
            for phase in PHASE_ORDER
        }

        acc = {(phase, metric): [] for phase in PHASE_ORDER for metric in METRIC_ORDER}

        for run_index in range(1, max_runs + 1):
            R_fix = build_R_fixation_robust(
                online_data,
                offline_data,
                run_index,
                Rrest_global,
                alpha=alpha,
                lamb=lamb,
                keep_frac=keep_frac,
            )

            for phase in PHASE_ORDER:
                X_pos = pick3d(online_pos[phase], run_index)
                X_neg = offline_neg[phase]
                if X_pos is None or X_neg is None:
                    continue

                R_task = ref_from_3d(X_pos)
                I = np.eye(P0[phase].shape[0])

                _, _, s_pos_id = distances_and_margins(X_pos, P0[phase], P1[phase], I)
                _, _, s_neg_id = distances_and_margins(X_neg, P0[phase], P1[phase], I)
                _, _, s_pos_task = distances_and_margins(X_pos, P0[phase], P1[phase], R_task)
                _, _, s_neg_task = distances_and_margins(X_neg, P0[phase], P1[phase], R_task)
                _, _, s_pos_fix = distances_and_margins(X_pos, P0[phase], P1[phase], R_fix)
                _, _, s_neg_fix = distances_and_margins(X_neg, P0[phase], P1[phase], R_fix)

                delta_pos = med_nan(s_pos_task) - med_nan(s_pos_id)
                delta_neg = med_nan(s_neg_task) - med_nan(s_neg_id)
                delta_sep = delta_pos - delta_neg

                acc[(phase, "Delta_pos")].append(delta_pos)
                acc[(phase, "Delta_neg")].append(delta_neg)
                acc[(phase, "Delta_sep")].append(delta_sep)

                auc_rows.extend(
                    [
                        {
                            "subj": subject_id,
                            "run": run_index,
                            "phase": phase,
                            "scheme": "TASK",
                            "auc": auc_from_scores(s_pos_task, s_neg_task),
                        },
                        {
                            "subj": subject_id,
                            "run": run_index,
                            "phase": phase,
                            "scheme": "FIX",
                            "auc": auc_from_scores(s_pos_fix, s_neg_fix),
                        },
                    ]
                )

        for phase in PHASE_ORDER:
            for metric in METRIC_ORDER:
                bias_rows.append(
                    {
                        "subj": subject_id,
                        "phase": phase,
                        "metric": metric,
                        "value": float(np.nanmean(acc[(phase, metric)])),
                    }
                )

    bias_df = pd.DataFrame(bias_rows).sort_values(["subj", "phase", "metric"]).reset_index(drop=True)
    auc_df = pd.DataFrame(auc_rows).sort_values(["subj", "run", "phase", "scheme"]).reset_index(drop=True)
    return bias_df, auc_df
