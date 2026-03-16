# Real-Time Onset and Offset BCI Exoskeleton

[Project Website](https://mitrakanishka.github.io/publications/real-time-decoding-movement-onset-offset/)

Repository for the paper: **Real-Time Decoding of Movement Onset and Offset for Brain-Controlled Rehabilitation Exoskeleton**

This repository contains:

- starter workflows for loading and preprocessing offline EEG recordings
- mirrored MATLAB and Python figure-generation pipelines that regenerate the paper results directly from the experiment data
- figure-ready outputs written to `fig_data/` and `figures/`

The repository assumes the experimental dataset folder `BCI_Harmony_ExperimentalData/` is placed at the repository root and kept outside version control.

## Share and Clone

If you want to link this repository from a website, the main repository URL is:

```text
https://github.com/mitrakanishka9732/bci_exo_startstop_ICRA_26
```

Users can clone it with:

```bash
git clone https://github.com/mitrakanishka9732/bci_exo_startstop_ICRA_26.git
cd bci_exo_startstop_ICRA_26
```

Important: the repository does not include the local experimental dataset folder `BCI_Harmony_ExperimentalData/`, so anyone using the full analysis pipeline will also need access to that dataset separately.

## Repository Overview

```text
bci_exo_startstop_ICRA_26/
├── BCI_Harmony_ExperimentalData/   # local experiment data used by the analyses
├── starter/                        # entry-point MATLAB and Python starter workflows
├── scripts/                        # mirrored MATLAB and Python figure-generation code
├── fig_data/                       # generated tables and metadata used by the figures
├── figures/                        # generated paper figures
├── environment.yml                 # Python environment definition
└── LICENSE
```

## 1. Environment Setup

### Python environment

Create the Conda environment with:

```bash
conda env create -f environment.yml
conda activate BCI
```

### MATLAB requirements

The MATLAB starter workflows were tested with:

- MATLAB R2025a
- Signal Processing Toolbox

The repository already includes the local helper code needed by the MATLAB starter workflow inside `starter/functions/`, including the GDF-loading utilities used by the starter scripts.

## 2. Experimental Data

The primary experimental dataset lives in:

```text
BCI_Harmony_ExperimentalData/
```

Many subfolders are present in the full dataset, but only a subset are required for the workflows currently implemented in this repository.

### Data hierarchy used in this repo

```text
BCI_Harmony_ExperimentalData/
├── offline_data/        # offline .gdf runs used by the starter workflows and Figure 2
├── online_python_log/   # online session logs used by Figures 3 and 4
├── epoched_data_of/     # offline epoched .mat files used by Figures 5 and 6
└── epoched_data_on/     # online epoched .mat files used by Figures 5 and 6
```

### Required subfolders

#### `offline_data/`

This folder contains the subject-level offline recording runs. Each run is stored as a `.gdf` recording together with its associated `.log` file. In this repository, these recordings are used for:

- the MATLAB starter preprocessing workflow
- the Python starter notebook
- the grand-average offline spectrogram in Figure 2

Typical organization looks like:

```text
offline_data/
└── Sub_1/
    └── Sub_1_run_1_..._offline/
        ├── Sub_1_run_1_offline.gdf
        └── Sub_1_run_1_offline.log
```

#### Offline trigger labels

The event markers in the offline `.gdf` files are one of the least intuitive parts of the dataset. These trigger labels define the run structure and are the main reference points used to align the offline analyses.

| Trigger | Meaning | Interpretation in the task |
| --- | --- | --- |
| `1000` | Start of run | Marks the beginning of an offline recording run |
| `300` | Start of countdown | Beginning of the pre-task countdown period |
| `100` | Begin start-MI | Onset of the movement-onset motor imagery period |
| `150` | Robot starts moving | Exoskeleton movement onset |
| `500` | Stop-MI | Onset of the movement-offset or stop-imagery period |
| `550` | Robot stops | Exoskeleton movement offset |
| `900` | Rest | Rest marker; in these data it occurs at the same timepoint as `550`, so it is redundant for timing |
| `950` | Robot returns home | Return phase back to the home position |
| `2000` | End of run | Marks the end of the offline recording run |

A typical run therefore follows the sequence:

```text
1000 -> 300 -> 100 -> 150 -> 500 -> 550/900 -> 950 -> 2000
```

In practice:

- `300` is the anchor used to align the countdown period before the task
- `100` and `500` separate the onset and offset motor-imagery phases
- `150`, `550`, and `950` capture the corresponding exoskeleton movement phases
- `1000` and `2000` mark the run boundaries

These event codes are especially important for understanding the starter preprocessing script and the grand-average spectrogram in Figure 2, because those analyses are organized around the temporal structure defined by these triggers.

#### `online_python_log/`

This folder contains the online-session Python logs. These logs are parsed directly to compute:

- command-delivery summaries for Figure 3
- online decoding-time statistics for Figure 4

#### `epoched_data_of/`

This folder contains subject-level offline epoched MATLAB `.mat` files. These are used for the offline reference and bias analyses that feed:

- Figure 5
- Figure 6

#### `epoched_data_on/`

This folder contains subject-level online epoched MATLAB `.mat` files. These are paired with the offline epoched files to evaluate online task-vs-fix behavior in:

- Figure 5
- Figure 6


## 3. Starter Code

The `starter/` folder is designed to give a new user a compact, readable entry point into the dataset before moving on to the full group-level figure scripts.

```text
starter/
├── offline_preprocessing_starter.m
├── offline_preprocessing_starter.ipynb
└── functions/
```

### `starter/offline_preprocessing_starter.m`

This MATLAB workflow opens one offline run and walks through the core preprocessing steps used throughout the project:

- loading a `.gdf` file
- selecting the EEG channels used in the analysis
- removing auxiliary or unwanted channels
- filtering the EEG
- regressing out EOG activity
- plotting quick quality-control views
- generating a C3-centered task spectrogram aligned to the offline trigger structure

Run it from MATLAB with:

```matlab
cd('/path/to/bci_exo_startstop_ICRA_26');
run('starter/offline_preprocessing_starter.m');
```

### `starter/offline_preprocessing_starter.ipynb`

This notebook provides the same overall entry point in Python. It is useful for users who want to inspect one offline recording step by step and understand how the basic preprocessing and visualization pipeline is organized before running the paper-level scripts.

### `starter/functions/`

This folder contains the MATLAB helper code used by the starter workflows. In most cases, a new user can treat it as support code and focus on the entry-point scripts above.

## 4. Reproducing the Figures

The figure pipelines in `scripts/python/` and `scripts/matlab/` are both fully data-driven and read directly from `BCI_Harmony_ExperimentalData/`. Both language implementations write their outputs to the same locations:

- `fig_data/`
- `figures/`

### Python

```bash
python scripts/python/run_all_figures.py
```

### MATLAB

Run from the repository root:

```matlab
addpath('scripts/matlab');
run_all_figures
```

### Mirrored figure scripts

Python:

- `scripts/python/figure2_grand_avg_spectrogram.py`
- `scripts/python/figure3_online_command_delivery.py`
- `scripts/python/figure4_online_decoding_time.py`
- `scripts/python/figure5_bias_shift_vs_identity.py`
- `scripts/python/figure6_auc_by_run_task_vs_fix.py`

MATLAB:

- `scripts/matlab/figure2_grand_avg_spectrogram.m`
- `scripts/matlab/figure3_online_command_delivery.m`
- `scripts/matlab/figure4_online_decoding_time.m`
- `scripts/matlab/figure5_bias_shift_vs_identity.m`
- `scripts/matlab/figure6_auc_by_run_task_vs_fix.m`

## Citation

<!-- Citation information to be added -->

## Questions

For any questions, contact Kanishka Mitra at [mitra819@mit.edu](mailto:mitra819@mit.edu).
