# BCI Harmony Start/Stop Repository

This repository contains starter code and analysis scripts for working with the BCI Harmony experimental dataset. The codebase is organized to support both:

- beginner-friendly offline preprocessing and visualization
- reproducible figure generation directly from the underlying experiment data

The repository expects the experimental dataset folder, `BCI_Harmony_ExperimentalData/`, to live at the repository root. That folder is kept outside version control, while this repository stores the analysis code, starter workflows, and generated outputs.

## Repository Overview

```text
bci_exo_startstop_ICRA_26/
├── BCI_Harmony_ExperimentalData/   # local experiment data (not tracked by git)
├── starter/                        # starter MATLAB and Python workflows
├── scripts/                        # reproducible analysis and figure scripts
├── data/                           # generated tables and metadata
├── figures/                        # generated figures
└── environment.yml                 # conda environment definition
```

## 1. Environment Setup

This repository uses a Conda environment named `BCI`.

### Create the environment

```bash
conda env create -f environment.yml
conda activate BCI
```

### Update an existing environment

If the `BCI` environment already exists, update it with:

```bash
conda env update -n BCI -f environment.yml --prune
conda activate BCI
```

### Verify the environment

```bash
python -V
```

The environment is configured around Python `3.10`.

### Launch JupyterLab

For the Python notebook starter workflow:

```bash
jupyter lab
```

## 2. Experimental Data

The main dataset folder used by this repository is:

```text
BCI_Harmony_ExperimentalData/
```

Only a subset of the available subfolders are required for the code currently included in this repository.

### Data hierarchy

```text
BCI_Harmony_ExperimentalData/
├── offline_data/        # offline .gdf recordings used by the starter workflows and Figure 2
├── online_python_log/   # online log files used for Figures 3 and 4
├── epoched_data_of/     # offline epoched .mat files used for Figures 5 and 6
└── epoched_data_on/     # online epoched .mat files used for Figures 5 and 6
```

### Required subfolders

#### `offline_data/`

This folder contains the offline recording runs for each subject. These recordings are stored as `.gdf` files, along with associated log files. In this repository, `offline_data/` is used for:

- the MATLAB starter preprocessing script
- the Python starter notebook
- the grand-average spectrogram analysis for Figure 2

Typical contents are organized by subject and run, for example:

```text
offline_data/
└── Sub_1/
    └── Sub_1_run_1_..._offline/
        ├── Sub_1_run_1_offline.gdf
        └── Sub_1_run_1_offline.log
```

#### `online_python_log/`

This folder contains the online experiment logs generated during online BCI sessions. In this repository, these logs are used to compute:

- online command-delivery summaries for Figure 3
- online decoding-time analysis for Figure 4

#### `epoched_data_of/`

This folder contains subject-level offline epoched data stored in MATLAB `.mat` files. In this repository, these files are used to compute offline reference and bias-related analyses that contribute to:

- Figure 5
- Figure 6

#### `epoched_data_on/`

This folder contains subject-level online epoched data stored in MATLAB `.mat` files. These files are paired with the offline epoched data and are used for:

- task-vs-fix bias analysis
- per-run AUC analysis

### Important note on `data/`

The repository also contains a top-level folder named `data/`, but this is different from `BCI_Harmony_ExperimentalData/`.

- `BCI_Harmony_ExperimentalData/` contains the source experiment data
- `data/` contains generated outputs such as CSV summaries and metadata created by the analysis scripts

## 3. Starter Code

The `starter/` folder is intended to provide simple, readable entry points for working with one subject's offline data before moving into the larger group-level analysis scripts.

### Contents of `starter/`

```text
starter/
├── offline_preprocessing_starter.m
├── offline_preprocessing_starter.ipynb
├── extract_offline_trigger_labels.m
└── functions/
```

### `offline_preprocessing_starter.m`

This is the main MATLAB starter script. It is written as a script, not a function, and walks through a simple offline EEG preprocessing workflow using one subject's offline `.gdf` file.

The script includes commented sections that explain each step:

1. locate the repository and dataset
2. open one offline `.gdf` file
3. identify EEG and auxiliary channels
4. remove non-EEG channels and selected rejected channels
5. apply filtering
6. remove EOG-related activity by regression
7. generate basic sanity-check plots
8. generate a task spectrogram using the same general analysis style as the older MATLAB workflow

To run it from MATLAB:

```matlab
cd('/path/to/bci_exo_startstop_ICRA_26');
run('starter/offline_preprocessing_starter.m');
```

### `offline_preprocessing_starter.ipynb`

This notebook is the fully independent Python version of the offline starter workflow. It mirrors the same general processing flow as the MATLAB starter script, but does not depend on MATLAB.

It is useful if you want to:

- inspect one offline run interactively
- understand the preprocessing pipeline step by step
- work entirely in Python

Open it with JupyterLab after activating the Conda environment:

```bash
jupyter lab
```

### `extract_offline_trigger_labels.m`

This MATLAB script scans all offline `.gdf` files and exports:

- a per-run trigger inventory
- a global list of unique trigger codes

The outputs are written to:

- `data/offline_trigger_inventory.csv`
- `data/offline_unique_trigger_labels.txt`

### `starter/functions/`

This folder contains MATLAB helper functions used by the starter MATLAB workflows. It includes the signal-processing utilities and supporting toolbox code needed by the scripts in `starter/`.

In most cases, new users do not need to edit anything inside `starter/functions/` to get started.

## Reproducing the Python Figures

All current figure scripts in `scripts/` are designed to run directly from the dataset and generate outputs into the repository-level `data/` and `figures/` folders.

Run all figure scripts with:

```bash
python scripts/run_all_figures.py
```

Current figure scripts include:

- `scripts/figure2_grand_avg_spectrogram.py`
- `scripts/figure3_online_command_delivery.py`
- `scripts/figure4_online_decoding_time.py`
- `scripts/figure5_bias_shift_vs_identity.py`
- `scripts/figure6_auc_by_run_task_vs_fix.py`

## Citation

<!-- Citation information to be added -->

## Questions

For any questions, contact Kanishka Mitra at [mitra819@mit.edu](mailto:mitra819@mit.edu).
