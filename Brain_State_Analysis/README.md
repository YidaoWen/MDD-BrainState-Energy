# Brain State Extraction and Temporal Dynamics Analysis

This directory contains the analysis pipeline for identifying recurrent brain states and characterizing their temporal dynamics in individuals with Major Depressive Disorder (MDD) and Healthy Controls (HC).

## Overview

The analysis is performed in two main steps:
1.  **State Extraction:** Defining discrete brain states based on the dominant activation of intrinsic networks (Yeo's 7 networks) at each time point.
2.  **Dynamics Calculation:** Computing temporal metrics (Fractional Occupancy, Dwell Time, Appearance Rate, Transition Probability) to quantify brain state stability and flexibility.

## Files
* `extract_brain_states.py`: Script to assign each fMRI time point to a dominant brain network state.
* `calculate_temporal_metrics.py`: Script to compute FO, DT, AR, and TP metrics from the extracted state sequences.

## Methodological Details

### 1. Brain State Extraction (`extract_brain_states.py`)
We utilized a **data-driven non-binary approach** tailored to individual brain activity.
* **Atlas:** Brainnetome Atlas (BNA) mapped to Yeo's 7 resting-state networks.
* **Method:** For each TR (Time Repetition), the BOLD signal is averaged within each of the 7 networks. The TR is assigned to the network with the highest mean amplitude.
* **Networks:** Visual (VIS), Somatomotor (SOM), Dorsal Attention (DAT), Ventral Attention (VAT), Limbic (LIM), Frontoparietal (FPN), Default Mode (DMN).

### 2. Temporal Metrics (`calculate_temporal_metrics.py`)
Four key metrics were calculated to characterize state dynamics:
* **Fractional Occupancy (FO):** Proportion of total scan time spent in a specific state.
* **Dwell Time (DT):** Average duration (in seconds) a state is maintained before switching.
* **Appearance Rate (AR):** Frequency of state occurrence per minute.
* **Transition Probability (TP):** Probability of switching from one state to another.

## Usage

### Prerequisites
* Python 3.x
* Required libraries: `numpy`, `pandas`, `scipy`
* Input data: Preprocessed fMRI BOLD time series (e.g., CSV files) and subject lists.

### Step 1: Extract States
Modify the paths in `extract_brain_states.py` to point to your data directory and run:
python extract_brain_states.py

### Step 2: Calculate Metrics

Once the dominance network labels are generated, run the second script to compute the dynamics:
python calculate_temporal_metrics.py

## References
- State extraction method adapted from: https://github.com/NeuroenergeticsLab/control_costs
- Dynamics calculation logic adapted from: https://github.com/singlesp/energy_landscape
