# Control Energy Analysis

This directory contains scripts to calculate brain state transition energies using **Network Control Theory (NCT)**. 

**All control energy calculations in this study were implemented using the `nctpy` toolbox.**
* **Repository:** [https://github.com/LindenParkesLab/nctpy](https://github.com/LindenParkesLab/nctpy)
* **Reference:** Parkes, L., et al. (2024).

## Pipeline & Files

The analysis is divided into **Subject-Level Calculation** (computationally intensive) and **Group-Level Aggregation**.

### 1. Subject-Level Scripts
Run these scripts for each individual subject.

* **`compute_control_energy.py`**
    * Calculates the baseline **Global Transition Energy** (aveTE) using uniform control weights.
    * **Output:** Energy matrices for transition pairs.

* **`compute_altered_control_energy_byAdd1.py`**
    * Calculates the energy required when a specific region's control weight is increased.
    * **Output:** Perturbed energy matrices used to determine **Regional Energy Regulation Capacity (rERC)**.

### 2. Group-Level Scripts
Run these scripts after subject-level calculations are complete to generate summary CSVs for statistics.

* **`organize_energy_results.py`**
    * Aggregates outputs from `compute_control_energy.py`.
    * **Output:** `df_aveTE_...csv` (Average Transition Energy) and `df_global_stability_...csv`.

* **`calculate_rERC.py`**
    * Computes the **rERC** metric by comparing baseline energy and perturbed energy.
    * **Formula:** $rERC = (Baseline - Perturbed) / Baseline$.
    * **Output:** `df_rERC_...csv`.
