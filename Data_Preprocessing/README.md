# Data Preprocessing Pipeline

This directory contains the scripts used for the preprocessing of fMRI and dMRI data, as well as the construction of structural and functional networks.

## Origin & Modifications

The pipeline is adapted from the **UKB-connectomics** pipeline originally developed by Sina Mansour et al.
* **Original Source:** [https://github.com/sina-mansour/UKB-connectomics](https://github.com/sina-mansour/UKB-connectomics)

We have modified the original scripts to specifically support the **Brainnetome Atlas (BNA)**, which integrates 246 cortical and subcortical regions into a unified parcellation scheme. Specific modifications include:
1.  **Unified Atlas Logic:** Simplified the workflow to handle the BNA as a single entity, removing the need to merge separate cortical and subcortical streams.
2.  **Custom Signal Extraction:** Added a Python script to extract BNA-based fMRI time series.
3.  **Tractography Optimization:** Adjusted MRtrix3 parameters for high-resolution probabilistic tractography (10M streamlines).

## File Descriptions

### 1. `preprocessing_pipeline.sh`
The master bash script that orchestrates the entire workflow for a single subject. It handles:
* Data organization and checking.
* Warping the BNA atlas from MNI space to individual native space.
* Calling the sub-scripts (listed below) to process fMRI and dMRI data.

### 2. `compute_BNA_fmri.py`
A Python script that:
* Loads the ICA-cleaned fMRI data.
* Loads the BNA atlas in native fMRI space.
* Extracts the mean time series for each of the 246 BNA regions.

### 3. `probabilistic_tractography_native_space.sh`
An MRtrix3-based script that performs:
* DWI preprocessing (conversion, denoising).
* Response function estimation (Dhollander algorithm).
* MSMT-CSD (Multi-Shell Multi-Tissue Constrained Spherical Deconvolution).
* Generation of 10 million probabilistic streamlines (iFOD2).
* SIFT2 re-weighting to improve quantitative accuracy.

### 4. `map_structural_connectivity.sh`
A script that:
* Maps the generated streamlines onto the native BNA atlas.
* Constructs connectivity matrices based on streamline counts, mean FA, and mean fiber length.
