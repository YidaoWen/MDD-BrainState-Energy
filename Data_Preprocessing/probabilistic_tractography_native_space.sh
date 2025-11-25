#!/bin/bash

# ==============================================================================
# Probabilistic Tractography Script
# ==============================================================================
# Modified from UKB-connectomics.
# This script performs DWI preprocessing and generates streamlines using MRtrix3.

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

main_dir=$1
ukb_subjects_dir=$2
ukb_subject_id=$3
ukb_instance=$4
streamlines=$5
num_threads=$6

# Ensure MRtrix bin is in PATH or define it here if necessary
# mrtrix_dir="/path/to/mrtrix3/bin" # Uncomment and set if not in environment path
# Command prefix (assumes mrtrix commands are in PATH, otherwise prepend $mrtrix_dir/)
MRTRIX_CMD="" 

script_dir="${main_dir}/scripts"
template_dir="${main_dir}/data/templates"
temporary_dir="${main_dir}/data/temporary"

fsaverage_dir="${template_dir}/freesurfer/fsaverage"
dmri_dir="${ukb_subjects_dir}/${ukb_subject_id}_${ukb_instance}/dMRI/dMRI"

# Set threading for MRtrix commands
threading="-nthreads ${num_threads}" 

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

echo -e "${GREEN}[INFO]${NC} `date`: Starting tractography for: ${ukb_subject_id}_${ukb_instance}"

cd "${dmri_dir}"

# Create a temporary directory to store files
tractography_dir="${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/tractography"
if [ ! -d "${tractography_dir}" ]; then
    mkdir -p "${tractography_dir}"
fi

# ------------------------------------------------------------------------------
# 1. Preprocessing & Conversion
# ------------------------------------------------------------------------------

# Convert the initial diffusion image to .mif
dwi_mif="${dmri_dir}/dwi.mif"
if [ ! -f ${dwi_mif} ]; then
    echo -e "${GREEN}[INFO]${NC} `date`: Converting dwi image to mif"
    mrconvert "${dmri_dir}/data_ud.nii.gz" "${dwi_mif}" \
              -fslgrad "${dmri_dir}/bvecs" "${dmri_dir}/bvals" \
              -datatype float32 -strides 0,0,0,1 ${threading} -info
fi

# Extract mean B0 image
dwi_meanbzero="${dmri_dir}/dwi_meanbzero.mif"
dwi_meanbzero_nii="${dmri_dir}/dwi_meanbzero.nii.gz"
if [ ! -f ${dwi_meanbzero} ]; then
    echo -e "${GREEN}[INFO]${NC} `date`: Extracting mean B0 image"
    dwiextract ${threading} -info "${dwi_mif}" -bzero - | mrmath ${threading} -info - mean -axis 3 "${dwi_meanbzero}"
    mrconvert "${dwi_meanbzero}" "${dwi_meanbzero_nii}" ${threading} -info
fi

# Create a dwi brain mask using FSL BET
dwi_meanbzero_brain="${dmri_dir}/dwi_meanbzero_brain.nii.gz"
dwi_meanbzero_brain_mask="${dmri_dir}/dwi_meanbzero_brain_mask.nii.gz"
if [ ! -f ${dwi_meanbzero_brain_mask} ]; then
    echo -e "${GREEN}[INFO]${NC} `date`: Computing dwi brain mask (BET)"
    bet "${dwi_meanbzero_nii}" "${dwi_meanbzero_brain}" -m -R -f 0.2 -g -0.05
fi

# ------------------------------------------------------------------------------
# 2. Response Function & CSD
# ------------------------------------------------------------------------------

# Estimate response function using dhollander method
wm_txt="${dmri_dir}/wm.txt"
gm_txt="${dmri_dir}/gm.txt"
csf_txt="${dmri_dir}/csf.txt"
if [ ! -f ${wm_txt} ] || [ ! -f ${gm_txt} ] || [ ! -f ${csf_txt} ]; then
    echo -e "${GREEN}[INFO]${NC} `date`: Estimation of response function using dhollander"
    dwi2response dhollander "${dwi_mif}" "${wm_txt}" "${gm_txt}" "${csf_txt}" \
                            -voxels "${dmri_dir}/voxels.mif" ${threading} -info
fi

# Multi-Shell, Multi-Tissue Constrained Spherical Deconvolution (MSMT-CSD)
wm_fod="${dmri_dir}/wmfod.mif"
gm_fod="${dmri_dir}/gmfod.mif"
csf_fod="${dmri_dir}/csffod.mif"
dwi_mask_dilated="${dmri_dir}/dwi_meanbzero_brain_mask_dilated_2.nii.gz"

if [ ! -f ${wm_fod} ]; then
    echo -e "${GREEN}[INFO]${NC} `date`: Running MSMT-CSD"
    # Create dilated mask
    maskfilter -npass 2 "${dwi_meanbzero_brain_mask}" dilate "${dwi_mask_dilated}" ${threading} -info
    # Perform CSD
    dwi2fod msmt_csd "${dwi_mif}" -mask "${dwi_mask_dilated}" "${wm_txt}" "${wm_fod}" \
            "${gm_txt}" "${gm_fod}" "${csf_txt}" "${csf_fod}" ${threading} -info
fi

# Intensity Normalisation (mtnormalise)
wm_fod_norm="${dmri_dir}/wmfod_norm.mif"
gm_fod_norm="${dmri_dir}/gmfod_norm.mif"
csf_fod_norm="${dmri_dir}/csffod_norm.mif"

if [ ! -f ${wm_fod_norm} ]; then
    echo -e "${GREEN}[INFO]${NC} `date`: Running multi-tissue log-domain intensity normalisation"
    # Create eroded mask
    maskfilter -npass 2 "${dwi_meanbzero_brain_mask}" erode "${dmri_dir}/dwi_meanbzero_brain_mask_eroded_2.nii.gz" ${threading} -info
    # Perform normalization
    mtnormalise "${wm_fod}" "${wm_fod_norm}" "${gm_fod}" "${gm_fod_norm}" "${csf_fod}" \
                "${csf_fod_norm}" -mask "${dmri_dir}/dwi_meanbzero_brain_mask_eroded_2.nii.gz" ${threading} -info
fi

# ------------------------------------------------------------------------------
# 3. 5TT & Registration
# ------------------------------------------------------------------------------

# Create a mask of white matter gray matter interface (5TT)
freesurfer_5tt_T1="${dmri_dir}/5tt.T1.freesurfer.mif"
freesurfer_5tt="${dmri_dir}/5tt.freesurfer.mif"
T1_brain="${ukb_subjects_dir}/${ukb_subject_id}_${ukb_instance}/T1/T1_brain.nii.gz"
gmwm_seed_T1="${dmri_dir}/gmwm_seed_T1.mif"
gmwm_seed="${dmri_dir}/gmwm_seed.mif"
transform_DWI_T1_FSL="${dmri_dir}/diff2struct_fsl.txt"
transform_DWI_T1="${dmri_dir}/diff2struct_mrtrix.txt"

if [ ! -f ${gmwm_seed} ]; then
    echo -e "${GREEN}[INFO]${NC} `date`: Running 5ttgen to get GM/WM interface mask"
    
    # Generate 5TT image using FSL algorithm
    5ttgen fsl -premasked ${T1_brain} "${freesurfer_5tt_T1}" -nocrop -sgm_amyg_hipp ${threading} -info
    
    # Generate boundary ribbon
    5tt2gmwmi "${freesurfer_5tt_T1}" "${gmwm_seed_T1}" ${threading} -info

    # Rigid body registration (T1 -> DWI)
    flirt -in "${dwi_meanbzero_brain}" -ref "${T1_brain}" \
          -cost normmi -dof 6 -omat "${transform_DWI_T1_FSL}"
          
    transformconvert "${transform_DWI_T1_FSL}" "${dwi_meanbzero_brain}" \
                     "${T1_brain}" flirt_import "${transform_DWI_T1}"

    # Apply transformation
    mrtransform "${freesurfer_5tt_T1}" "${freesurfer_5tt}" -linear "${transform_DWI_T1}" -inverse ${threading} -info
    mrtransform "${gmwm_seed_T1}" "${gmwm_seed}" -linear "${transform_DWI_T1}" -inverse ${threading} -info
fi

# ------------------------------------------------------------------------------
# 4. Streamline Generation & SIFT2
# ------------------------------------------------------------------------------

# Generate Streamlines (tckgen)
tracks="${dmri_dir}/tracks_${streamlines}.tck"
trackstats="${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/tractography/stats/tracks_${streamlines}_stats.json"
mkdir -p "${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/tractography/stats"

if [ ! -f ${tracks} ]; then
    echo -e "${GREEN}[INFO]${NC} `date`: Running probabilistic tractography (iFOD2)"
    tckgen -seed_gmwmi "${gmwm_seed}" -act "${freesurfer_5tt}" -seeds "${streamlines}" \
           -maxlength 250 -cutoff 0.1 ${threading} "${wm_fod_norm}" "${tracks}" -power 0.5 \
           -info -samples 3 -select 0                
    tckstats "${tracks}" -output "${trackstats}"
fi

# Resample endpoints
endpoints="${dmri_dir}/tracks_${streamlines}_endpoints.tck"
if [ ! -f ${endpoints} ]; then
    echo -e "${GREEN}[INFO]${NC} `date`: Resampling streamline endpoints"
    tckresample ${threading} -info -endpoints "${tracks}" "${endpoints}"
fi

# SIFT2 Weighting
sift_weights_npy="${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/tractography/metrics/sift_weights.npy"
sift_stats="${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/tractography/stats/sift_stats.csv"
mkdir -p "${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/tractography/metrics"

if [ ! -f ${sift_weights_npy} ]; then
    echo -e "${GREEN}[INFO]${NC} `date`: Running SIFT2"
    tcksift2 ${threading} -info "${tracks}" -act "${freesurfer_5tt}" -config NPYFloatMaxSavePrecision 16 \
             -csv "${sift_stats}" "${wm_fod_norm}" "${sift_weights_npy}"
fi

# Sample metrics along streamlines (FA, Length, etc.)
streamline_length_npy="${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/tractography/metrics/streamline_metric_length.npy"
streamline_mean_fa_npy="${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/tractography/metrics/streamline_metric_FA_mean.npy"

if [ ! -f ${streamline_length_npy} ]; then
    echo -e "${GREEN}[INFO]${NC} `date`: Sampling metrics along tracks"
    # Length
    tckstats ${threading} -info -config NPYFloatMaxSavePrecision 16 -dump \
             "${streamline_length_npy}" "${tracks}"
    # FA
    tcksample ${threading} -precise -info -config NPYFloatMaxSavePrecision 16 -stat_tck mean \
              "${tracks}" "${dmri_dir}/dti_FA.nii.gz" "${streamline_mean_fa_npy}"
fi

# Convert endpoints to NPY binaries
endpoints_npy="${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/tractography/endpoints/tracks_${streamlines}_endpoints.npy"
mkdir -p "${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/tractography/endpoints"

if [ ! -f ${endpoints_npy} ]; then
    echo -e "${GREEN}[INFO]${NC} `date`: Converting endpoints to .npy binaries"
    python3 "${script_dir}/python/save_endpoints_as_npy.py" "${endpoints}" "${endpoints_npy}"
fi

echo -e "${GREEN}[INFO]${NC} `date`: Finished tractography for: ${ukb_subject_id}_${ukb_instance}"
