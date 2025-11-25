#!/bin/bash

# ==============================================================================
# UK Biobank Connectivity Mapping Pipeline
# ==============================================================================
# This script is adapted from the UKB-connectomics pipeline:
# https://github.com/sina-mansour/UKB-connectomics/blob/main/scripts/bash/UKB_connectivity_mapping_pipeline.sh
#
# Description:
# This pipeline processes a single UK Biobank subject to generate:
# 1. Individual volumetric atlases (warped from MNI space)
# 2. Regional functional time-series (mapped to the Brainnetome Atlas)
# 3. Structural connectivity matrices (via probabilistic tractography)
#
# Usage:
#   ./preprocessing_pipeline.sh <codes_dir> <subjects_dir> <subject_id_idx> <n_threads> <delete_download>
#
# Arguments:
#   <codes_dir>      : Main directory containing scripts and data templates.
#   <subjects_dir>   : Directory where subject data is stored.
#   <subject_id_idx> : Line number in the subject list file to process (for array jobs).
#   <n_threads>      : Number of threads for tractography (default: 0, uses all available).
#   <delete_download>: Delete raw downloaded data after processing? (yes/no, default: yes).
# ==============================================================================

# ------------------------------------------------------------------------------
# Configuration & Setup
# ------------------------------------------------------------------------------

# Input arguments
main_dir=$1
ukb_subjects_dir=$2
line_number=$3
n_threads=${4:-0}
delete_download=${5:-yes}

# Define internal paths
working_dir=`pwd`
script_dir="${main_dir}/scripts"
template_dir="${main_dir}/data/templates"
temporary_dir="${main_dir}/data/temporary"
output_dir="${main_dir}/data/output"

# Logging colors
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

# ------------------------------------------------------------------------------
# 1. Subject Identification & Data Preparation
# ------------------------------------------------------------------------------

# Read subject information from the combined bulk file.
# Note: The input file list should be prepared in advance in the temporary directory.
# Format expected: <subject_id> <instance_id>
# Example: subject_instance=(`sed "${line_number}q;d" "${temporary_dir}/bulk/subject_list.combined"`)

# [Placeholder] Ensure the correct subject list file is referenced here.
subject_list_file="${temporary_dir}/bulk/subject_list.combined" 

if [ ! -f "$subject_list_file" ]; then
    echo -e "${RED}[ERROR] Subject list file not found: ${subject_list_file}${NC}"
    exit 1
fi

subject_instance=(`sed "${line_number}q;d" "${subject_list_file}"`)
ukb_subject_id=${subject_instance[0]}
ukb_instance=${subject_instance[1]}

echo -e "${GREEN}[INFO] Processing Subject: ${ukb_subject_id} (Instance: ${ukb_instance})${NC}"

# Execute download script (copies data from central repository to working directory)
"${script_dir}/bash/download_subject_data.sh" "${main_dir}" "${ukb_subjects_dir}" "${ukb_subject_id}" "${ukb_instance}" "${working_dir}"

# ------------------------------------------------------------------------------
# 2. Check for Previous Computations (Resume capability)
# ------------------------------------------------------------------------------
# Unzip previously computed results if available to avoid redundant processing.

# Check Atlases
out_atlases_zip="${output_dir}/subjects/${ukb_subject_id}_${ukb_instance}/${ukb_subject_id}_atlases_${ukb_instance}.zip"
if [ -f "${out_atlases_zip}" ]; then
    unzip -n "${out_atlases_zip}" -d "${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/"
else
    echo -e "${GREEN}[INFO] No previous computations found for atlases.${NC}"
fi

# Check fMRI data
out_fmri_zip="${output_dir}/subjects/${ukb_subject_id}_${ukb_instance}/${ukb_subject_id}_fMRI_${ukb_instance}.zip"
if [ -f "${out_fmri_zip}" ]; then
    unzip -n "${out_fmri_zip}" -d "${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/"
else
    echo -e "${GREEN}[INFO] No previous computations found for fMRI.${NC}"
fi

# Check Tractography data
out_tractography_zip="${output_dir}/subjects/${ukb_subject_id}_${ukb_instance}/${ukb_subject_id}_tractography_${ukb_instance}.zip"
if [ -f "${out_tractography_zip}" ]; then
    unzip -n "${out_tractography_zip}" -d "${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/"
else
    echo -e "${GREEN}[INFO] No previous computations found for tractography.${NC}"
fi

cd "${working_dir}"

# ------------------------------------------------------------------------------
# 3. Generate Native Volumetric Atlases
# ------------------------------------------------------------------------------
# Maps the Brainnetome Atlas (BNA) from MNI space to individual native space.

echo -e "${GREEN}[INFO] Mapping MNI subcortical atlases to native volume.${NC}"

mkdir -p "${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/transforms"
mkdir -p "${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/atlases"

# Compute inverse warp (MNI -> T1 Native)
echo -e "${GREEN}[INFO] Computing the inverse warp.${NC}"
inverse_warp="${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/transforms/MNI_to_T1_inversed_warp_coef.nii.gz"
if [ ! -f ${inverse_warp} ]; then
    invwarp --ref="${ukb_subjects_dir}/${ukb_subject_id}_${ukb_instance}/T1/T1_brain.nii.gz" \
            --warp="${ukb_subjects_dir}/${ukb_subject_id}_${ukb_instance}/T1/transforms/T1_to_MNI_warp_coef.nii.gz" \
            --out="${inverse_warp}"
fi

echo -e "${GREEN}[INFO] Mapping BNA atlas to native volume.${NC}"

atlas_name="BNA"
# [CONFIGURATION REQUIRED] Path to the standard atlas template
atlas_path="${main_dir}/data/templates/BNA_Atlas" 
atlas_location="${atlas_path}/tpl-MNI152NLin6Asym_res-01_atlas-Brainnetome_desc-246Parcels_dseg.nii.gz"

# Apply warp to map atlas to native T1 space
native_atlas_location="${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/atlases/native.${atlas_name}.nii.gz"
if [ ! -f ${native_atlas_location} ]; then
    applywarp --ref="${ukb_subjects_dir}/${ukb_subject_id}_${ukb_instance}/T1/T1_brain.nii.gz" \
              --in="${atlas_location}" \
              --warp="${inverse_warp}" \
              --out="${native_atlas_location}" \
              --interp=nn
fi

# Resample atlas to fMRI resolution
native_atlas_location_fMRI="${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/atlases/native.fMRI_space.${atlas_name}.nii.gz"
if [ ! -f ${native_atlas_location_fMRI} ]; then
    mri_vol2vol --mov "${native_atlas_location}" \
                --targ "${ukb_subjects_dir}/${ukb_subject_id}_${ukb_instance}/fMRI/rfMRI.ica/mean_func.nii.gz" \
                --interp nearest \
                --regheader \
                --o "${native_atlas_location_fMRI}"
fi

echo -e "${GREEN}[INFO] Native volumetric atlases generated.${NC}"

# ------------------------------------------------------------------------------
# 4. Map Functional Time-series
# ------------------------------------------------------------------------------

# Extract fMRI time-series for BNA regions (Cortical + Subcortical)
echo -e "${GREEN}[INFO] Mapping fMRI on ${atlas_name} atlas.${NC}"

fMRI_data="${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/fMRI/fMRI.${atlas_name}.csv.gz"
if [ ! -f ${fMRI_data} ]; then
    python3 "${script_dir}/python/compute_BNA_fmri.py" "${main_dir}" "${ukb_subjects_dir}" "${ukb_subject_id}" "${ukb_instance}" "${atlas_name}"
fi

# Extract Global Signal
echo -e "${GREEN}[INFO] Mapping fMRI for global signal.${NC}"

global_signal_fmri="${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/fMRI/fMRI.global_signal.csv.gz"
if [ ! -f ${global_signal_fmri} ]; then
    python3 "${script_dir}/python/compute_global_fmri.py" "${main_dir}" "${ukb_subjects_dir}" "${ukb_subject_id}" "${ukb_instance}"
fi

echo -e "${GREEN}[INFO] All fMRI time-series generated.${NC}"

# ------------------------------------------------------------------------------
# 5. Map Structural Connectivity
# ------------------------------------------------------------------------------

# Ensure fsaverage templates are available
if [ ! -d "${template_dir}/freesurfer/fsaverage" ]; then
    mkdir -p "${template_dir}/freesurfer/"
    cp -r "${FREESURFER_HOME}/subjects/fsaverage" "${template_dir}/freesurfer/"
fi

# Run probabilistic tractography
echo -e "${GREEN}[INFO] Running probabilistic tractography pipeline.${NC}"

streamlines="10M" # Number of streamlines

"${script_dir}/bash/probabilistic_tractography_native_space.sh" "${main_dir}" "${ukb_subjects_dir}" "${ukb_subject_id}" "${ukb_instance}" "${streamlines}" "${n_threads}"

echo -e "${GREEN}[INFO] Mapping structural connectomes on ${atlas_name}.${NC}"

# Map connectivity to the BNA atlas
"${script_dir}/bash/map_structural_connectivity.sh" "${main_dir}" "${ukb_subjects_dir}" "${ukb_subject_id}" "${ukb_instance}" "${streamlines}" "${atlas_name}" "${n_threads}"

# ------------------------------------------------------------------------------
# 6. Cleanup and Archiving
# ------------------------------------------------------------------------------

# Remove raw downloaded data if requested (to save space)
if [ -d "${ukb_subjects_dir}/${ukb_subject_id}_${ukb_instance}" ] && [ ${delete_download} == "yes" ]; then
    echo -e "${GREEN}[INFO] Deleting downloaded files.${NC}"
    rm -r "${ukb_subjects_dir}/${ukb_subject_id}_${ukb_instance}/"
fi

# Compress outputs
echo -e "${GREEN}[INFO] Compressing computed files.${NC}"
mkdir -p "${output_dir}/subjects/${ukb_subject_id}_${ukb_instance}/"
cd "${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/"

# Archive Atlases and Transforms
zip -urvTm "${output_dir}/subjects/${ukb_subject_id}_${ukb_instance}/${ukb_subject_id}_atlases_${ukb_instance}.zip" \
         "./atlases" "./transforms"       
         
# Archive Functional MRI data
zip -urvTm "${output_dir}/subjects/${ukb_subject_id}_${ukb_instance}/${ukb_subject_id}_fMRI_${ukb_instance}.zip" \
         "./fMRI"
         
# Archive Tractography data
zip -urvTm "${output_dir}/subjects/${ukb_subject_id}_${ukb_instance}/${ukb_subject_id}_tractography_${ukb_instance}.zip" \
         "./tractography"

echo -e "${GREEN}[INFO] Pipeline completed successfully for subject ${ukb_subject_id}.${NC}"
