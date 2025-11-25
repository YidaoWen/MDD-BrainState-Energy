#!/bin/bash

# ==============================================================================
# Structural Connectivity Mapping Script
# ==============================================================================
# This script maps the generated streamlines to the BNA atlas to construct connectomes.

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

main_dir=$1
ukb_subjects_dir=$2
ukb_subject_id=$3
ukb_instance=$4
streamlines=$5
cortical_atlas_name=$6 # Represents the unified BNA atlas
num_threads=$7

# Ensure MRtrix bin is in PATH
MRTRIX_CMD="" 

script_dir="${main_dir}/scripts"
template_dir="${main_dir}/data/templates"
temporary_dir="${main_dir}/data/temporary"

dmri_dir="${ukb_subjects_dir}/${ukb_subject_id}_${ukb_instance}/dMRI/dMRI"

threading="-nthreads ${num_threads}" 

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

# Define atlas path
cortical_atlas_file="${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/atlases/native.${cortical_atlas_name}.nii.gz"

echo -e "${GREEN}[INFO]${NC} `date`: Starting connectivity mapping for: ${ukb_subject_id}_${ukb_instance} on ${cortical_atlas_name}."

# ------------------------------------------------------------------------------
# 1. Transform Atlas to DWI Space
# ------------------------------------------------------------------------------

cortical_atlas_dwi="${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/tractography/atlases/native.dMRI_space.${cortical_atlas_name}.nii.gz"
transform_DWI_T1="${dmri_dir}/diff2struct_mrtrix.txt"
mkdir -p "${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/tractography/atlases"

if [ ! -f ${cortical_atlas_dwi} ]; then
    echo -e "${GREEN}[INFO]${NC} `date`: Transforming atlas to dMRI space"
    # Using nearest neighbor interpolation for label preservation
    mrtransform "${cortical_atlas_file}" "${cortical_atlas_dwi}" -linear "${transform_DWI_T1}" -inverse -interp nearest \
                -datatype uint32 ${threading} -info
fi

# Copy to combinations folder (maintaining structure for compatibility)
combined_atlas_dwi="${dmri_dir}/atlases/combinations/native.dMRI_space.${cortical_atlas_name}.nii.gz"
mkdir -p "${dmri_dir}/atlases/combinations"
if [ ! -f ${combined_atlas_dwi} ]; then
    cp "${cortical_atlas_dwi}" "${combined_atlas_dwi}"
fi

# ------------------------------------------------------------------------------
# 2. Compute Connectomes
# ------------------------------------------------------------------------------

# Define inputs
endpoints="${dmri_dir}/tracks_${streamlines}_endpoints.tck"
sift_weights="${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/tractography/metrics/sift_weights.npy"
streamline_length="${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/tractography/metrics/streamline_metric_length.npy"
streamline_mean_fa="${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/tractography/metrics/streamline_metric_FA_mean.npy"

# Define outputs
streamline_count="${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/tractography/connectomes/${cortical_atlas_name}/connectome_streamline_count_${streamlines}.csv"
sift2_fbc="${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/tractography/connectomes/${cortical_atlas_name}/connectome_sift2_fbc_${streamlines}.csv"
mean_length="${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/tractography/connectomes/${cortical_atlas_name}/connectome_mean_length_${streamlines}.csv"
mean_fa="${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/tractography/connectomes/${cortical_atlas_name}/connectome_mean_FA_${streamlines}.csv"

mkdir -p "${temporary_dir}/subjects/${ukb_subject_id}_${ukb_instance}/tractography/connectomes/${cortical_atlas_name}/"

if [ ! -f ${streamline_count} ]; then
    
    # 1. Streamline Count (NOS)
    echo -e "${GREEN}[INFO]${NC} `date`: Computing connectomes - Streamline Count"
    tck2connectome ${threading} -info -symmetric -assignment_radial_search 4 \
                   "${endpoints}" "${combined_atlas_dwi}" "${streamline_count}"
                   
    # 2. SIFT2 Fiber Bundle Capacity (FBC)
    echo -e "${GREEN}[INFO]${NC} `date`: Computing connectomes - SIFT2 FBC"
    tck2connectome ${threading} -info -symmetric -assignment_radial_search 4 -tck_weights_in \
                   "${sift_weights}" "${endpoints}" "${combined_atlas_dwi}" "${sift2_fbc}"
                   
    # 3. Mean Length
    echo -e "${GREEN}[INFO]${NC} `date`: Computing connectomes - Mean Length"
    tck2connectome ${threading} -info -symmetric -assignment_radial_search 4 -tck_weights_in \
                   "${sift_weights}" -scale_file "${streamline_length}" -stat_edge mean \
                   "${endpoints}" "${combined_atlas_dwi}" "${mean_length}"
                   
    # 4. Mean Fractional Anisotropy (FA)
    echo -e "${GREEN}[INFO]${NC} `date`: Computing connectomes - Mean FA"
    tck2connectome ${threading} -info -symmetric -assignment_radial_search 4 -tck_weights_in \
                   "${sift_weights}" -scale_file "${streamline_mean_fa}" -stat_edge mean \
                   "${endpoints}" "${combined_atlas_dwi}" "${mean_fa}"
fi

echo -e "${GREEN}[INFO]${NC} `date`: Finished connectivity mapping for: ${ukb_subject_id}_${ukb_instance}."
