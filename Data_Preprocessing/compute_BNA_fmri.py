# script made from the notebook codes
import os
import sys
import numpy as np
import pandas as pd
import nibabel as nib

def ensure_dir(file_name):
    os.makedirs(os.path.dirname(file_name), exist_ok=True)
    return file_name

if __name__ == '__main__':
    # Parse command line arguments
    # Expected args: main_dir, ukb_subjects_dir, ukb_subject_id, ukb_instance, atlas_name
    if len(sys.argv) < 6:
        print("Usage: python compute_BNA_fmri.py <main_dir> <subjects_dir> <subject_id> <instance> <atlas_name>")
        sys.exit(1)

    main_dir, ukb_subjects_dir, ukb_subject_id, ukb_instance, atlas_name = sys.argv[1:]

    template_dir = "{}/data/templates".format(main_dir)
    temporary_dir = "{}/data/temporary".format(main_dir)
    output_dir = "{}/data/output".format(main_dir)

    # 1. Load the volumetric atlas (BNA, including cortical and subcortical regions)
    atlas_path = '{}/subjects/{}_{}/atlases/native.fMRI_space.{}.nii.gz'.format(temporary_dir, ukb_subject_id, ukb_instance, atlas_name)
    cortical_atlas = nib.load(atlas_path)

    # 2. Load names of all labels from the color lookup table
    # Note: Ensure this LUT file exists in your template directory
    lut_path = os.path.join(template_dir, 'BNA_Atlas', 'BN_Atlas_246_LUT.txt')
    
    cortical_labels = pd.DataFrame(
        np.genfromtxt(lut_path, dtype='str'),
        columns=['index', 'label_name', 'R', 'G', 'B', 'A'],
    ).astype(
        dtype={
            "index": "int",
            "label_name": "str",
            "R": "int",
            "G": "int",
            "B": "int",
            "A": "int",
        }
    )

    # 3. Load the ICA-cleaned fMRI data
    clean_fmri = nib.load('{}/{}_{}/fMRI/rfMRI.ica/filtered_func_data_clean.nii.gz'.format(ukb_subjects_dir, ukb_subject_id, ukb_instance))

    # 4. Extract label names
    cortical_atlas_fmri = cortical_labels[['index', 'label_name']][cortical_labels['label_name'] != '???'].copy()

    # 5. Calculate average fMRI signal over every label from the atlas
    # Resulting shape: (N_regions, N_timepoints)
    cortical_atlas_fmri = pd.concat(
        [
            cortical_atlas_fmri['label_name'],
            pd.DataFrame(
                np.array(
                    [
                        np.mean(clean_fmri.get_fdata()[cortical_atlas.get_fdata() == x['index']], axis=0)
                        for (i, x) in cortical_atlas_fmri.iterrows()
                    ]
                ),
                index=cortical_atlas_fmri.index,
                columns=['timepoint_{}'.format(x) for x in range(clean_fmri.shape[-1])],
            )
        ],
        axis=1
    )

    # 6. Write out the resulting time-series to a CSV file
    output_file = '{}/subjects/{}_{}/fMRI/fMRI.{}.csv.gz'.format(temporary_dir, ukb_subject_id, ukb_instance, atlas_name)
    cortical_atlas_fmri.to_csv(ensure_dir(output_file), index=False)
