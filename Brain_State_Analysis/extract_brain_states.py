import warnings
warnings.filterwarnings("ignore")
import os
import glob
import numpy as np
import pandas as pd
import scipy.stats as stats

def tr_net_labels(df, nlabels, activity_thr):
    """
    Determine the dominant network for a single time point.
    
    Averages amplitude network-wise for each time point and assigns the label 
    of the network with the highest amplitude. If no network surpasses the 
    activity_thr, assigns the extra label NOTA ('None Of The Above').
    
    Returns:
        int: Index of the dominant network (1-based). Returns nlabels (NOTA) if threshold not met.
    """
    net_means = []
    # Loop through networks (1 to nlabels-1, excluding NOTA)
    for i in range(1, nlabels):
        mean_val = df[df['network_id'] == i]['value'].mean()
        net_means.append(mean_val)
    
    net_means = np.array(net_means)
    idx = np.argmax(net_means)
    
    if net_means[idx] >= activity_thr:
        return idx + 1
    else:
        return nlabels

# ------- Configuration -------
# Define dataset groups
datasets = ['HC', 'MDD']

# Configuration paths (Modify these for your environment)
# Use relative paths or environment variables for portability
ROOT_DIR = '../../ukb_raw_data/' 
OUTPUT_ROOT = '../results/state_maps_results/'
ATLAS_ORDER_FILE = '../resources/BNA_210_yeo_order.csv'

# Parameters
ACTIVITY_THR = 0 # Threshold for activity

# Define network labels
# NOTA is included for noise handling (thresholding)
LABELS = ['Vis', 'SomMot', 'DorsAttn', 'VentAttn', 'Limbic', 'Frontoparietal', 'Default', 'NOTA']

# Load Atlas Mapping
# The atlas_order file should contain columns: 'ROI' (1-based index) and 'network_id'
if not os.path.exists(ATLAS_ORDER_FILE):
    raise FileNotFoundError(f"Atlas order file not found: {ATLAS_ORDER_FILE}")

atlas_order = pd.read_csv(ATLAS_ORDER_FILE)
nlabels = len(LABELS)
nrois = len(atlas_order) 

# Mappings
id2net = dict(zip(np.arange(nlabels) + 1, LABELS))
roi2id = dict(zip(atlas_order['ROI'], atlas_order['network_id']))

# ------- Main Processing Loop -------

for dataset in datasets:
    print(f"Processing dataset: {dataset}")
    
    dataset_dir = os.path.join(ROOT_DIR, f'rfMRI_bold_{dataset}')
    output_dir = os.path.join(OUTPUT_ROOT, dataset, f'thr_{ACTIVITY_THR}')
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    # Load subject list
    subject_list_path = os.path.join(ROOT_DIR, f'{dataset}_list.txt')
    subjects = np.loadtxt(subject_list_path).astype(int)
    
    # DataFrame to store group-level state maps
    group_states = pd.DataFrame(columns=['subject', 'state', 'value'])

    for subj_id in subjects:
        print(f'  Starting subject {subj_id}')
        
        # Load time series
        ttss_path = os.path.join(dataset_dir, f'{subj_id}_bold.csv')
        if not os.path.exists(ttss_path):
            print(f"    Warning: File not found {ttss_path}, skipping.")
            continue
            
        ttss_246 = pd.read_csv(ttss_path, header=None).values # Shape: (Time, ROIs)
        
        # Use first 210 ROIs (excluding subcortical regions if consistent with BNA_210_yeo_order)
        ttss = ttss_246[:, :210] 
        
        # Z-score normalization along time axis
        ttss = stats.zscore(ttss, axis=0)

        # Prepare dataframe for processing
        # Flatten array: [Time, ROI] -> [ROI, Time] (transposed) -> Flattened
        ttss_vec = ttss.T.flatten() 
        n_timepoints = ttss.shape[0]
        n_rois_curr = ttss.shape[1]
        
        roi_order = np.repeat(np.arange(n_rois_curr), n_timepoints) + 1 
        time_order = np.tile(np.arange(n_timepoints), n_rois_curr)

        ttss_df = pd.DataFrame({
            'value': pd.Series(ttss_vec, dtype='float'),
            'roi': pd.Series(roi_order, dtype='int32'),
            'tr': pd.Series(time_order, dtype='int32'),
        })
        
        # Map ROIs to Networks
        ttss_df['network_id'] = ttss_df['roi'].map(roi2id)
        # ttss_df['network'] = ttss_df['network_id'].map(id2net) # Optional for debug

        # Calculate dominant network for each TR
        time_labels = np.zeros(n_timepoints)
        for tr in range(n_timepoints):
            temp_label = tr_net_labels(ttss_df[ttss_df.tr == tr], nlabels, ACTIVITY_THR)
            time_labels[tr] = temp_label

        # Save dominance labels (The state sequence)
        np.save(os.path.join(output_dir, f'{subj_id}_dominance_network_labels.npy'), time_labels)
        
        # --- Create State Maps (Mean activity during dominance) ---
        state_masks = pd.DataFrame({
            lab: (time_labels == lab_idx + 1) for lab_idx, lab in enumerate(LABELS)
        }).values

        for i, mask in enumerate(state_masks.T):
            # Calculate mean activity pattern for this state
            if np.sum(mask) > 0:
                state_map = np.mean(ttss[mask], axis=0)
            else:
                state_map = np.zeros(n_rois_curr) # Handle case where state never appears

            tmp_df = atlas_order[['ROI', 'network_id']].copy()
            tmp_df['subject'] = f'{subj_id}'
            tmp_df['state'] = LABELS[i]
            tmp_df['value'] = state_map
            
            # Note: In newer pandas versions, append is deprecated. Using concat is preferred.
            group_states = pd.concat([group_states, tmp_df], ignore_index=True)

    # Clean up types
    group_states['ROI'] = group_states['ROI'].fillna(0).astype(int)
    group_states['network_id'] = group_states['network_id'].fillna(0).astype(int)

    # Save group-level state maps
    group_states.to_csv(os.path.join(output_dir, 'state-maps_subject-level.csv'), index=False)
    print(f"Finished processing {dataset}. Results saved.")
