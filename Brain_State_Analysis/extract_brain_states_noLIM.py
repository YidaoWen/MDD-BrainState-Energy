import warnings
warnings.filterwarnings("ignore")
import os
import numpy as np
import pandas as pd
import scipy.stats as stats

def tr_net_labels(df, nlabels, activity_thr):
    """
    Determine the dominant network for a single time point.
    
    Averages amplitude network-wise for each time point.
    If the highest amplitude >= activity_thr, returns the network ID (1-based).
    Otherwise, returns nlabels (which corresponds to NOTA).
    """
    net_means = []
    # Loop through networks (1 to nlabels-1). 
    # nlabels acts as the index for NOTA.
    for i in range(1, nlabels):
        val = df[df['network_id'] == i]['value'].mean()
        net_means.append(val)
    
    net_means = np.array(net_means)
    idx = np.argmax(net_means)
    
    # Check threshold logic (Standard for NoLIM analysis)
    if net_means[idx] >= activity_thr:
        return idx + 1
    else:
        return nlabels

# ==========================================
# Configuration
# ==========================================
datasets = ['HC', 'MDD']
activity_thr = 0

# Paths (Relative paths for portability)
ROOT_DIR = '../../ukb_raw_data/'
ATLAS_ORDER_FILE = '../resources/BNA_210_yeo_order.csv'
OUTPUT_ROOT = '../results/state_maps_results_NoLIM/'

# ==========================================
# Main Processing
# ==========================================

# 1. Setup Network Configuration (Removing Limbic)
original_labels = ['Vis', 'SomMot', 'DorsAttn', 'VentAttn', 'Limbic', 'Frontoparietal', 'Default', 'NOTA']

if not os.path.exists(ATLAS_ORDER_FILE):
    raise FileNotFoundError(f"Atlas order file not found: {ATLAS_ORDER_FILE}")
original_atlas_order = pd.read_csv(ATLAS_ORDER_FILE)

# Identify Limbic ID
limbic_label_name = 'Limbic'
try:
    limbic_id = original_labels.index(limbic_label_name) + 1
    print(f"Identified Limbic network original ID: {limbic_id}")
except ValueError:
    raise ValueError("Limbic label not found in configuration.")

# Filter Atlas: Remove rows belonging to Limbic
atlas_order = original_atlas_order[original_atlas_order['network_id'] != limbic_id].copy()

# Identify ROI indices to keep (0-based) from the data matrix
# Assumes 'ROI' column in CSV is 1-based
rois_to_keep_indices = atlas_order['ROI'].values - 1

# Update Labels list
labels = [label for label in original_labels if label != limbic_label_name]
nlabels = len(labels) # e.g., 7 (6 networks + NOTA)

print(f"Remaining networks: {labels[:-1]}")
print(f"Total labels (incl. NOTA): {nlabels}")

# Remap Network IDs (e.g., 1,2,4,5 -> 1,2,3,4)
unique_old_ids = sorted(atlas_order['network_id'].unique())
old_to_new_id_map = {old_id: new_id for new_id, old_id in enumerate(unique_old_ids, 1)}
atlas_order['network_id'] = atlas_order['network_id'].map(old_to_new_id_map)

nrois = len(atlas_order)
print(f"Remaining ROIs: {nrois}")

# Create Mappings
id2net = dict(zip(np.arange(nlabels) + 1, labels)) 
roi2id = dict(zip(atlas_order['ROI'], atlas_order['network_id'])) 


# 2. Process Subjects
for dataset in datasets:
    print(f"\nProcessing Group: {dataset}")
    
    dataset_dir = os.path.join(ROOT_DIR, f'rfMRI_bold_{dataset}')
    output_dir = os.path.join(OUTPUT_ROOT, dataset, f'thr_{activity_thr}')
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    subject_list_path = os.path.join(ROOT_DIR, f'{dataset}_list.txt')
    if not os.path.exists(subject_list_path):
        print(f"Subject list not found: {subject_list_path}")
        continue
    subjects = np.loadtxt(subject_list_path).astype(int)
    
    all_states_list = []

    for subj_id in subjects:
        # Load Time Series
        ttss_path = os.path.join(dataset_dir, f'{subj_id}_bold.csv')
        if not os.path.exists(ttss_path):
            continue
            
        ttss_246 = pd.read_csv(ttss_path, header=None).values
        
        # Cut to 210 and then filter for NoLIM ROIs
        ttss_full = ttss_246[:, :210]
        ttss = ttss_full[:, rois_to_keep_indices]

        # Z-score normalization
        ttss = stats.zscore(ttss, axis=0)
        
        # Flatten
        ttss_vec = ttss.T.flatten()
        n_timepoints = ttss.shape[0]
        curr_n_rois = ttss.shape[1]
        
        roi_order = np.repeat(np.arange(curr_n_rois), n_timepoints) + 1
        time_order = np.tile(np.arange(n_timepoints), curr_n_rois)

        ttss_df = pd.DataFrame({
            'value': pd.Series(ttss_vec, dtype='float'),
            'roi': pd.Series(roi_order, dtype='int32'),
            'tr': pd.Series(time_order, dtype='int32'),
        })
        
        # Correctly map continuous ROI index (1..N) back to Original ROI ID
        continuous_roi_to_original_roi = dict(enumerate(atlas_order['ROI'], 1))
        ttss_df['original_roi'] = ttss_df['roi'].map(continuous_roi_to_original_roi)
        
        # Map to (New) Network ID
        ttss_df['network_id'] = ttss_df['original_roi'].map(roi2id)
        # ttss_df['network'] = ttss_df['network_id'].map(id2net)

        # Calculate Labels
        time_labels = np.zeros(n_timepoints)
        for tr in range(n_timepoints):
            temp_label = tr_net_labels(ttss_df[ttss_df.tr==tr], nlabels, activity_thr)
            time_labels[tr] = temp_label

        np.save(os.path.join(output_dir, f'{subj_id}_dominance_network_labels.npy'), time_labels)
        
        # Generate State Maps (using new labels list)
        state_masks = pd.DataFrame({
            lab: (time_labels == lab_idx + 1) for lab_idx, lab in enumerate(labels)
        }).values

        for i, mask in enumerate(state_masks.T):
            if np.sum(mask) > 0:
                state_mean = np.mean(ttss[mask], axis=0)
            else:
                state_mean = np.zeros(curr_n_rois)
            
            if not np.all(np.isnan(state_mean)):
                tmp_df = atlas_order[['ROI', 'network_id']].copy()
                tmp_df['subject'] = f'{subj_id}'
                tmp_df['state'] = labels[i]
                tmp_df['value'] = state_mean
                all_states_list.append(tmp_df)
        
        # print(f"  Processed {subj_id}")

    # Save Group Results
    if all_states_list:
        group_states = pd.concat(all_states_list, ignore_index=True)
        group_states['ROI'] = group_states['ROI'].astype(int)
        group_states['network_id'] = group_states['network_id'].astype(int)
        
        save_path = os.path.join(output_dir, 'state-maps_subject-level.csv')
        group_states.to_csv(save_path, index=False)
        print(f"Saved NoLIM state maps to: {save_path}")
    else:
        print("No state data generated.")
