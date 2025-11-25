import os
import numpy as np
import scipy.io

def calculate_dynamics_for_group(subject_list, data_dir, n_states, tr_duration_s, total_timepoints):
    """
    Calculate temporal dynamic metrics (FO, DT, AR, TP) for a group of subjects.
    
    Args:
        subject_list (list): List of subject IDs.
        data_dir (str): Directory containing the .npy label files.
        n_states (int): Number of states.
        tr_duration_s (float): TR duration in seconds.
        total_timepoints (int): Total number of time points in the scan.
        
    Returns:
        tuple: (fo_all, dt_all, ar_all, tp_all, tp_nopersist_all)
    """
    scan_length_min = (total_timepoints * tr_duration_s) / 60.0
    
    # Initialize lists to store results
    fo_list, dt_list, ar_list, tp_list, tp_nopersist_list = [], [], [], [], []

    for subject_id in subject_list:
        file_path = os.path.join(data_dir, f'{subject_id}_dominance_network_labels.npy')
        
        if not os.path.exists(file_path):
            print(f"Warning: File not found {file_path}, skipping subject {subject_id}.")
            continue
            
        # Load state time series (1-based labels)
        time_labels = np.load(file_path).astype(int)
        
        # --- 1. Fractional Occupancy (FO) ---
        fo = np.array([np.sum(time_labels == i) for i in range(1, n_states + 1)]) / total_timepoints
        fo_list.append(fo)
        
        # --- 2. Dwell Time (DT) & Appearance Rate (AR) ---
        dt_subject = np.zeros(n_states)
        ar_subject = np.zeros(n_states)
        
        # Find transition points
        change_points = np.where(np.diff(time_labels) != 0)[0] + 1
        block_starts = np.insert(change_points, 0, 0)
        block_ends = np.append(change_points, total_timepoints)
        
        # Calculate block lengths and corresponding states
        block_lengths = block_ends - block_starts
        block_states = time_labels[block_starts] - 1 # Convert to 0-based index
        
        for i in range(n_states):
            current_state_block_lengths = block_lengths[block_states == i]
            
            if len(current_state_block_lengths) > 0:
                # Dwell Time = Mean block length * TR (s)
                dt_subject[i] = np.mean(current_state_block_lengths) * tr_duration_s
                # Appearance Rate = Number of blocks / Total scan time (min)
                ar_subject[i] = len(current_state_block_lengths) / scan_length_min
            else:
                dt_subject[i] = 0
                ar_subject[i] = 0
        
        dt_list.append(dt_subject)
        ar_list.append(ar_subject)
        
        # --- 3. Transition Probability (TP) ---
        tp_matrix = np.zeros((n_states, n_states))
        states_0based = time_labels - 1 
        
        for t in range(total_timepoints - 1):
            from_state = states_0based[t]
            to_state = states_0based[t+1]
            # Ensure valid state range
            if 0 <= from_state < n_states and 0 <= to_state < n_states:
                tp_matrix[from_state, to_state] += 1
        
        # Calculate TP with self-transitions (TP2D)
        row_sums = tp_matrix.sum(axis=1, keepdims=True)
        tp2d = np.divide(tp_matrix, row_sums, out=np.zeros_like(tp_matrix), where=row_sums!=0)
        
        # Calculate TP without self-transitions (TP2D_NoPersist)
        tp_nopersist = tp_matrix.copy()
        np.fill_diagonal(tp_nopersist, 0)
        row_sums_nopersist = tp_nopersist.sum(axis=1, keepdims=True)
        tp2d_nopersist = np.divide(tp_nopersist, row_sums_nopersist, out=np.zeros_like(tp_nopersist), where=row_sums_nopersist!=0)
        
        tp_list.append(tp2d)
        tp_nopersist_list.append(tp2d_nopersist)
        
        # print(f"Processed subject: {subject_id}") # Uncomment for verbose output

    # Stack lists into arrays and Transpose to match (n_states, n_subjects) format
    return (np.stack(fo_list).T, np.stack(dt_list).T, np.stack(ar_list).T, 
            np.stack(tp_list), np.stack(tp_nopersist_list))


# --- Main Execution ---
if __name__ == "__main__":
    
    # Parameters
    TOTAL_TIMEPOINTS = 490
    TR_DURATION_S = 0.735
    
    # Paths (Modify as needed)
    ROOT_DIR = '../../ukb_raw_data/'
    # Directory where outputs from step 1 (extract_brain_states.py) are stored
    DATA_SOURCE_DIR = '../results/state_maps_results_NoLIM/' 
    OUTPUT_DIR = '../results/states_metrics_NoLIM/'
    
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    # Load subject lists
    mdd_list = np.loadtxt(os.path.join(ROOT_DIR, 'MDD_list.txt')).astype(int)
    hc_list = np.loadtxt(os.path.join(ROOT_DIR, 'HC_list.txt')).astype(int)

    # Analysis Loop
    for thre_type in ['thr_0', 'nonthr']:
        if thre_type == 'nonthr':
            n_states = 6
        elif thre_type == 'thr_0':
            n_states = 7 # Includes NOTA initially

        print(f"\n--- Processing Condition: {thre_type} ---")
        
        # Process MDD Group
        print("  Processing MDD group...")
        mdd_data_dir = os.path.join(DATA_SOURCE_DIR, 'MDD', thre_type)
        MDDfo, MDDdt, MDDar, MDD_TP2D, MDD_TP2D_NoPersist = calculate_dynamics_for_group(
            mdd_list, mdd_data_dir, n_states, TR_DURATION_S, TOTAL_TIMEPOINTS
        )

        # Process HC Group
        print("  Processing HC group...")
        hc_data_dir = os.path.join(DATA_SOURCE_DIR, 'HC', thre_type)
        HCfo, HCdt, HCar, HC_TP2D, HC_TP2D_NoPersist = calculate_dynamics_for_group(
            hc_list, hc_data_dir, n_states, TR_DURATION_S, TOTAL_TIMEPOINTS
        )

        # --- Consistency Check: Remove NOTA state for 'thr_0' ---
        if thre_type == 'thr_0':
            print("  Removing NOTA state data to maintain consistency (7 -> 6 states)...")
            # Slice to keep only the first 6 states
            MDDfo = MDDfo[:-1, :]
            MDDdt = MDDdt[:-1, :]
            MDDar = MDDar[:-1, :]
            HCfo = HCfo[:-1, :]
            HCdt = HCdt[:-1, :]
            HCar = HCar[:-1, :]

            # Slice TP matrices (n_subs, n_states, n_states) -> keep 6x6
            MDD_TP2D = MDD_TP2D[:, :-1, :-1]
            MDD_TP2D_NoPersist = MDD_TP2D_NoPersist[:, :-1, :-1]
            HC_TP2D = HC_TP2D[:, :-1, :-1]
            HC_TP2D_NoPersist = HC_TP2D_NoPersist[:, :-1, :-1]
        
        # --- Save Results to .mat file ---
        output_mat_path = os.path.join(OUTPUT_DIR, f'states_{thre_type}_fo_dt_ar_tp.mat')
        
        mat_dict = {
            'MDDfo': MDDfo, 'MDDdt': MDDdt, 'MDDar': MDDar,
            'HCfo': HCfo, 'HCdt': HCdt, 'HCar': HCar,
            'MDD_TP2D': MDD_TP2D, 'MDD_TP2D_NoPersist': MDD_TP2D_NoPersist,
            'HC_TP2D': HC_TP2D, 'HC_TP2D_NoPersist': HC_TP2D_NoPersist
        }
        
        scipy.io.savemat(output_mat_path, mat_dict)
        print(f"  Results saved to: {output_mat_path}")

    print("\nAnalysis completed successfully.")
