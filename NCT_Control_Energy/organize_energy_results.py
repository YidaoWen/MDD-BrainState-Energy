import argparse
import os
import numpy as np
import pandas as pd
import scipy.io as sio

def load_subject_data(root_dir, condition, threshold, sub_id, T_value):
    """
    Load energy matrix for a single subject from CSV or MAT.
    Assumes output from compute_control_energy.py
    """
    # Construct expected file path
    # e.g., ../results/energy/HC/nonthr/123456_OCE_T1_nonthr.csv
    file_name = f'{sub_id}_OCE_T{T_value}_{threshold}.csv'
    file_path = os.path.join(root_dir, 'energy', condition, threshold, file_name)
    
    if os.path.exists(file_path):
        # Read CSV (no header)
        return pd.read_csv(file_path, header=None).values
    else:
        print(f"Warning: File not found {file_path}")
        return None

def organize_energy_metrics(root_dir, output_dir, mdd_list, hc_list, threshold, T_value):
    """
    Aggregates energy matrices and computes summary metrics (TE, PE, Stability).
    """
    print(f"\nProcessing: Threshold={threshold}, T={T_value}")
    
    # Initialize lists to store row data
    data_all_te = [] # For 7x7 flat (49 rows per sub)
    data_ave_te = [] # For average TE (1 row per sub)
    data_ave_pe = [] # For average PE (1 row per sub)
    data_stability = [] # For global stability (1 row per sub)
    
    # Process both groups
    groups = [('MDD', mdd_list), ('HC', hc_list)]
    
    for condition, sub_list in groups:
        print(f"  Loading {condition} group ({len(sub_list)} subjects)...")
        
        for sub_id in sub_list:
            matrix = load_subject_data(root_dir, condition, threshold, sub_id, T_value)
            
            if matrix is None:
                continue
                
            n_states = matrix.shape[0]
            n_transitions = n_states * n_states
            
            # --- 1. Flattened TE (All 7x7 pairs) ---
            # Flatten column-wise (Fortran order) to match MATLAB behavior if needed, 
            # or row-wise (C order). Let's use Row-major (C) for standard Python/Pandas logic.
            # If you specifically need column-major (F), change to 'F'.
            flat_matrix = matrix.flatten(order='C') 
            
            for k in range(n_transitions):
                data_all_te.append({
                    'subject_id': sub_id,
                    'condition': condition,
                    'state_k': f'TE_{k+1}', # 1-based index
                    'TE': flat_matrix[k]
                })

            # --- 2. Average Transition Energy (Ave TE) ---
            # Off-diagonal elements
            mask_off_diag = ~np.eye(n_states, dtype=bool)
            ave_te = np.mean(matrix[mask_off_diag])
            
            data_ave_te.append({
                'subject_id': sub_id,
                'condition': condition,
                'state_k': 'ave_TE',
                'ave_TE': ave_te
            })
            
            # --- 3. Average Persistence Energy (Ave PE) ---
            # Diagonal elements
            diag_elements = np.diag(matrix)
            ave_pe = np.mean(diag_elements)
            
            data_ave_pe.append({
                'subject_id': sub_id,
                'condition': condition,
                'state_k': 'ave_PE',
                'ave_PE': ave_pe
            })
            
            # --- 4. Global Stability ---
            # Defined as 1 / log10(Ave_PE)
            # Handle potential log(0) or very small numbers
            if ave_pe > 0:
                stability = 1 / np.log10(ave_pe)
            else:
                stability = np.nan
                
            data_stability.append({
                'subject_id': sub_id,
                'condition': condition,
                'state_k': 'global_stability',
                'stability': stability
            })

    # --- Save DataFrames ---
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Full 7x7 TE
    df_all_te = pd.DataFrame(data_all_te)
    path_all_te = os.path.join(output_dir, f'df_oce_T{T_value}_{threshold}.csv')
    df_all_te.to_csv(path_all_te, index=False)
    print(f"    Saved Full TE to: {path_all_te}")

    # 2. Average TE
    df_ave_te = pd.DataFrame(data_ave_te)
    path_ave_te = os.path.join(output_dir, f'df_aveTE_T{T_value}_{threshold}.csv')
    df_ave_te.to_csv(path_ave_te, index=False)
    print(f"    Saved Ave TE to: {path_ave_te}")

    # 3. Average PE
    df_ave_pe = pd.DataFrame(data_ave_pe)
    path_ave_pe = os.path.join(output_dir, f'df_avePE_T{T_value}_{threshold}.csv')
    df_ave_pe.to_csv(path_ave_pe, index=False)
    print(f"    Saved Ave PE to: {path_ave_pe}")
    
    # 4. Global Stability
    df_stability = pd.DataFrame(data_stability)
    path_stability = os.path.join(output_dir, f'df_global_stability_T{T_value}_{threshold}.csv')
    df_stability.to_csv(path_stability, index=False)
    print(f"    Saved Stability to: {path_stability}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Organize and aggregate energy metrics for statistical analysis.")
    parser.add_argument("--root_dir", default='../results/', help="Root directory containing 'energy' folder")
    parser.add_argument("--list_dir", default='../../ukb_raw_data/', help="Directory containing subject list text files")
    parser.add_argument("--output_dir", default='../results/df_FoDtAr_TP_E/', help="Directory to save aggregated CSVs")
    
    args = parser.parse_args()
    
    # Load Subject Lists
    mdd_list = np.loadtxt(os.path.join(args.list_dir, 'MDD_list.txt')).astype(int)
    hc_list = np.loadtxt(os.path.join(args.list_dir, 'HC_list.txt')).astype(int)
    
    # Loop through conditions
    # thresholds = ['thr_0', 'nonthr']
    # T_values = [1, 3]
    
    # Example: Run for thr_0 and T=3 as in your snippet
    for threshold in ['thr_0']:
        for T_val in [3]:
            organize_energy_metrics(
                args.root_dir, 
                args.output_dir, 
                mdd_list, 
                hc_list, 
                threshold, 
                T_val
            )
