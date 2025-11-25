import argparse
import os
import sys
import numpy as np
import pandas as pd
import scipy.io as sio

# Import nctpy functions (Ensure nctpy is installed or in PYTHONPATH)
try:
    from nctpy.pipelines import ComputeControlEnergy
    from nctpy.utils import normalize_state
except ImportError:
    print("Error: nctpy toolbox not found. Please install it from https://github.com/LindenParkesLab/nctpy")
    sys.exit(1)

def assemble_control_tasks(states_x0, states_xf, control_set, n_nodes):
    """
    Assemble a list of control tasks for calculating transition energy between all pairs of states.
    
    Args:
        states_x0 (np.ndarray): Initial states matrix (Nodes x States).
        states_xf (np.ndarray): Target states matrix (Nodes x States).
        control_set (np.ndarray): Control input matrix B (Nodes x Nodes).
        n_nodes (int): Number of brain nodes.
        
    Returns:
        list: A list of dictionary objects, each defining a control task (x0 -> xf).
    """
    n_states = states_x0.shape[1]
    control_tasks = []

    # Define state trajectory constraints (constrain full state trajectory equally)
    trajectory_constraints = np.eye(n_nodes)

    # Define mixing parameter rho (weighting for energy cost vs. state deviation)
    rho = 1

    # Loop through all initial states
    for initial_idx in range(n_states):
        initial_state = normalize_state(states_x0[:, initial_idx])
        
        # Loop through all target states
        for target_idx in range(n_states):
            target_state = normalize_state(states_xf[:, target_idx])

            control_task = {
                "x0": initial_state,
                "xf": target_state,
                "B": control_set,
                "S": trajectory_constraints,
                "rho": rho
            }
            control_tasks.append(control_task)        
            
    return control_tasks   

def get_control_energy(subject_id, condition_type, threshold, root_dir, output_root, states_root, n_roi=210):
    """
    Calculate optimal control energy for a single subject across all state transitions.
    
    Args:
        subject_id (int/str): Subject Identifier.
        condition_type (str): Group condition ('HC' or 'MDD').
        threshold (str): Thresholding type (e.g., 'nonthr' or 'thr_0').
        root_dir (str): Path to raw data.
        output_root (str): Path to save energy results.
        states_root (str): Path to extracted brain states.
        n_roi (int): Number of ROIs to include (default 210 for BNA subcortical cut).
    """
    
    # 1. Setup Output Directory
    results_dir = os.path.join(output_root, 'energy', condition_type, threshold)
    if not os.path.exists(results_dir):
        os.makedirs(results_dir, exist_ok=True)
    
    # 2. Load Brain States
    # State file path logic matches the output from Brain_State_Analysis
    states_file = os.path.join(states_root, condition_type, threshold, 'state-maps_subject-level.csv')
    
    if not os.path.exists(states_file):
        print(f"Error: State file not found at {states_file}")
        return

    states_all_df = pd.read_csv(states_file) 
    
    # Filter for current subject
    # Check if 'subject' column exists, handle potential type mismatch (int vs str)
    states_all_df['subject'] = states_all_df['subject'].astype(str)
    subj_states_df = states_all_df[states_all_df['subject'] == str(subject_id)]
    
    if subj_states_df.empty:
        print(f"Warning: No states found for subject {subject_id}")
        return

    # Define Labels (Order matters!)
    # Note: If running NoLIM analysis, ensure the input states CSV does not contain Limbic, 
    # or handle the label list dynamically. Here we assume standard 7 networks for simplicity.
    labels = ['Vis', 'SomMot', 'DorsAttn', 'VentAttn', 'Limbic', 'Frontoparietal', 'Default']
    # Filter out labels not present in the dataframe (handles NoLIM case implicitly if label is missing)
    present_labels = [l for l in labels if l in subj_states_df['state'].unique()]
    
    n_states = len(present_labels)
    states_matrix = np.zeros((n_roi, n_states))

    for idx, label in enumerate(present_labels):
        temp_state = subj_states_df[subj_states_df['state'] == label]['value'].values
        # Ensure dimensionality matches
        if len(temp_state) == n_roi:
            states_matrix[:, idx] = temp_state
        else:
             print(f"Error: Dimension mismatch for subject {subject_id}, state {label}. Expected {n_roi}, got {len(temp_state)}")
             return

    # 3. Load Structural Connectivity (A Matrix)
    # Assumes file naming: {subject_id}_FA_sc.csv
    adjacency_path = os.path.join(root_dir, f'sc_{condition_type}', f'{subject_id}_FA_sc.csv')
    
    if not os.path.exists(adjacency_path):
        print(f"Error: SC file not found at {adjacency_path}")
        return
        
    adjacency_temp = pd.read_csv(adjacency_path, header=None).values
    # Slice to ROI count (e.g., 210)
    adjacency_temp = adjacency_temp[:n_roi, :n_roi]   
    np.fill_diagonal(adjacency_temp, 0) # Remove self-connections
            
    # 4. Prepare Control Tasks
    control_set = np.eye(n_roi) # Full control (B = Identity)
    control_tasks = assemble_control_tasks(states_matrix, states_matrix, control_set, n_roi)

    print(f'Processing Subject: {subject_id} | Condition: {condition_type} | States: {n_states}', flush=True)

    # 5. Compute Control Energy (T=1)
    # T represents the time horizon for the transition
    compute_energy_t1 = ComputeControlEnergy(
        A=adjacency_temp, control_tasks=control_tasks, system="continuous", c=1, T=1
    )
    compute_energy_t1.run()        
    oce_t1 = np.reshape(compute_energy_t1.E, (n_states, n_states))
    del compute_energy_t1

    # 6. Compute Control Energy (T=3) - Optional based on your research need
    # compute_energy_t3 = ComputeControlEnergy(
    #     A=adjacency_temp, control_tasks=control_tasks, system="continuous", c=1, T=3
    # )
    # compute_energy_t3.run()        
    # oce_t3 = np.reshape(compute_energy_t3.E, (n_states, n_states))
    # del compute_energy_t3
                
    # 7. Save Results
    print(f'Saving results for {subject_id}...')
    
    df_t1 = pd.DataFrame(oce_t1)
    df_t1.to_csv(os.path.join(results_dir, f'{subject_id}_OCE_T1_{threshold}.csv'), header=False, index=False)
    
    # df_t3 = pd.DataFrame(oce_t3)
    # df_t3.to_csv(os.path.join(results_dir, f'{subject_id}_OCE_T3_{threshold}.csv'), header=False, index=False)

if __name__ == "__main__":
    # Argument Parser for CLI usage (e.g., for cluster jobs)
    parser = argparse.ArgumentParser(description="Calculate Control Energy for a single subject.")
    parser.add_argument("--sub_id", required=True, help="Subject ID")
    parser.add_argument("--condition", required=True, choices=['HC', 'MDD'], help="Group Condition")
    parser.add_argument("--threshold", default='nonthr', help="State extraction threshold type")
    parser.add_argument("--root_dir", default='../../ukb_raw_data/', help="Path to raw data")
    parser.add_argument("--output_dir", default='../results/TCE/', help="Path to save results")
    parser.add_argument("--states_dir", default='../results/state_maps_results/', help="Path to extracted states")
    
    args = parser.parse_args()
    
    get_control_energy(
        subject_id=args.sub_id,
        condition_type=args.condition,
        threshold=args.threshold,
        root_dir=args.root_dir,
        output_root=args.output_dir,
        states_root=args.states_dir
    )
