import argparse
import os
import sys
import numpy as np
import pandas as pd
import scipy.io as sio

# Import nctpy functions
try:
    from nctpy.pipelines import ComputeControlEnergy
    from nctpy.utils import normalize_state
except ImportError:
    print("Error: nctpy toolbox not found. Please ensure it is installed.")
    sys.exit(1)

def assemble_control_tasks(states_x0, states_xf, control_set, n_nodes):
    """
    Assemble control tasks with a specific control input matrix B.
    """
    n_states = states_x0.shape[1]
    control_tasks = []

    # Constraints and Parameters
    trajectory_constraints = np.eye(n_nodes)
    rho = 1

    for initial_idx in range(n_states):
        initial_state = normalize_state(states_x0[:, initial_idx])
        for target_idx in range(n_states):
            target_state = normalize_state(states_xf[:, target_idx])

            control_task = {
                "x0": initial_state,
                "xf": target_state,
                "B": control_set,  # This is the customized B matrix
                "S": trajectory_constraints,
                "rho": rho
            }
            control_tasks.append(control_task)        
            
    return control_tasks   

def get_regional_control_energy(subject_id, condition_type, threshold, root_dir, output_root, states_root, n_roi=210):
    """
    Calculate regional control energy (perturbation analysis) for a single subject.
    Iterates through each ROI, increases its control weight by 1, and computes transition energy.
    
    Args:
        subject_id (int/str): Subject Identifier.
        condition_type (str): Group condition ('HC' or 'MDD').
        threshold (str): Thresholding type.
        root_dir (str): Path to raw data.
        output_root (str): Path to save results.
        states_root (str): Path to extracted brain states.
        n_roi (int): Number of ROIs (Nodes).
    """
    
    # 1. Setup Output Directory
    results_dir = os.path.join(output_root, 'energy', 'B_add1', condition_type, threshold)
    if not os.path.exists(results_dir):
        os.makedirs(results_dir, exist_ok=True)
    
    # 2. Load Brain States
    states_file = os.path.join(states_root, condition_type, threshold, 'state-maps_subject-level.csv')
    if not os.path.exists(states_file):
        print(f"Error: State file not found at {states_file}")
        return

    states_all_df = pd.read_csv(states_file) 
    
    # Filter for current subject
    states_all_df['subject'] = states_all_df['subject'].astype(str)
    subj_states_df = states_all_df[states_all_df['subject'] == str(subject_id)]
    
    if subj_states_df.empty:
        print(f"Warning: No states found for subject {subject_id}")
        return

    # Define Labels (Dynamically detected)
    # Assumes standard labels, but works if subset is present
    standard_labels = ['Vis', 'SomMot', 'DorsAttn', 'VentAttn', 'Limbic', 'Frontoparietal', 'Default']
    present_labels = [l for l in standard_labels if l in subj_states_df['state'].unique()]
    
    n_states = len(present_labels)
    states_matrix = np.zeros((n_roi, n_states))

    for idx, label in enumerate(present_labels):
        temp_state = subj_states_df[subj_states_df['state'] == label]['value'].values
        if len(temp_state) == n_roi:
            states_matrix[:, idx] = temp_state
        else:
             print(f"Error: Dimension mismatch. Expected {n_roi}, got {len(temp_state)}")
             return

    # 3. Load Structural Connectivity (A Matrix)
    adjacency_path = os.path.join(root_dir, f'sc_{condition_type}', f'{subject_id}_FA_sc.csv')
    if not os.path.exists(adjacency_path):
        print(f"Error: SC file not found at {adjacency_path}")
        return
        
    adjacency_temp = pd.read_csv(adjacency_path, header=None).values
    adjacency_temp = adjacency_temp[:n_roi, :n_roi]   
    np.fill_diagonal(adjacency_temp, 0)

    # 4. Regional Perturbation Loop
    # Store results: (n_ROI, n_states, n_states)
    optimal_control_energy_t1 = np.zeros((n_roi, n_states, n_states))
    # optimal_control_energy_t3 = np.zeros((n_roi, n_states, n_states)) # Uncomment if needed

    print(f'Starting Regional Analysis for Subject: {subject_id} ({condition_type})', flush=True)

    for k in range(n_roi):
        # Construct heterogeneous control set
        # B matrix is Identity, but the k-th diagonal element is increased to 2
        control_set = np.eye(n_roi) 
        control_set[k, k] = 2 

        if k % 20 == 0: # Log progress every 20 nodes
            print(f'  Processing ROI {k}/{n_roi}...', flush=True)                             
        
        control_tasks_temp = assemble_control_tasks(states_matrix, states_matrix, control_set, n_roi) 
            
        # Compute Control Energy (T=1)
        compute_energy_t1 = ComputeControlEnergy(
            A=adjacency_temp, control_tasks=control_tasks_temp, system="continuous", c=1, T=1
        )
        compute_energy_t1.run()        
        optimal_control_energy_t1[k, :, :] = np.reshape(compute_energy_t1.E, (n_states, n_states))
        del compute_energy_t1

        # Compute Control Energy (T=3) - Optional
        # compute_energy_t3 = ComputeControlEnergy(
        #     A=adjacency_temp, control_tasks=control_tasks_temp, system="continuous", c=1, T=3
        # )
        # compute_energy_t3.run()        
        # optimal_control_energy_t3[k, :, :] = np.reshape(compute_energy_t3.E, (n_states, n_states))
        # del compute_energy_t3
              
    # 5. Save Results
    print(f'Saving regional energy results for {subject_id}...')
    
    # Save as .mat file (compatible with MATLAB for downstream analysis if needed)
    sio.savemat(os.path.join(results_dir, f'{subject_id}_B_add1_OCE_individualSC_T1_{threshold}.mat'),
                {'oce': optimal_control_energy_t1})
                
    # sio.savemat(os.path.join(results_dir, f'{subject_id}_B_add1_OCE_individualSC_T3_{threshold}.mat'),
    #             {'oce': optimal_control_energy_t3})

# CLI Entry Point
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate Regional Control Energy (rERC analysis).")
    parser.add_argument("--threshold", default='nonthr', help="State extraction threshold type")
    parser.add_argument("--condition", required=True, choices=['HC', 'MDD'], help="Group Condition")
    parser.add_argument("--sub_id", required=True, help="Subject ID")
    parser.add_argument("--root_dir", default='../../ukb_raw_data/', help="Path to raw data")
    parser.add_argument("--output_dir", default='../results/TCE/', help="Path to save results")
    parser.add_argument("--states_dir", default='../results/state_maps_results/', help="Path to extracted states")
    
    args = parser.parse_args()
    print(f"Arguments: {args}")

    get_regional_control_energy(
        subject_id=args.sub_id,
        condition_type=args.condition,
        threshold=args.threshold,
        root_dir=args.root_dir,
        output_root=args.output_dir,
        states_root=args.states_dir
    )
