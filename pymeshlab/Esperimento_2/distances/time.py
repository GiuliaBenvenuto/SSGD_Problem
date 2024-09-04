import os
import numpy as np
import pandas as pd

def read_distances(file_path):
    with open(file_path, 'r') as file:
        data = file.readlines()[4:]  # Skip the first four lines with metadata
    distances = [float(line.strip()) for line in data]
    return distances

def calculate_smape(gt, est):
    gt = np.array(gt)
    est = np.array(est)
    if len(gt) == 0 or len(est) == 0:
        return 0.0
    denom = np.abs(gt) + np.abs(est)
    nonzero = denom != 0
    smape = np.abs(gt[nonzero] - est[nonzero]) / denom[nonzero]
    return np.mean(smape) * 100  # Convert to percentage

def process_folder(folder_path):
    results = []
    files = os.listdir(folder_path)
    mesh_methods = {}

    # Organize files by mesh and method
    for file in files:
        if file.endswith('.txt'):
            parts = file.split('_')
            mesh_name = parts[1]  # Corrected to extract mesh name from the second part of the filename
            method = parts[2].split('.')[0]  # Extract method from the third part of the filename
            if mesh_name not in mesh_methods:
                mesh_methods[mesh_name] = {}
            mesh_methods[mesh_name][method] = os.path.join(folder_path, file)
    
    # Compute SMAPE for each mesh between VTP and other methods
    for mesh, methods in mesh_methods.items():
        if 'VTP' in methods:
            gt_distances = read_distances(methods['VTP'])
            row_result = {'MeshName': mesh}
            for method, path in methods.items():
                if method != 'VTP':
                    est_distances = read_distances(path)
                    smape = calculate_smape(gt_distances, est_distances)
                    row_result[f'SMAPE_{method}'] = smape
            results.append(row_result)
    
    # Create DataFrame and save to CSV
    df = pd.DataFrame(results)
    df.to_csv(os.path.join(folder_path, 'smape_results.csv'), index=False)

# Replace 'path_to_folder' with the path to your folder containing the txt files
process_folder('path_to_folder')
