import os
import numpy as np
import pandas as pd

def calculate_smape(gt, est):
    gt = np.array(gt)
    est = np.array(est)
    if len(gt) == 0 or len(est) == 0:
        return 0.0
    min_length = min(len(gt), len(est))
    gt = gt[:min_length]
    est = est[:min_length]
    denom = np.abs(gt) + np.abs(est)
    nonzero = denom != 0
    smape = np.abs(gt[nonzero] - est[nonzero]) / denom[nonzero]
    smape = np.mean(smape) * 100
    return round(smape, 5)

def read_distances(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    data = []
    vertices = 0
    for line in lines:
        if line.startswith('#'):
            if 'vertices' in line:
                vertices = int(line.split(',')[0].split(':')[-1].strip())
        else:
            try:
                data.append(float(line.strip()))
            except ValueError:
                pass
    return data, vertices

def process_directory(directory, ground_truth_file):
    files = [f for f in os.listdir(directory) if f.endswith('.txt')]
    ground_truth_data, _ = read_distances(os.path.join(directory, ground_truth_file))

    results = {}
    for file in files:
        mesh_name = file.split('_')[0] + '_' + file.split('_')[1]  # Assumes format '0_blub_Edge_663.txt'
        method_name = 'SMAPE_' + file.split('_')[2]
        data, num_vertices = read_distances(os.path.join(directory, file))
        smape = calculate_smape(ground_truth_data, data)
        
        if mesh_name not in results:
            results[mesh_name] = {'MeshName': mesh_name.split('_')[1], 'NumVertices': num_vertices}
        results[mesh_name][method_name] = smape

    df = pd.DataFrame(list(results.values()))
    df.to_csv(os.path.join(".", 'DRAGON_SMAPE_16.csv'), index=False)

# Specify the directory where your files are located and the ground truth file name
directory = '../DISTANCES_dragon_reordered'  # Adjust the path as needed
ground_truth_file = 'M7_reordered_VTP_16.txt'  # Adjust the file name as needed

process_directory(directory, ground_truth_file)
