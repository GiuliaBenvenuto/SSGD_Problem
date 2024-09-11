import os
import numpy as np
import pandas as pd

def read_distances_and_vertices(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    metadata = lines[1]  # Get the second line that contains vertices information
    num_vertices = int(metadata.split(':')[1].split(',')[0].strip())  # Parse number of vertices
    distances = [float(line.strip()) for line in lines[4:]]  # Skip the first four lines with metadata
    return distances, num_vertices


def calculate_smape(gt, est):
    # Convert inputs to numpy arrays
    gt = np.array(gt)
    est = np.array(est)

    # If either array is empty, return 0.0
    if len(gt) == 0 or len(est) == 0:
        return 0.0

    # Truncate both arrays to the length of the shorter array
    min_length = min(len(gt), len(est))
    gt = gt[:min_length]
    est = est[:min_length]

    # Calculate SMAPE
    denom = (np.abs(gt) + np.abs(est)) / 2
    nonzero = denom != 0
    smape = np.abs(gt[nonzero] - est[nonzero]) / denom[nonzero]
    smape = np.mean(smape) * 100

    # Truncate SMAPE to five decimal places
    return round(smape, 5)


def process_folder(folder_path):
    results = {}
    files = os.listdir(folder_path)
    mesh_methods = {}
    vertices_count = {}

    # Organize files by mesh and method
    for file in files:
        if file.endswith('.txt'):
            parts = file.split('_')
            mesh_name = parts[1]  # Extract the mesh name from after the first underscore
            method = parts[2].split('.')[0]  # Extract method name, assuming method is the third segment in the filename
            if mesh_name not in mesh_methods:
                mesh_methods[mesh_name] = {}
            file_path = os.path.join(folder_path, file)
            distances, num_vertices = read_distances_and_vertices(file_path)
            mesh_methods[mesh_name][method] = distances
            vertices_count[mesh_name] = num_vertices  # Store the number of vertices for each mesh
    
    # Compute SMAPE for each mesh between VTP and other methods
    for mesh, methods in mesh_methods.items():
        if 'VTP' in methods:
            gt_distances = methods['VTP']
            row_result = {'MeshName': mesh, 'NumVertices': vertices_count[mesh]}
            for method, distances in methods.items():
                if method != 'VTP':
                    smape = calculate_smape(gt_distances, distances)
                    row_result[f'SMAPE_{method}'] = smape
            results[mesh] = row_result  # Store results for each mesh in a dictionary

    # Convert results dictionary to DataFrame and save to CSV
    df = pd.DataFrame(list(results.values()))
    df.to_csv(os.path.join("./GEO_EXT_FM_corretti", 'SMAPE_organic_2519.csv'), index=False)

# Replace 'path_to_folder' with the path to your folder containing the txt files
process_folder('../DISTANCES_OK/2519')
