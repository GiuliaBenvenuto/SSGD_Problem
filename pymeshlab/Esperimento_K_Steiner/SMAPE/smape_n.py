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
    denom = (np.abs(gt) + np.abs(est)) / 2
    nonzero = denom != 0
    smape = np.abs(gt[nonzero] - est[nonzero]) / denom[nonzero]
    smape = np.mean(smape) * 100
    return round(smape, 5)

def read_distances(file_path):
    with open(file_path, 'r') as file:
        data = file.readlines()
    distances = []
    for line in data:
        stripped_line = line.strip()
        if stripped_line and not stripped_line.startswith('#'):
            try:
                distance = float(stripped_line)
                distances.append(distance)
            except ValueError:
                continue  # Skip any lines that cannot be converted to float
    return distances

# Specify the directory containing the files
directory = '.'  # Update this to the correct directory path where your files are stored

# Read distances from the VTP file
vtp_file = [f for f in os.listdir(directory) if 'VTP' in f][0]
vtp_distances = read_distances(os.path.join(directory, vtp_file))

# Prepare to collect SMAPE results
results = []

# Iterate over files and calculate SMAPE
for file in os.listdir(directory):
    if 'VTP' not in file and file.endswith('.txt'):
        n_value = file.split('_')[3][1:]  # Extract N value, assuming it's formatted like N1, N3, etc.
        other_distances = read_distances(os.path.join(directory, file))
        smape = calculate_smape(vtp_distances, other_distances)
        results.append({'N': n_value, 'File': file.split('.')[0], 'SMAPE': smape})

# Convert results to DataFrame and save as CSV
results_df = pd.DataFrame(results)
results_df.to_csv(os.path.join(directory, 'N_SMAPE_Fertility.csv'), index=False)

print("Comparison complete. Results saved to N_SMAPE_Fertility.csv.")
