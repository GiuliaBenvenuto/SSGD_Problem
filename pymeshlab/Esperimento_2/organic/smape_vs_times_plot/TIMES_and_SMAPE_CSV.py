import os
import re
import pandas as pd
import numpy as np

def read_file_data(file_path):
    try:
        with open(file_path, 'r') as file:
            content = file.read()
        
        # Extract number of vertices and faces
        match = re.search(r"Number of vertices: (\d+), Number of faces: (\d+)", content)
        if not match:
            print(f"Warning: Couldn't extract vertices and faces from file {file_path}")
            return None
        vertices = int(match.group(1))
        faces = int(match.group(2))
        
        # Extract preprocessing and query times
        time_match = re.search(r"Preprocessing time: ([\d\.e-]+)\n# Query time: ([\d\.e-]+)", content)
        if not time_match:
            print(f"Warning: Couldn't extract times from file {file_path}")
            return None
        preprocessing_time = float(time_match.group(1))
        query_time = float(time_match.group(2))
        
        # Extract distances
        distances = list(map(float, re.findall(r"^[\d\.e-]+$", content, re.MULTILINE)))
        
        if not distances:
            print(f"Warning: No distances found in file {file_path}")
            return None

        return vertices, faces, preprocessing_time, query_time, distances
    except Exception as e:
        print(f"Error processing file {file_path}: {str(e)}")
        return None

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

# Define the path to the folder where the script should work
base_folder = "../DISTANCES_OK/93"  # Modify this to your desired path

results = []

# Iterate through each file in the base folder
for filename in os.listdir(base_folder):
    if filename.endswith(".txt"):
        parts = filename.split('_')
        if len(parts) < 3:
            print(f"Warning: Filename {filename} does not match expected format. Skipping.")
            continue
        
        # Extract mesh name and method
        mesh_name = f"{parts[0]}_{parts[1]}"
        method = parts[2].split('.')[0]
        
        # If the method is VTP, set it as ground truth
        if method == "VTP":
            gt_file_path = os.path.join(base_folder, filename)
            gt_data = read_file_data(gt_file_path)
            if gt_data is None:
                print(f"Failed to read ground truth file {filename}. Skipping.")
                continue
            
            gt_vertices, gt_faces, gt_preprocessing_time, gt_query_time, gt_distances = gt_data
            
            # Append the VTP method as its own line in the results
            results.append({
                "Mesh": mesh_name,
                "Method": "VTP",
                "Vertices": gt_vertices,
                "Faces": gt_faces,
                "Preprocessing Time": gt_preprocessing_time,
                "Query Time": gt_query_time,
                "SMAPE": 0.0  # SMAPE with itself is zero
            })
            
            # Now compare this ground truth with all other methods for the same mesh
            for compare_filename in os.listdir(base_folder):
                if compare_filename.startswith(mesh_name) and compare_filename != filename:
                    compare_file_path = os.path.join(base_folder, compare_filename)
                    compare_method = compare_filename.split('_')[2].split('.')[0]
                    
                    compare_data = read_file_data(compare_file_path)
                    if compare_data is None:
                        continue
                    
                    vertices, faces, preprocessing_time, query_time, distances = compare_data
                    
                    # Calculate SMAPE
                    smape = calculate_smape(gt_distances, distances)
                    
                    # Append results
                    results.append({
                        "Mesh": mesh_name,
                        "Method": compare_method,
                        "Vertices": vertices,
                        "Faces": faces,
                        "Preprocessing Time": preprocessing_time,
                        "Query Time": query_time,
                        "SMAPE": smape
                    })

print(f"Generated results for {len(results)} files")

# Convert results to a DataFrame
df = pd.DataFrame(results)

# Save the DataFrame with SMAPE values
output_csv = "smape_times_with_VTP.csv"
df.to_csv(output_csv, index=False)

print(f"CSV file has been saved as {output_csv} with {len(df)} rows.")
if len(df) == 0:
    print("Warning: The CSV file is empty. No valid data was processed.")
