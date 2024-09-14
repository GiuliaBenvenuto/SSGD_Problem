# import os
# import re
# import pandas as pd
# import numpy as np

# # Define the path to the main folder (current directory since the script is inside the "distances" folder)
# base_folder = "."

# # Define a pattern to extract relevant data from the files, including vertices and faces
# time_pattern = r"# Number of vertices: (\d+), Number of faces: (\d+)\n# Preprocessing time: ([\d\.e-]+)\n# Query time: ([\d\.e-]+)"

# # Initialize a list to hold the results
# results = []

# # Iterate through each folder in the base folder
# for folder in os.listdir(base_folder):
#     folder_path = os.path.join(base_folder, folder)

#     if os.path.isdir(folder_path):
#         # Iterate through each file in the folder
#         for filename in os.listdir(folder_path):
#             file_path = os.path.join(folder_path, filename)

#             if filename.endswith(".txt"):
#                 with open(file_path, 'r') as file:
#                     content = file.read()

#                     # Extract the mesh name, method, and times
#                     mesh_name = filename.split('_')[1]
#                     method = filename.split('_')[2]
#                     match = re.search(time_pattern, content)

#                     if match:
#                         vertices = int(match.group(1))
#                         faces = int(match.group(2))
#                         preprocessing_time = float(match.group(3))
#                         query_time = float(match.group(4))

#                         # Approximate very small values to 0
#                         if preprocessing_time < 1e-06:
#                             preprocessing_time = 0.0
#                         if query_time < 1e-06:
#                             query_time = 0.0

#                         # Extract distances from file
#                         distances = list(map(float, re.findall(r"distance: ([\d\.e-]+)", content)))

#                         results.append({
#                             "Mesh": mesh_name,
#                             "Method": method,
#                             "Vertices": vertices,
#                             "Faces": faces,
#                             "Preprocessing Time": preprocessing_time,
#                             "Query Time": query_time,
#                             "Source Vertex": folder,
#                             "Distances": distances
#                         })
#                     else:
#                         print(f"No match found for {filename}")

# # Convert results to a DataFrame
# df = pd.DataFrame(results)

# # Write the detailed times data to a CSV file
# df.to_csv("times_mech_ORIGINAL.csv", index=False)

# # Calculate the mean preprocessing and query times for each mesh and method
# # Keep vertices and faces as additional columns by taking the first value (assuming they are the same for each mesh-method combination)
# mean_times = df.groupby(['Mesh', 'Method', 'Vertices', 'Faces'], as_index=False).mean(numeric_only=True)

# # Ensure that even zero values are included correctly in the output
# mean_times.fillna(0, inplace=True)

# # Save the resulting mean times DataFrame to another CSV file
# output_csv = "mean_times_mech_ORIGINAL.csv"
# mean_times.to_csv(output_csv, index=False)

# # Load ground truth distances from the specific file
# gt_file = os.path.join(base_folder, "663/5_blub_VTP_663.txt")
# with open(gt_file, 'r') as gt_file:
#     gt_content = gt_file.read()
#     gt_distances = list(map(float, re.findall(r"distance: ([\d\.e-]+)", gt_content)))

# # Function to calculate SMAPE
# def calculate_smape(gt, est):
#     # Convert inputs to numpy arrays
#     gt = np.array(gt)
#     est = np.array(est)

#     # If either array is empty, return 0.0
#     if len(gt) == 0 or len(est) == 0:
#         return 0.0

#     # Truncate both arrays to the length of the shorter array
#     min_length = min(len(gt), len(est))
#     gt = gt[:min_length]
#     est = est[:min_length]

#     # Calculate SMAPE
#     denom = (np.abs(gt) + np.abs(est)) / 2
#     nonzero = denom != 0
#     smape = np.abs(gt[nonzero] - est[nonzero]) / denom[nonzero]
#     smape = np.mean(smape) * 100

#     # Truncate SMAPE to five decimal places
#     return round(smape, 5)

# # Calculate SMAPE for each row in the DataFrame
# df['SMAPE'] = df['Distances'].apply(lambda est: calculate_smape(gt_distances, est))

# # Save the DataFrame with SMAPE values
# df.to_csv("times_mech_with_smape.csv", index=False)


# print(f"CSV file has been saved as {output_csv} with SMAPE calculations.")



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
    
    denom = (np.abs(gt) + np.abs(est))/2
    nonzero = denom != 0
    smape = np.abs(gt[nonzero] - est[nonzero]) / denom[nonzero]
    smape = np.mean(smape) * 100
    
    return round(smape, 5)

# Define the path to the main folder
base_folder = "."

results = []

# Iterate through each folder in the base folder
for folder in os.listdir(base_folder):
    folder_path = os.path.join(base_folder, folder)

    if os.path.isdir(folder_path):
        print(f"Processing folder: {folder}")
        
        # First, read the ground truth file
        gt_file_path = os.path.join(folder_path, "5_blub_VTP_663.txt")
        if not os.path.exists(gt_file_path):
            print(f"Ground truth file not found in folder {folder}. Skipping this folder.")
            continue
        
        gt_data = read_file_data(gt_file_path)
        if gt_data is None:
            print(f"Failed to read ground truth file in folder {folder}. Skipping this folder.")
            continue
        
        gt_vertices, gt_faces, gt_preprocessing_time, gt_query_time, gt_distances = gt_data
        
        # Now process all other .txt files in the folder
        for filename in os.listdir(folder_path):
            if filename.endswith(".txt") and filename != "5_blub_VTP_663.txt":
                file_path = os.path.join(folder_path, filename)
                print(f"Processing file: {filename}")
                
                # Extract mesh name and method from filename
                parts = filename.split('_')
                if len(parts) < 3:
                    print(f"Warning: Filename {filename} does not match expected format. Skipping.")
                    continue
                # mesh name schould be part 0 and part 1
                mesh_name = parts[0] + "_" + parts[1]
                # mesh_name = parts[1]
                method = parts[2].split('.')[0]
                
                # Read file data
                file_data = read_file_data(file_path)
                if file_data is None:
                    continue
                
                vertices, faces, preprocessing_time, query_time, distances = file_data
                
                # Calculate SMAPE
                smape = calculate_smape(gt_distances, distances)
                
                # Append results
                results.append({
                    "Folder": folder,
                    "Mesh": mesh_name,
                    "Method": method,
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
output_csv = "smape_times_comparison.csv"
df.to_csv(output_csv, index=False)

print(f"CSV file has been saved as {output_csv} with {len(df)} rows.")
if len(df) == 0:
    print("Warning: The CSV file is empty. No valid data was processed.")
