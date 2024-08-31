# Python script to extract the preprocessing and query times from the output files
# and calculate the mean times for each mesh and method of times respect to the 5 source vertices

import os
import re
import pandas as pd

# Define the path to the main folder (current directory since the script is inside the "distances" folder)
base_folder = "."

# Define a pattern to extract relevant data from the files, including vertices and faces
time_pattern = r"# Number of vertices: (\d+), Number of faces: (\d+)\n# Preprocessing time: ([\d\.e-]+)\n# Query time: ([\d\.e-]+)"

# Initialize a list to hold the results
results = []

# Iterate through each folder in the base folder
for folder in os.listdir(base_folder):
    folder_path = os.path.join(base_folder, folder)
    
    if os.path.isdir(folder_path):
        # Iterate through each file in the folder
        for filename in os.listdir(folder_path):
            file_path = os.path.join(folder_path, filename)
            
            if filename.endswith(".txt"):
                with open(file_path, 'r') as file:
                    content = file.read()
                    
                    # Extract the mesh name, method, and times
                    mesh_name = filename.split('_')[1]
                    method = filename.split('_')[2]
                    # print(f"Processing {filename}...")  # Debug: Show the file being processed
                    match = re.search(time_pattern, content)
                    
                    if match:
                        vertices = int(match.group(1))
                        faces = int(match.group(2))
                        preprocessing_time = float(match.group(3))
                        query_time = float(match.group(4))
                        
                        # Approximate very small values to 0
                        if preprocessing_time < 1e-06:
                            preprocessing_time = 0.0
                        if query_time < 1e-06:
                            query_time = 0.0
                        
                        results.append({
                            "Mesh": mesh_name,
                            "Method": method,
                            "Vertices": vertices,
                            "Faces": faces,
                            "Preprocessing Time": preprocessing_time,
                            "Query Time": query_time,
                            "Source Vertex": folder
                        })
                    else:
                        print(f"No match found for {filename}")  # Debug: If the pattern doesn't match

# Convert results to a DataFrame
df = pd.DataFrame(results)

# Write the detailed times data to a CSV file
df.to_csv("times.csv", index=False)

# Calculate the mean preprocessing and query times for each mesh and method
# Keep vertices and faces as additional columns by taking the first value (assuming they are the same for each mesh-method combination)
mean_times = df.groupby(['Mesh', 'Method', 'Vertices', 'Faces'], as_index=False).mean(numeric_only=True)

# Ensure that even zero values are included correctly in the output
mean_times.fillna(0, inplace=True)

# Save the resulting mean times DataFrame to another CSV file
output_csv = "mean_times.csv"
mean_times.to_csv(output_csv, index=False)

print(f"CSV file has been saved as {output_csv}")
