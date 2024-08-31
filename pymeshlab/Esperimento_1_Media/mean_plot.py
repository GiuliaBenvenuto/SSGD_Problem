# import os
# import numpy as np
# import matplotlib.pyplot as plt

# def read_distances_from_file(file_path):
#     distances = []
#     num_vertices = 0
#     num_faces = 0
#     with open(file_path, 'r') as file:
#         for line in file:
#             if line.startswith('#'):
#                 if "Number of vertices" in line and "Number of faces" in line:
#                     parts = line.strip().split(',')
#                     num_vertices = int(parts[0].split(':')[1].strip())
#                     num_faces = int(parts[1].split(':')[1].strip())
#             else:
#                 distances.append(float(line.strip()))
#     return np.array(distances), num_vertices, num_faces

# def process_directory(directory_path, key_string, reference_filename):
#     file_data = {}
#     mesh_resolutions = {}
#     reference_included = False
    
#     for filename in sorted(os.listdir(directory_path)):
#         if filename.endswith('.txt') and (key_string in filename or filename == reference_filename):
#             file_path = os.path.join(directory_path, filename)
#             distances, num_vertices, num_faces = read_distances_from_file(file_path)
#             file_data[filename] = distances
#             mesh_resolutions[filename] = num_faces  # Use num_faces as resolution
#             if filename == reference_filename:
#                 reference_included = True
    
#     # Ensure the reference file is always included if not already
#     if not reference_included:
#         reference_file_path = os.path.join(directory_path, reference_filename)
#         if os.path.exists(reference_file_path):
#             distances, num_vertices, num_faces = read_distances_from_file(reference_file_path)
#             file_data[reference_filename] = distances
#             mesh_resolutions[reference_filename] = num_faces
#         else:
#             raise FileNotFoundError(f"Reference file {reference_filename} not found in the directory.")
    
#     return file_data, mesh_resolutions

# import numpy as np

# def compute_mean_percentage_errors(data, reference_key):
#     if reference_key not in data:
#         raise ValueError(f"Reference file {reference_key} not found in the processed data.")
    
#     mean_percentage_errors = {}
#     reference_data = data[reference_key]
    
#     for key, values in data.items():
#         if key == reference_key:
#             continue  # Skip the reference file
        
#         min_length = min(len(values), len(reference_data))
#         sliced_values = values[:min_length]
#         sliced_reference_data = reference_data[:min_length]
        
#         # Calculate percentage error using the modified SMAPE formula
#         numerator = sliced_values - sliced_reference_data
#         denominator = (np.abs(sliced_values) + np.abs(sliced_reference_data)) / 2
        
#         with np.errstate(divide='ignore', invalid='ignore'):
#             percentage_error = 100 * numerator / denominator
        
#         # Filter out invalid values (such as division by zero or NaN)
#         valid_mask = np.isfinite(percentage_error)
#         valid_percentage_error = percentage_error[valid_mask]
        
#         if len(valid_percentage_error) > 0:
#             # Calculate mean of the valid percentage errors
#             mean_error = np.mean(np.abs(valid_percentage_error))
#             mean_percentage_errors[key] = mean_error
    
#     return mean_percentage_errors


# def plot_mean_errors(mean_errors, mesh_resolutions):
#     fig, ax = plt.subplots(figsize=(10, 6))
#     methods = set(k.split('_')[0] for k in mean_errors.keys())
    
#     for method in methods:
#         x_values = []
#         y_values = []
#         for key, mean_error in mean_errors.items():
#             if key.startswith(method):
#                 x_values.append(mesh_resolutions[key])
#                 y_values.append(mean_error)
        
#         # Sort the data points by x_values (resolution)
#         sorted_indices = np.argsort(x_values)
#         x_values = np.array(x_values)[sorted_indices]
#         y_values = np.array(y_values)[sorted_indices]

#         # Debugging output to check data
#         print(f"Method: {method}, x_values: {x_values}, y_values: {y_values}")

#         # Plot with lines connecting the dots
#         ax.plot(x_values, y_values, marker='o', label=method, linestyle='-', markersize=5)

#         # Annotate the mean error next to the data points
#         for x, y in zip(x_values, y_values):
#             ax.annotate(f"{y:.2f}", xy=(x, y), xytext=(5, 5), textcoords="offset points", ha='center')
    
#     ax.set_title('Mean Error vs Mesh Resolution')
#     ax.set_xlabel('Number of Faces')
#     ax.set_ylabel('Mean Error')
#     ax.set_xscale('log')
#     ax.legend(title="Methods")
#     ax.grid(True)
    
#     plt.show()


# # Directory and reference setup
# directory_path = 'data/distances_dragon_subdiv'
# key_string = 'Edge'
# reference_key = 'DS7_VTP_16.txt'

# # Process directory, ensuring reference file is included
# file_data, mesh_resolutions = process_directory(directory_path, key_string, reference_key)

# # Compute mean errors and plot them
# if reference_key in file_data:
#     mean_errors = compute_mean_percentage_errors(file_data, reference_key)
#     plot_mean_errors(mean_errors, mesh_resolutions)
# else:
#     print(f"Reference file {reference_key} was not found.")


import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

def read_distances_from_file(file_path):
    distances = []
    num_vertices = 0
    num_faces = 0
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                if "Number of vertices" in line and "Number of faces" in line:
                    parts = line.strip().split(',')
                    num_vertices = int(parts[0].split(':')[1].strip())
                    num_faces = int(parts[1].split(':')[1].strip())
            else:
                distances.append(float(line.strip()))
    return np.array(distances), num_vertices, num_faces

def process_directory(directory_path, key_string, reference_filename):
    file_data = {}
    mesh_resolutions = {}
    reference_included = False
    
    for filename in sorted(os.listdir(directory_path)):
        if filename.endswith('.txt') and (key_string in filename or filename == reference_filename):
            file_path = os.path.join(directory_path, filename)
            distances, num_vertices, num_faces = read_distances_from_file(file_path)
            file_data[filename] = distances
            mesh_resolutions[filename] = num_faces  # Use num_faces as resolution
            if filename == reference_filename:
                reference_included = True
    
    # Ensure the reference file is always included if not already
    if not reference_included:
        reference_file_path = os.path.join(directory_path, reference_filename)
        if os.path.exists(reference_file_path):
            distances, num_vertices, num_faces = read_distances_from_file(reference_file_path)
            file_data[reference_filename] = distances
            mesh_resolutions[reference_filename] = num_faces
        else:
            raise FileNotFoundError(f"Reference file {reference_filename} not found in the directory.")
    
    return file_data, mesh_resolutions

import numpy as np

def compute_mean_percentage_errors(data, reference_key):
    if reference_key not in data:
        raise ValueError(f"Reference file {reference_key} not found in the processed data.")
    
    mean_percentage_errors = {}
    reference_data = data[reference_key]
    
    for key, values in data.items():
        if key == reference_key:
            continue  # Skip the reference file
        
        min_length = min(len(values), len(reference_data))
        sliced_values = values[:min_length]
        sliced_reference_data = reference_data[:min_length]
        
        # Calculate percentage error using the modified SMAPE formula
        numerator = sliced_values - sliced_reference_data
        denominator = (np.abs(sliced_values) + np.abs(sliced_reference_data)) / 2
        
        with np.errstate(divide='ignore', invalid='ignore'):
            percentage_error = 100 * numerator / denominator
        
        # Filter out invalid values (such as division by zero or NaN)
        valid_mask = np.isfinite(percentage_error)
        valid_percentage_error = percentage_error[valid_mask]
        
        if len(valid_percentage_error) > 0:
            # Calculate mean of the valid percentage errors
            mean_error = np.mean(np.abs(valid_percentage_error))
            mean_percentage_errors[key] = mean_error
    
    return mean_percentage_errors


# def plot_mean_errors(mean_errors, mesh_resolutions):
#     # Prepare the data in a format suitable for Seaborn
#     data = []
#     methods = set(k.split('_')[0] for k in mean_errors.keys())
    
#     for method in methods:
#         for key, mean_error in mean_errors.items():
#             if key.startswith(method):
#                 data.append({
#                     'Method': method,
#                     'Resolution': mesh_resolutions[key],
#                     'Mean Error': mean_error
#                 })
    
#     # Convert the data to a pandas DataFrame
#     df = pd.DataFrame(data)

#     # Debugging output to check data
#     print(df)

#     # Set up the plot
#     plt.figure(figsize=(10, 6))
    
#     # Create a line plot
#     sns.lineplot(data=df, x='Resolution', y='Mean Error', hue='Method', marker='o')

#     # Annotate the mean error next to the data points
#     for method in methods:
#         subset = df[df['Method'] == method]
#         for i in range(subset.shape[0]):
#             plt.text(
#                 subset['Resolution'].iloc[i], subset['Mean Error'].iloc[i],
#                 f"{subset['Mean Error'].iloc[i]:.2f}", 
#                 ha='center', va='bottom', fontsize=9, color='black'
#             )

#     # Set the scale for x-axis to logarithmic
#     plt.xscale('log')
    
#     # Set title and labels
#     plt.title('Mean Error vs Mesh Resolution')
#     plt.xlabel('Number of Faces')
#     plt.ylabel('Mean Error')
    
#     # Display the legend and grid
#     plt.legend(title='Methods')
#     plt.grid(True)
    
#     # Show the plot
#     plt.show()

def plot_mean_errors(mean_errors, mesh_resolutions):
    data = [{'Method': key, 'Resolution': mesh_resolutions[key], 'Mean Error': mean_errors[key]} for key in mean_errors]
    df = pd.DataFrame(data)
    df = df.sort_values(by='Resolution')

    plt.figure(figsize=(15, 8))
    plot = sns.lineplot(data=df, x='Resolution', y='Mean Error', hue='Method', marker='o', markersize=10)
    plt.xscale('log')
    
    # Adjust the title and position
    plt.title('Mean Error vs Mesh Resolution', y=1.05, fontsize=16, fontweight='bold')
    
    plt.xlabel('Number of Faces')
    plt.ylabel('Mean Error')
    plt.grid(True)

    # Enhancements for annotation
    for i, row in df.iterrows():
        x, y = row['Resolution'], row['Mean Error']
        annotation_text = f"$\\mathbf{{Mean\ Error:}}$ {y:.3f}\n$\\mathbf{{Resolution:}}$ {x}"

        vertical_alignment = 'bottom' if i % 2 == 0 else 'top'
        horizontal_alignment = 'right' if x > df['Resolution'].median() else 'left'
        offset = (5, 10) if vertical_alignment == 'bottom' else (-5, -10)

        plt.annotate(annotation_text, (x, y), textcoords="offset points", xytext=offset, ha=horizontal_alignment,
                     va=vertical_alignment, fontsize=10, color='black',  # Increased font size to 10
                     arrowprops=dict(arrowstyle='->', color='gray'),
                     bbox=dict(boxstyle="round,pad=0.3", facecolor='white', edgecolor='black', alpha=1.0))

    plt.tight_layout()  # Adjust layout to make room for legend and title
    plt.show()



# Directory and reference setup
directory_path = 'data/distances_dragon_subdiv'
key_string = 'Edge'
reference_key = 'DS7_VTP_16.txt'

# Process directory, ensuring reference file is included
file_data, mesh_resolutions = process_directory(directory_path, key_string, reference_key)

# Compute mean errors and plot them
if reference_key in file_data:
    mean_errors = compute_mean_percentage_errors(file_data, reference_key)
    plot_mean_errors(mean_errors, mesh_resolutions)
else:
    print(f"Reference file {reference_key} was not found.")