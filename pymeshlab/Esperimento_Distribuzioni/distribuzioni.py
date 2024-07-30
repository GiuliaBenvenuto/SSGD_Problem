# import os
# import numpy as np

# def read_distances_from_file(file_path):
#     distances = []
#     with open(file_path, 'r') as file:
#         for line in file:
#             if not line.startswith('#'):
#                 distances.append(float(line.strip()))
#     return np.array(distances)

# def process_directory(directory_path, key_string):
#     file_data = {}
    
#     for filename in sorted(os.listdir(directory_path)):
#         if filename.endswith('.txt') and key_string in filename:
#             file_path = os.path.join(directory_path, filename)
#             distances = read_distances_from_file(file_path)
#             file_data[filename] = distances
    
#     return file_data

# def compute_differences(data, reference_key):
#     differences = {}
#     reference_data = data[reference_key]
    
#     for key, values in data.items():
#         # Determine the minimum length to use for slicing
#         min_length = min(len(values), len(reference_data))
#         # Slice both arrays to the minimum length before subtracting
#         sliced_values = values[:min_length]
#         sliced_reference_data = reference_data[:min_length]
#         differences[key] = sliced_values - sliced_reference_data
    
#     return differences

# # Replace this with the path to your directory
# directory_path = 'data/distances_bunny_500f'

# # Process only files containing 'VTP'
# file_data = process_directory(directory_path, 'VTP')

# # Assuming the last VTP file is your reference ('Mk')
# reference_key = 'bunny_500f_6_VTP_100.txt'

# # Compute differences
# differences = compute_differences(file_data, reference_key)

# # Output differences
# for filename, diff in differences.items():
#     print(f"Differences for {filename}:")
#     print(diff)
#     print()






# import os
# import numpy as np
# import matplotlib.pyplot as plt

# def read_distances_from_file(file_path):
#     distances = []
#     with open(file_path, 'r') as file:
#         for line in file:
#             if not line.startswith('#'):
#                 distances.append(float(line.strip()))
#     return np.array(distances)

# def process_directory(directory_path, key_string, reference_filename):
#     file_data = {}
#     reference_included = False
    
#     for filename in sorted(os.listdir(directory_path)):
#         if filename.endswith('.txt') and (key_string in filename or filename == reference_filename):
#             file_path = os.path.join(directory_path, filename)
#             distances = read_distances_from_file(file_path)
#             file_data[filename] = distances
#             if filename == reference_filename:
#                 reference_included = True
    
#     # Ensure the reference file is always included if not already
#     if not reference_included:
#         reference_file_path = os.path.join(directory_path, reference_filename)
#         if os.path.exists(reference_file_path):
#             distances = read_distances_from_file(reference_file_path)
#             file_data[reference_filename] = distances
#         else:
#             raise FileNotFoundError(f"Reference file {reference_filename} not found in the directory.")
    
#     return file_data

# def compute_differences(data, reference_key):
#     if reference_key not in data:
#         raise ValueError(f"Reference file {reference_key} not found in the processed data.")
#     differences = {}
#     reference_data = data[reference_key]
#     for key, values in data.items():
#         min_length = min(len(values), len(reference_data))
#         sliced_values = values[:min_length]
#         sliced_reference_data = reference_data[:min_length]
#         differences[key] = sliced_values - sliced_reference_data
#     return differences

# def plot_differences(differences, reference_key):
#     fig, ax = plt.subplots(figsize=(10, 6))
#     global_min = min(np.min(diff) for diff in differences.values())
#     global_max = max(np.max(diff) for diff in differences.values())
#     x_range = np.linspace(global_min, global_max, 300)
#     for key, diff in differences.items():
#         if key == reference_key:
#             ax.axvline(0, color='red', linewidth=2, linestyle='-', label=f"{key} (reference, all zeros)")
#         else:
#             mean, std = np.mean(diff), np.std(diff)
#             p = 1 / (std * np.sqrt(2 * np.pi)) * np.exp(-((x_range - mean) ** 2) / (2 * std**2))
#             ax.plot(x_range, p, linewidth=2, label=key)
#     ax.axhline(0, color='black', linewidth=0.5, linestyle='--')
#     ax.axvline(0, color='black', linewidth=0.5, linestyle='--')
#     ax.set_title('Distribution of Differences')
#     ax.set_xlabel('Difference')
#     ax.set_ylabel('Density')
#     ax.legend()
#     plt.show()

# # Directory and reference setup
# directory_path = 'data/distances_bunny_500f'
# key_string = 'Extended'  # Filter criteria for non-reference files
# reference_key = 'bunny_500f_6_VTP_100.txt'  # Reference file

# # Process directory, ensuring reference file is included
# file_data = process_directory(directory_path, key_string, reference_key)

# # Compute and plot differences
# if reference_key in file_data:
#     differences = compute_differences(file_data, reference_key)
#     for filename, diff in differences.items():
#         print(f"Differences for {filename}:")
#         print(diff)
#         print()
#     plot_differences(differences, reference_key)
# else:
#     print(f"Reference file {reference_key} was not found.")

# plot_differences(differences, reference_key)



import os
import numpy as np
import matplotlib.pyplot as plt

def read_distances_from_file(file_path):
    distances = []
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                distances.append(float(line.strip()))
    return np.array(distances)

def process_directory(directory_path, key_string, reference_filename):
    file_data = {}
    reference_included = False
    
    for filename in sorted(os.listdir(directory_path)):
        if filename.endswith('.txt') and (key_string in filename or filename == reference_filename):
            file_path = os.path.join(directory_path, filename)
            distances = read_distances_from_file(file_path)
            file_data[filename] = distances
            if filename == reference_filename:
                reference_included = True
    
    # Ensure the reference file is always included if not already
    if not reference_included:
        reference_file_path = os.path.join(directory_path, reference_filename)
        if os.path.exists(reference_file_path):
            distances = read_distances_from_file(reference_file_path)
            file_data[reference_filename] = distances
        else:
            raise FileNotFoundError(f"Reference file {reference_filename} not found in the directory.")
    
    return file_data

def compute_differences(data, reference_key):
    if reference_key not in data:
        raise ValueError(f"Reference file {reference_key} not found in the processed data.")
    differences = {}
    reference_data = data[reference_key]
    for key, values in data.items():
        min_length = min(len(values), len(reference_data))
        sliced_values = values[:min_length]
        sliced_reference_data = reference_data[:min_length]
        differences[key] = sliced_values - sliced_reference_data
    return differences


fig, ax = None, None
def plot_differences(differences, reference_key):
    global fig, ax  # Refer to the global variables inside the function
    fig, ax = plt.subplots(figsize=(10, 6))
    lines = []
    global_min = min(np.min(diff) for diff in differences.values())
    global_max = max(np.max(diff) for diff in differences.values())
    x_range = np.linspace(global_min, global_max, 300)

    for key, diff in differences.items():
        mean, std = np.mean(diff), np.std(diff)
        p = 1 / (std * np.sqrt(2 * np.pi)) * np.exp(-((x_range - mean) ** 2) / (2 * std**2))
        line, = ax.plot(x_range, p, linewidth=2, label=key)
        lines.append(line)

    ax.axhline(0, color='black', linewidth=0.5, linestyle='--')
    ax.axvline(0, color='black', linewidth=0.5, linestyle='--')
    ax.set_title('Distribution of Differences')
    ax.set_xlabel('Difference')
    ax.set_ylabel('Density')
    leg = ax.legend(fancybox=True, shadow=True)
    lined = {}
    for legline, origline in zip(leg.get_lines(), lines):
        legline.set_picker(5)
        lined[legline] = origline

    fig.canvas.mpl_connect('pick_event', lambda event: toggle_visibility(event, fig, lined))
    plt.show()

def toggle_visibility(event, fig, lined):
    legline = event.artist
    origline = lined[legline]
    visible = not origline.get_visible()
    origline.set_visible(visible)
    legline.set_alpha(1.0 if visible else 0.2)
    fig.canvas.draw()

# Directory and reference setup
directory_path = 'data/distances_bunny_500f'
key_string = 'VTP'  # Filter criteria for non-reference files
reference_key = 'bunny_500f_6_VTP_100.txt'  # Reference file

# Process directory, ensuring reference file is included
file_data = process_directory(directory_path, key_string, reference_key)

# Compute and plot differences
if reference_key in file_data:
    differences = compute_differences(file_data, reference_key)
    for filename, diff in differences.items():
        print(f"Differences for {filename}:")
        print(diff)
        print()
    plot_differences(differences, reference_key)
else:
    print(f"Reference file {reference_key} was not found.")
