# --------- VERSIONE CON ERRORE ASSOLUTO ------------
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


# fig, ax = None, None
# def plot_differences(differences, reference_key):
#     global fig, ax  # Refer to the global variables inside the function
#     fig, ax = plt.subplots(figsize=(10, 6))
#     lines = []
#     global_min = min(np.min(diff) for diff in differences.values())
#     global_max = max(np.max(diff) for diff in differences.values())
#     x_range = np.linspace(global_min, global_max, 300)

#     for key, diff in differences.items():
#         mean, std = np.mean(diff), np.std(diff)
#         p = 1 / (std * np.sqrt(2 * np.pi)) * np.exp(-((x_range - mean) ** 2) / (2 * std**2))
#         line, = ax.plot(x_range, p, linewidth=2, label=key)
#         lines.append(line)

#     ax.axhline(0, color='black', linewidth=0.5, linestyle='--')
#     ax.axvline(0, color='black', linewidth=0.5, linestyle='--')
#     ax.set_title('Distribution of Differences')
#     ax.set_xlabel('Difference')
#     ax.set_ylabel('Density')
#     leg = ax.legend(fancybox=True, shadow=True)
#     lined = {}
#     for legline, origline in zip(leg.get_lines(), lines):
#         legline.set_picker(5)
#         lined[legline] = origline

#     fig.canvas.mpl_connect('pick_event', lambda event: toggle_visibility(event, fig, lined))
#     plt.show()

# def toggle_visibility(event, fig, lined):
#     legline = event.artist
#     origline = lined[legline]
#     visible = not origline.get_visible()
#     origline.set_visible(visible)
#     legline.set_alpha(1.0 if visible else 0.2)
#     fig.canvas.draw()

# # Directory and reference setup
# directory_path = 'data/distances_bob_500f'
# key_string = 'Extended'  # Filter criteria for non-reference files
# reference_key = 'bob_500f_6_VTP_100.txt'  # Reference file

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




# ----------- VERSIONE CON ERRORE PERCENTUALE ------------
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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


def compute_percentage_errors(data, reference_key):
    if reference_key not in data:
        raise ValueError(f"Reference file {reference_key} not found in the processed data.")
    percentage_errors = {}
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
        
        # Filter out invalid values
        valid_mask = np.isfinite(percentage_error)
        percentage_errors[key] = percentage_error[valid_mask]
        
        print(f"File: {key}")
        print(f"Valid errors: {np.sum(valid_mask)} out of {len(percentage_error)}")
        print(f"Invalid errors: {np.sum(~valid_mask)}")
        print(f"Mean: {np.mean(percentage_errors[key]):.2f}%, Std: {np.std(percentage_errors[key]):.2f}%")
        print()
    
    return percentage_errors


# # Gaussian distributions of the percentage errors
def plot_percentage_errors(percentage_errors, reference_key):
    global fig, ax
    fig, ax = plt.subplots(figsize=(14, 8))
    lines = []
    x_range = np.linspace(-100, 100, 300)  # Error range from -100% to 100%


    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']

    for i, (key, errors) in enumerate(percentage_errors.items()):
        if key != reference_key and len(errors) > 0:
            mean, std = np.mean(errors), np.std(errors)
            # p: probability density function (PDF) of a normal (Gaussian) distribution.
            # It's calculated using the formula for a normal distribution: p(x) = (1 / (σ * √(2π))) * e^(-(x - μ)^2 / (2σ^2)) 
            # where μ is the mean and σ is the standard deviation.
            # This formula creates a bell-shaped curve centered at the mean, with its width determined by the standard deviation.
            p = 1 / (std * np.sqrt(2 * np.pi)) * np.exp(-((x_range - mean) ** 2) / (2 * std**2))

            line, = ax.plot(x_range, p, linewidth=2, label=f"{key} (μ={mean:.2f}, σ={std:.2f})", color=colors[i % len(colors)])
            lines.append(line)

    if not lines:
        print("No valid data to plot.")
        return

    # Draw the x and y axis lines prominently
    ax.axhline(0, color='black', linewidth=1)
    ax.axvline(0, color='black', linewidth=1)

    ax.set_title('Dragon Mesh Set Created with Loop Subdivision\n'
             'Distribution of Percentage Errors for "Dragon" meshes with VTP Method',
              fontweight='bold')
    ax.set_xlabel('Percentage Error going from -100% to 100%\n (Zoom between -10% and 10%)')
    ax.set_ylabel('Density')
    ax.set_xlim(-10, 10)

    # Add a grey grid
    ax.grid(True, which='both', linestyle='--', linewidth=0.5, color='grey')

    leg = ax.legend(fancybox=True, shadow=True)
    lined = {}
    for legline, origline in zip(leg.get_lines(), lines):
        legline.set_picker(5)
        lined[legline] = origline

    fig.canvas.mpl_connect('pick_event', lambda event: toggle_visibility(event, fig, lined))
    plt.show()


# WITH SEABORN
# def plot_percentage_errors(percentage_errors, reference_key):
#     # Set the style and color palette
#     sns.set_style("whitegrid")
#     sns.set_palette("husl", 8)
    
#     fig, ax = plt.subplots(figsize=(16, 10))
    
#     x_range = np.linspace(-100, 100, 300)
    
#     for key, errors in percentage_errors.items():
#         if key != reference_key and len(errors) > 0:
#             mean, std = np.mean(errors), np.std(errors)
#             p = 1 / (std * np.sqrt(2 * np.pi)) * np.exp(-((x_range - mean) ** 2) / (2 * std**2))
            
#             sns.lineplot(x=x_range, y=p, label=f"{key} (μ={mean:.2f}, σ={std:.2f})", linewidth=2.5)
    
#     # Customize the plot
#     ax.set_title('Dragon Mesh Set Created with Loop Subdivision', fontsize=20, fontweight='bold')
#     ax.set_title('Distribution of Percentage Errors for "Dragon" meshes with Edge Method', fontsize=20, fontweight='bold')
#     ax.set_xlabel('Percentage Error going from -100% to +100%\n (Zoom between -100% and +100%)', fontsize=16)
#     # ax.set_xlabel('Percentage Error going from -100% to +100%', fontsize=16)

#     ax.set_ylabel('Density', fontsize=16)
#     ax.set_xlim(-100, 100)
    
#     # Add a shaded area for better readability
#     ax.axvspan(-10, 10, alpha=0.1, color='gray')
    
#     # Customize the legend
#     leg = ax.legend(title="Mesh Comparisons", fontsize=12, title_fontsize=14)
#     leg.get_frame().set_alpha(0.8)
    
#     # Improve tick labels
#     ax.tick_params(axis='both', which='major', labelsize=12)

#     # Add dotted grid lines
#     ax.grid(True, linestyle=':', alpha=0.7)
    
#     # Adjust layout and display
#     plt.tight_layout()
#     plt.show()



def toggle_visibility(event, fig, lined):
    legline = event.artist
    origline = lined[legline]
    visible = not origline.get_visible()
    origline.set_visible(visible)
    legline.set_alpha(1.0 if visible else 0.2)
    fig.canvas.draw()


# Directory and reference setup
# directory_path = 'data/distances_bob_500f'
# directory_path = 'data/distances_bob_500f'
# directory_path = 'data/distances_asian_dragon'
directory_path = 'data/distances_dragon_subdiv'

key_string = 'Trettner'  # Filter criteria for non-reference files
# reference_key = 'bob_500f_6_VTP_100.txt'  # Reference file - GROUND TRUTH
reference_key = 'DS7_VTP_16.txt'


# Process directory, ensuring reference file is included
file_data = process_directory(directory_path, key_string, reference_key)


# Compute and plot percentage errors
if reference_key in file_data:
    percentage_errors = compute_percentage_errors(file_data, reference_key)
    for filename, errors in percentage_errors.items():
        if filename != reference_key:
            print(f"Percentage Errors for {filename}:")
            print(errors)
            print(f"Mean: {np.mean(errors):.2f}%, Std: {np.std(errors):.2f}%")
            print()
    plot_percentage_errors(percentage_errors, reference_key)
else:
    print(f"Reference file {reference_key} was not found.")