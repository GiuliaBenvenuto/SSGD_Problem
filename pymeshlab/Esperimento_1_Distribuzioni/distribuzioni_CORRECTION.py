# ----------- VERSIONE CON ERRORE PERCENTUALE ------------
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm 

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


## QUELLO CON CUI HO FATTO GLI ULTIMI GRAFICI
# Gaussian distributions of the percentage errors
# def plot_percentage_errors(percentage_errors, reference_key):
#     global fig, ax
#     fig, ax = plt.subplots(figsize=(14, 8))
#     lines = []
    
#     x_range = np.linspace(-100, 100, 100000)  # Error range from -100% to 100%

#     colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']

#     for i, (key, errors) in enumerate(percentage_errors.items()):
#         if key != reference_key and len(errors) > 0:
#             mean, std = np.mean(errors), np.std(errors)
#             # p: probability density function (PDF) of a normal (Gaussian) distribution.
#             # It's calculated using the formula for a normal distribution: p(x) = (1 / (σ * √(2π))) * e^(-(x - μ)^2 / (2σ^2)) 
#             # where μ is the mean and σ is the standard deviation.
#             # This formula creates a bell-shaped curve centered at the mean, with its width determined by the standard deviation.
#             p = 1 / (std * np.sqrt(2 * np.pi)) * np.exp(-((x_range - mean) ** 2) / (2 * std**2))

#             line, = ax.plot(x_range, p, linewidth=2, label=f"{key} (μ={mean:.2f}, σ={std:.2f})", color=colors[i % len(colors)])
#             lines.append(line)

#     if not lines:
#         print("No valid data to plot.")
#         return

#     # Draw the x and y axis lines prominently
#     ax.axhline(0, color='black', linewidth=1)
#     ax.axvline(0, color='black', linewidth=1)

#     ax.set_title(f'Dragon Mesh Set Created with Decimation\n'
#              f'Distribution of Percentage Errors for "Dragon" meshes with {key_string} Method',
#              fontweight='bold')
#     ax.set_xlabel('Percentage Error going from -20% to 20%\n (Zoom between -100% and 100%)')
#     ax.set_ylabel('Density')
#     ax.set_xlim(-20, 20)

#     # Add a grey grid
#     ax.grid(True, which='both', linestyle='--', linewidth=0.5, color='grey')

#     leg = ax.legend(fancybox=True, shadow=True)
#     lined = {}
#     for legline, origline in zip(leg.get_lines(), lines):
#         legline.set_picker(5)
#         lined[legline] = origline

#     fig.canvas.mpl_connect('pick_event', lambda event: toggle_visibility(event, fig, lined))
#     plt.show()



# CON LIMITE SULLE Y
import matplotlib.pyplot as plt
import numpy as np

def plot_percentage_errors(percentage_errors, reference_key):
    global fig, ax
    fig, ax = plt.subplots(figsize=(14, 8))
    lines = []
    
    x_range = np.linspace(-100, 100, 100000)  # Error range from -100% to 100%
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']


    for i, (key, errors) in enumerate(percentage_errors.items()):
        # JUMP M0
        if i == 0:
            continue  # This skips the first iteration

        if key != reference_key and len(errors) > 0:
            mean, std = np.mean(errors), np.std(errors)
            # Calculate PDF of a normal (Gaussian) distribution
            p = 1 / (std * np.sqrt(2 * np.pi)) * np.exp(-((x_range - mean) ** 2) / (2 * std**2))
            line, = ax.plot(x_range, p, linewidth=2, label=f"{key} (μ={mean:.2f}, σ={std:.2f})", color=colors[i % len(colors)])
            lines.append(line)

    if not lines:
        print("No valid data to plot.")
        return

    # Draw the x and y axis lines prominently
    ax.axhline(0, color='black', linewidth=1)
    ax.axvline(0, color='black', linewidth=1)

    # ax.set_title(f'Dragon Mesh Set Created with Decimation\n'
    #              f'Distribution of Percentage Errors for "Dragon" meshes with {key_string} Method',
    #              fontweight='bold')
    ax.set_title(f'Dragon Mesh Set Created with Decimation\n'
                f'Distribution of Percentage Errors for "Dragon" meshes with {key_string} Method',
                fontweight='bold')
    ax.set_xlabel('Percentage Error going from -100% to 100%\n (Zoom between -15% and 15%)')
    ax.set_ylabel('Density')
    ax.set_xlim(-15, 15)
    
    # Set visualization limits for y-axis
    ax.set_ylim(0, 5)  # Limits visualization to y-values between 0 and 3

    # Add a grey grid
    ax.grid(True, which='both', linestyle='--', linewidth=0.5, color='grey')

    leg = ax.legend(fancybox=True, shadow=True)
    lined = {}
    for legline, origline in zip(leg.get_lines(), lines):
        legline.set_picker(5)
        lined[legline] = origline

    fig.canvas.mpl_connect('pick_event', lambda event: toggle_visibility(event, fig, lined))
    plt.show()

# To use this function, ensure that 'percentage_errors' and 'reference_key' are properly defined.


# To use this function, make sure to define percentage_errors and reference_key appropriately.


## PROVA
# def plot_percentage_errors(percentage_errors, reference_key):
#     fig, ax = plt.subplots(figsize=(12, 8))
    
#     # Set up colors and lines
#     colors = plt.cm.viridis(np.linspace(0, 1, len(percentage_errors)))
#     color_idx = 0
    
#     # Loop through each file's errors
#     for key, errors in percentage_errors.items():
#         if key == reference_key:
#             continue  # Skip plotting the reference dataset

#         # Calculate statistics
#         mean_error = np.mean(errors)
#         std_error = np.std(errors)

#         # Create x values for plotting the Gaussian distribution
#         x_values = np.linspace(mean_error - 4*std_error, mean_error + 4*std_error, 1000)
#         y_values = norm.pdf(x_values, mean_error, std_error)

#         # tell me the x value for which y is maximum
#         max_y = np.max(y_values)
#         max_x = x_values[np.argmax(y_values)]
#         print(f"Max y value for {key}: {max_y} at x={max_x}")

#         # Plot the Gaussian distribution as a line
#         ax.plot(x_values, y_values, label=f'{key} (μ={mean_error:.2f}%, σ={std_error:.2f}%)', color=colors[color_idx])
#         color_idx += 1

#     # Formatting the plot
#     ax.set_title('Estimated Gaussian Distributions of Percentage Errors')
#     ax.set_xlabel('Percentage Error (%)')
#     ax.set_ylabel('Probability Density')
#     ax.legend(title='Datasets', loc='upper right')
#     ax.grid(True)

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
directory_path = 'data/DISTANCES_dragon_reordered'
# directory_path = 'data/blub500f_distances'

key_string = 'VTP'  # Filter criteria for non-reference files
# reference_key = 'bob_500f_6_VTP_100.txt'  # Reference file - GROUND TRUTH
reference_key = 'M7_reordered_VTP_16.txt'
# reference_key = '5_blub_VTP_663.txt'


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