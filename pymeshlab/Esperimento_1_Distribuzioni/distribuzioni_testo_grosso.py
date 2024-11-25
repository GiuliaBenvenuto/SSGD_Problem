# ----------- VERSIONE CON ERRORE PERCENTUALE ------------
# ----------- NUMERO VERTICI NELLA LEGENDA ------------
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm 

def read_distances_from_file(file_path):
    distances = []
    num_vertices = None  # Initialize to store the number of vertices
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                if 'Number of vertices:' in line:
                    num_vertices = int(line.split('Number of vertices:')[1].split(',')[0].strip())
            else:
                distances.append(float(line.strip()))
    return np.array(distances), num_vertices


def process_directory(directory_path, key_string, reference_filename):
    file_data = {}
    mesh_info = {}
    reference_included = False
    
    for filename in sorted(os.listdir(directory_path)):
        if filename.endswith('.txt') and (key_string in filename or filename == reference_filename):
            file_path = os.path.join(directory_path, filename)
            distances, num_vertices = read_distances_from_file(file_path)
            file_data[filename] = distances
            mesh_info[filename] = num_vertices
            if filename == reference_filename:
                reference_included = True
    
    # Ensure the reference file is always included if not already
    if not reference_included:
        reference_file_path = os.path.join(directory_path, reference_filename)
        if os.path.exists(reference_file_path):
            distances = read_distances_from_file(reference_file_path)
            file_data[reference_filename] = distances
            mesh_info[reference_filename] = num_vertices
        else:
            raise FileNotFoundError(f"Reference file {reference_filename} not found in the directory.")
    
    return file_data, mesh_info



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


# CON LIMITE SULLE Y
def plot_percentage_errors(percentage_errors, reference_key, mesh_info, plot_name="plot_result"):
    global fig, ax
    fig, ax = plt.subplots(figsize=(12, 8))
    lines = []
    
    x_range = np.linspace(-100, 100, 100000)  # Error range from -100% to 100%
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']

    for i, (key, errors) in enumerate(percentage_errors.items()):
        if i == 0:
            continue  # This skips the first iteration

        if key != reference_key and len(errors) > 0:
            mean, std = np.mean(errors), np.std(errors)
            # Calculate PDF of a normal (Gaussian) distribution
            p = 1 / (std * np.sqrt(2 * np.pi)) * np.exp(-((x_range - mean) ** 2) / (2 * std**2))
            
            vertices_label = f"#V: {mesh_info[key]:,}" if key in mesh_info and mesh_info[key] is not None else ""
            label = f"(μ={mean:.2f}, σ={std:.2f}) {vertices_label}"
            
            line, = ax.plot(x_range, p, linewidth=2, label=label, color=colors[i % len(colors)])
            lines.append(line)

    if not lines:
        print("No valid data to plot.")
        return

    # Draw the x and y axis lines prominently
    ax.axhline(0, color='black', linewidth=2)
    ax.axvline(0, color='black', linewidth=2)

    # Set title with increased font size
    ax.set_title(f'Thai Statue Mesh Set \n Distribution of Percentage Errors with k-Ring Graph Method',
                 fontsize=20, fontweight='bold', pad = 10)

    # Set axis labels with larger font size
    ax.set_xlabel('Percentage Error going from -100% to 100%\n(Zoom between -15% and 15%)', fontsize=17)
    ax.set_ylabel('Density', fontsize=17)

    # Adjust axis limits for better visualization
    ax.set_xlim(-15, 15)
    ax.set_ylim(0, 2)

    # Increase tick label size
    ax.tick_params(axis='both', labelsize=16)

    # Add a grey grid with finer lines
    ax.grid(True, which='both', linestyle='--', linewidth=1, color='grey')

    # Modify legend with larger text
    leg = ax.legend(fancybox=True, shadow=True, fontsize=15)
    lined = {}
    for legline, origline in zip(leg.get_lines(), lines):
        legline.set_picker(5)
        lined[legline] = origline

    fig.canvas.mpl_connect('pick_event', lambda event: toggle_visibility(event, fig, lined))
    
    # Determine save path and create directory if needed
    script_directory = os.path.dirname(os.path.realpath(__file__))  # Get the directory of the script
    save_directory = os.path.join(script_directory, 'thai_plot')
    if not os.path.exists(save_directory):
        os.makedirs(save_directory)

    save_path = os.path.join(save_directory, f'{plot_name}.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"Plot saved to {save_path}")



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
# directory_path = 'data/DIST_THAI_DECIMATION'
# directory_path = 'data/blub500f_distances'
# directory_path = 'data/DIST_BLUB/663'
# directory_path = 'data/DIST_BUNNY_500'
# directory_path = 'data/DIST_DRAGO_DECIMATION'
directory_path = 'RIFACCIO_PER_TESI/THAI/DIST_THAI_DECIMATION'

key_string = 'Extended'  # Filter criteria for non-reference files

# reference_key = 'bob_500f_6_VTP_100.txt'  # Reference file - GROUND TRUTH
# reference_key = 'M7_reordered_VTP_16.txt'
# reference_key = '5_blub_VTP_663.txt'
# reference_key = '6_bunny500_VTP_100.txt'
reference_key = '8_thai_VTP_2.txt'


# Process directory, ensuring reference file is included
file_data, mesh_info = process_directory(directory_path, key_string, reference_key)


# Compute and plot percentage errors
if reference_key in file_data:
    percentage_errors = compute_percentage_errors(file_data, reference_key)
    for filename, errors in percentage_errors.items():
        if filename != reference_key:
            print(f"Percentage Errors for {filename}:")
            print(errors)
            print(f"Mean: {np.mean(errors):.2f}%, Std: {np.std(errors):.2f}%")
            print()
    plot_percentage_errors(percentage_errors, reference_key, mesh_info, plot_name=f"thai_{key_string}")
else:
    print(f"Reference file {reference_key} was not found.")