import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Load the data
csv_path = "times_and_smape_32.csv"
df = pd.read_csv(csv_path)

# Calculate Query Time per Vertex
df['Query Time per Vertex'] = df['Query Time'] / df['Vertices']

# Define the color mapping for methods
color_mapping = {
    'Edge': '#1f77b4',
    'Extended': '#ff7f0e',
    'Fast Marching': '#2ca02c',
    'Geotangle': '#d62728',
    'Heat': '#9467bd',
    'Lanthier': '#8c564b',
    'Trettner': '#e377c2',
    'VTP': '#7f7f7f'
}

# Map colors to methods
df['Color'] = df['Method'].map(color_mapping)

# Plotting
plt.figure(figsize=(14, 10))

# Scatter plot for non-VTP methods
for method, group in df[df['Method'] != 'VTP'].groupby('Method'):
    plt.scatter(group['Query Time per Vertex'], group['SMAPE'], 
                color=group['Color'].iloc[0], edgecolors='white', 
                label=method, s=100, marker='o', linewidth=0.5)

# Scatter plot for VTP (using a small positive value instead of zero)
vtp_data = df[df['Method'] == 'VTP']
small_positive = 1e-3  # Adjust this value as needed
plt.scatter(vtp_data['Query Time per Vertex'], 
            np.full_like(vtp_data['Query Time per Vertex'], small_positive),
            color=color_mapping['VTP'], edgecolors='white', 
            label='VTP', s=100, marker='o', linewidth=0.5)

plt.xlabel('Query Time / Number of Vertices', fontsize=14)
plt.ylabel('SMAPE', fontsize=14)
plt.xscale('log')
plt.yscale('log')
plt.title('Scatter Plot of Query Time per Vertex vs SMAPE (Logarithmic Scales)', fontsize=16, fontweight='bold')
plt.grid(True, linestyle='--', alpha=0.6, which='both')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(title='Method', fontsize=12, title_fontsize=14, loc='upper right')

# Add a text annotation to indicate where y=0 would be
# plt.text(plt.xlim()[0], small_positive, 'VTP (SMAPE ≈ 0)', 
#          verticalalignment='bottom', horizontalalignment='left',
#          fontsize=10, color='black', bbox=dict(facecolor='white', edgecolor='none', alpha=0.7))

plt.tight_layout()

# Save the plot
output_path = "query_time_vs_smape_plot_updated.png"
plt.savefig(output_path, dpi=300, bbox_inches='tight')

# Show the plot
plt.show()

print(f"Updated plot saved as {output_path}")