# import matplotlib.pyplot as plt
# import pandas as pd

# # Load the data
# csv_path = "times_and_smape_93.csv"
# df = pd.read_csv(csv_path)

# # Calculate Query Time per Vertex
# df['Query Time per Vertex'] = df['Query Time'] / df['Vertices']

# # Define the color mapping for methods
# color_mapping = {
#     'Edge': '#1f77b4',
#     'Extended': '#ff7f0e',
#     'Fast Marching': '#2ca02c',
#     'Geotangle': '#d62728',
#     'Heat': '#9467bd',
#     'Lanthier': '#8c564b',
#     'Trettner': '#e377c2',
#     'VTP': '#7f7f7f'
# }

# # Map colors to methods
# df['Color'] = df['Method'].map(color_mapping)

# # Plotting the scatter plot with circular markers with white borders
# plt.figure(figsize=(12, 8))

# # Scatter plot with circular markers, white borders, and logarithmic scales
# for method, group in df.groupby('Method'):
#     plt.scatter(group['Query Time per Vertex'], group['SMAPE'], 
#                 color=group['Color'].iloc[0], edgecolors='white', label=method, s=100, marker='o', linewidth=0.5)

# plt.xlabel('Query Time / Number of Vertices')
# plt.ylabel('SMAPE')
# plt.xscale('log')
# plt.yscale('log')
# plt.title('Scatter Plot of Query Time per Vertex vs SMAPE (Logarithmic Scales)')
# plt.grid(True, linestyle='--', alpha=0.7)
# plt.legend(title='Method')
# plt.tight_layout()


# # Save the plot as an image file in the same directory
# output_path = "plot_query_time_vs_smape.png"
# plt.savefig(output_path, dpi=300, bbox_inches='tight')

# # Show the plot
# plt.show()


import matplotlib.pyplot as plt
import pandas as pd

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

# Plotting the scatter plot with circular markers with white borders
plt.figure(figsize=(14, 10))

# Scatter plot with circular markers, white borders, and logarithmic scales
for method, group in df.groupby('Method'):
    plt.scatter(group['Query Time per Vertex'], group['SMAPE'], 
                color=group['Color'].iloc[0], edgecolors='white', label=method, s=100, marker='o', linewidth=0.5)

    # Adding labels to each point
    # for _, row in group.iterrows():
    #     plt.text(row['Query Time per Vertex'], row['SMAPE'], 
    #              row['Mesh'], fontsize=8, ha='right', va='bottom', alpha=0.7)

plt.xlabel('Query Time / Number of Vertices', fontsize=14)
plt.ylabel('SMAPE', fontsize=14)
plt.xscale('log')
plt.yscale('log')
plt.title('Scatter Plot of Query Time per Vertex vs SMAPE (Logarithmic Scales)', fontsize=16, fontweight='bold')
plt.grid(True, linestyle='--', alpha=0.6, which='both')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(title='Method', fontsize=12, title_fontsize=14, loc='upper left')
plt.tight_layout()

# Save the plot as an image file in the same directory
output_path = "plot_query_time_vs_smape.png"
plt.savefig(output_path, dpi=300, bbox_inches='tight')

# Show the plot
plt.show()

print(f"Plot saved as {output_path}")
