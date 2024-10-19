import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the CSV file
mean_times_df = pd.read_csv('../DISTANCES_OK/mean_times_FINAL_MOD.csv')  # Make sure to set the correct path to your CSV file

# Calculate the size of the mesh
mean_times_df['Mesh Size'] = mean_times_df['Vertices']


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

# Create a scatter plot for Query Time with logarithmic scales on both axes
plt.figure(figsize=(10, 6))
sns.scatterplot(data=mean_times_df, x='Mesh Size', y='Query Time', hue='Method', palette=color_mapping, s=50, marker='o')
plt.xscale('log')  # Set the x-axis to a logarithmic scale
plt.yscale('log')  # Set the y-axis to a logarithmic scale
plt.title('Mesh Size vs Query Time')
plt.grid(True, linestyle='--', alpha=0.6, which='both')
plt.xlabel('Mesh Size (Vertices)')
plt.ylabel('Query Time')
plt.legend(title='Method', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Save the plot
plt.savefig('mesh_size_vs_query_FINAL.png', dpi=300)

# Show the plot
plt.show()

# Filter rows where preprocessing is zero
filtered_df = mean_times_df[mean_times_df['Preprocessing Time'] > 0]

# Create a scatter plot for Preprocessing Time with logarithmic scales on both axes
plt.figure(figsize=(10, 6))
sns.scatterplot(data=filtered_df, x='Mesh Size', y='Preprocessing Time', hue='Method', palette=color_mapping, s=50, marker='o')
plt.xscale('log')  # Set the x-axis to a logarithmic scale
plt.yscale('log')  # Set the y-axis to a logarithmic scale
plt.title('Mesh Size vs Preprocessing Time')
plt.grid(True, linestyle='--', alpha=0.6, which='both')
plt.xlabel('Mesh Size (Vertices)')
plt.ylabel('Preprocessing Time')
plt.legend(title='Method', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Save the plot
plt.savefig('mesh_size_vs_preprocessing_FINAL.png', dpi=300)

# Show the plot
plt.show()
