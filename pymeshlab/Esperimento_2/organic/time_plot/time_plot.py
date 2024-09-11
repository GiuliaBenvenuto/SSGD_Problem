import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the CSV file
mean_times_df = pd.read_csv('../DISTANCES_OK/mean_times_FINAL_MOD.csv')  # Make sure to set the correct path to your CSV file

# Calculate the size of the mesh
mean_times_df['Mesh Size'] = mean_times_df['Vertices']

# Set up the plotting environment
sns.set(style="whitegrid")

# Create a scatter plot for Query Time with logarithmic scales on both axes
plt.figure(figsize=(10, 6))
sns.scatterplot(data=mean_times_df, x='Mesh Size', y='Query Time', hue='Method', palette='tab10', s=50, marker='o')
plt.xscale('log')  # Set the x-axis to a logarithmic scale
plt.yscale('log')  # Set the y-axis to a logarithmic scale
plt.title('Mesh Size vs Query Time')
plt.xlabel('Mesh Size (Vertices)')
plt.ylabel('Query Time')
plt.legend(title='Method', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Save the plot
plt.savefig('mesh_size_vs_query_FINAL.png')

# Show the plot
plt.show()

# Create a scatter plot for Preprocessing Time with logarithmic scales on both axes
plt.figure(figsize=(10, 6))
sns.scatterplot(data=mean_times_df, x='Mesh Size', y='Preprocessing Time', hue='Method', palette='tab10', s=50, marker='o')
plt.xscale('log')  # Set the x-axis to a logarithmic scale
plt.yscale('log')  # Set the y-axis to a logarithmic scale
plt.title('Mesh Size vs Preprocessing Time')
plt.xlabel('Mesh Size (Vertices)')
plt.ylabel('Preprocessing Time')
plt.legend(title='Method', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Save the plot
plt.savefig('mesh_size_vs_preprocessing_FINAL.png')

# Show the plot
plt.show()
