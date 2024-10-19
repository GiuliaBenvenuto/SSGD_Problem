import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the CSV file
#file_path = 'avg_edge_length_dragon_16.csv'
file_path = 'avg_edge_length_66.csv'
df = pd.read_csv(file_path)

# Extract AvgEdgeLength and SMAPE values for plotting
avg_edge_length = df['AvgEdgeLength']
methods = ['SMAPE_VTP', 'SMAPE_Trettner', 'SMAPE_FastMarching', 'SMAPE_Heat', 
           'SMAPE_Geotangle', 'SMAPE_Edge', 'SMAPE_Extended', 'SMAPE_Lanthier']

# Define the color mapping as requested
color_mapping = {
    'VTP': '#7f7f7f',
    'Trettner': '#e377c2',
    'Fast Marching': '#2ca02c',
    'Heat': '#9467bd',
    'Geotangle': '#d62728',
    'Edge': '#1f77b4',
    'Extended': '#ff7f0e',
    'Lanthier': '#8c564b'
}

# Corresponding keys for the columns in the data
method_labels = {
    'SMAPE_VTP': 'VTP',
    'SMAPE_Trettner': 'Trettner',
    'SMAPE_FastMarching': 'Fast Marching',
    'SMAPE_Heat': 'Heat',
    'SMAPE_Geotangle': 'Geotangle',
    'SMAPE_Edge': 'Edge',
    'SMAPE_Extended': 'Extended',
    'SMAPE_Lanthier': 'Lanthier'
}

# Plotting with specified colors and dots
sns.set(style="whitegrid")
plt.figure(figsize=(12, 8))

for method, label in method_labels.items():
    # scatter plot
    sns.scatterplot(
        x=avg_edge_length, y=df[method], label='_nolegend_', 
        color=color_mapping[label], s=50, marker='o', edgecolor='white', linewidth=0.5
    )
    # line plot
    sns.lineplot(
        x=avg_edge_length, y=df[method], label=label, 
        color=color_mapping[label], linewidth=1.2
    )

plt.xlabel('Average Edge Length')
plt.ylabel('SMAPE Value')
plt.xscale('log')
plt.yscale('log')
plt.title('SMAPE vs. Average Edge Length for Dragon Mesh Set', fontweight='bold')
plt.legend(title='Method')
plt.tight_layout()
plt.grid(True, which="both", ls="--")

# Save the figure
plt.savefig('smape_vs_avg_edge_length_66.png', dpi=300, bbox_inches='tight')

plt.show()



