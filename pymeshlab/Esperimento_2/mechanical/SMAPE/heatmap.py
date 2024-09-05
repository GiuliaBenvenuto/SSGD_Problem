import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the data
data = pd.read_csv('SMAPE_32.csv')  # Adjust the path to your CSV file

# Prepare the data for plotting by melting and then pivoting for the heatmap
data_melted = data.melt(id_vars=['MeshName', 'NumVertices'], var_name='Method', value_name='SMAPE')
heatmap_data = data_melted.pivot_table(index='MeshName', columns='Method', values='SMAPE', aggfunc='mean')

# Plotting the heatmap
plt.figure(figsize=(12, 10))
ax = sns.heatmap(heatmap_data, annot=True, fmt=".2f", cmap="coolwarm", linewidths=.5)
ax.set_title('Heatmap of SMAPE Errors by Mesh and Method')
ax.set_xlabel('Method')
ax.set_ylabel('Mesh')
plt.xticks(rotation=45)  # Rotate x-axis labels to 45 degrees
plt.tight_layout()  # Adjust layout to make room for label rotation
plt.show()