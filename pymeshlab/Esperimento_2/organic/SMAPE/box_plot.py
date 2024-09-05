import pandas as pd
import matplotlib.pyplot as plt

# Load the data
data = pd.read_csv('SMAPE_93_Heat.csv')  # Make sure to specify the correct path to your CSV file

# Prepare the data for plotting by melting the dataframe for easier plotting
data_melted = data.melt(id_vars=['MeshName', 'NumVertices'], var_name='Method', value_name='SMAPE')

# Plotting boxplot per method, grouped by meshes using matplotlib
plt.figure(figsize=(14, 8))
ax = plt.boxplot(
    [data_melted[data_melted['Method'] == method]['SMAPE'].dropna() for method in data_melted['Method'].unique()],
    labels=data_melted['Method'].unique(),
    vert=True, patch_artist=True
)
plt.xticks(rotation=45)
plt.title('Distribution of SMAPE Errors Per Method Across Different Meshes')
plt.xlabel('Method')
plt.ylabel('SMAPE Error')
plt.grid(True)
plt.tight_layout()
plt.show()
