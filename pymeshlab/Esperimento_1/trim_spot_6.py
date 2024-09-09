import pandas as pd

# Load the data from the CSV file
data = pd.read_csv('./data/GROUND_TRUTH_GIUSTE/SPOT_GT.csv')

# Extract the 'vertex_174' column
vertex_174_data = data['vertex_174']

# Save to a text file, one value per line
vertex_174_data.to_csv('./vertex_174_data.txt', index=False, header=False)
