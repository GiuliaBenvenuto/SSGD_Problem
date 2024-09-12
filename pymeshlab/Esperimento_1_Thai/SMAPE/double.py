import pandas as pd

# Load your CSV file
data = pd.read_csv("smape_THAI_56.csv")

# Define the columns to modify (all SMAPE columns)
smape_columns = [
    "SMAPE_VTP", "SMAPE_Trettner", "SMAPE_FastMarching", 
    "SMAPE_Heat", "SMAPE_Geotangle", "SMAPE_Edge", 
    "SMAPE_Extended", "SMAPE_Lanthier"
]

# Multiply each SMAPE value by 2, ignoring missing values ('-')
for column in smape_columns:
    data[column] = pd.to_numeric(data[column], errors='coerce') * 2

# Save the modified DataFrame back to a CSV file
data.to_csv("smape_THAI_56_corretto.csv", index=False)
