import pandas as pd
import matplotlib.pyplot as plt

# Function to plot the data from a given CSV file
def plot_smape_error(csv_file_path):
    # Read the CSV file
    data = pd.read_csv(csv_file_path)
    
    # Plotting the data
    plt.figure(figsize=(14, 8))
    plt.plot(data['K'], data['SMAPE'], marker='o')
    
    # Adding labels for each point
    for k, smape in zip(data['K'], data['SMAPE']):
        plt.text(k + 0.1, smape + 0.004, f'{smape:.4f}', fontsize=9, ha='center')  # Adjusted position
    
    # Setting labels and title
    plt.xlabel('K: Number of Rings ')
    plt.ylabel('SMAPE Error')
    plt.title('SMAPE Error vs K Value for the Fertility Organic Mesh', fontweight='bold')
    plt.grid(True)
    
    # Save the plot as a PNG file in the same directory as the input file
    output_path = csv_file_path.rsplit('.', 1)[0] + '_smape_plot.png'
    plt.savefig(output_path)
    plt.show()

# Example usage with a hypothetical file path
plot_smape_error('../DIST/k/K_SMAPE_Fertility.csv')
