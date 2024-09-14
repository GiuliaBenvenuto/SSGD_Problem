import pandas as pd
import matplotlib.pyplot as plt

# Function to plot the data from a given CSV file with "N" and "SMAPE"
def plot_smape_error_n(csv_file_path):
    # Read the CSV file
    data = pd.read_csv(csv_file_path)
    
    # Plotting the data
    plt.figure(figsize=(14, 8))
    plt.plot(data['N'], data['SMAPE'], marker='o')
    
    # Adding labels for each point with adjusted position
    for n, smape in zip(data['N'], data['SMAPE']):
        plt.text(n + 0.2, smape + 0.006, f'{smape:.4f}', fontsize=9, ha='center')
    
    # Setting labels and title
    plt.xlabel('N: Number of Steiner Points per Edge')
    plt.ylabel('SMAPE Error')
    plt.title('SMAPE Error vs N Value for the Fertility Organic Mesh', fontweight='bold')
    plt.grid(True)
    
    # Save the plot as a PNG file in the same directory as the input file
    output_path = csv_file_path.rsplit('.', 1)[0] + '_smape_plot.png'
    plt.savefig(output_path)
    plt.show()

# Example usage with a hypothetical file path
plot_smape_error_n('../DIST/n/N_SMAPE_Fertility.csv')
