import pandas as pd
import matplotlib.pyplot as plt

# Function to plot the data from a given CSV file with "N" and "SMAPE"
def plot_smape_error_n(csv_file_path):
    # Read the CSV file
    data = pd.read_csv(csv_file_path)
    
    # Plotting the data
    plt.figure(figsize=(14, 8))
    plt.plot(data['N'], data['SMAPE'], marker='o')
    
    # Adjust tick label font size
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    # Adding labels for each point with adjusted position
    for n, smape in zip(data['N'], data['SMAPE']):
        plt.text(n + 0.2, smape + 0.02, f'{smape:.4f}', fontsize=12, ha='center')
    
    # Setting labels and title
    plt.xlabel('N: Number of Steiner Points per Edge',fontsize=16)
    plt.ylabel('SMAPE Error', fontsize=16)
    plt.title('SMAPE Error vs N Value for the Fertility Organic Mesh', fontweight='bold', fontsize=20)
    plt.grid(True)
    
    # Save the plot as a PNG file in the same directory as the input file
    output_path = csv_file_path.rsplit('.', 1)[0] + '_smape_plot_PROVA.png'
    plt.savefig(output_path, dpi=300)
    plt.show()

# Example usage with a hypothetical file path
plot_smape_error_n('../DIST/n/N_SMAPE_Fertility.csv')




