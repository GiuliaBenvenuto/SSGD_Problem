import matplotlib.pyplot as plt

# Data from the dataset
K = [3, 4, 5, 6, 7, 8, 9]
query_time = [0.00614158, 0.00984837, 0.0119733, 0.0164816, 0.0186132, 0.0199324, 0.0212537]
pre_time = [6.47139, 16.8154, 35.6599, 68.1973, 119.679, 200.875, 318.182]

# Plot for Preprocessing Time
plt.figure(figsize=(8, 6))
plt.plot(K, pre_time, marker='o', color='blue', label="Preprocessing Time")
for i, txt in enumerate(pre_time):
    plt.text(K[i] - 0.1, pre_time[i], f"{txt:.2f}", fontsize=12, ha='right', va='bottom')
plt.xlabel("K: Number of Rings",fontsize=14)
plt.ylabel("Preprocessing Time (seconds)",fontsize=14)
plt.title("Preprocessing Time vs. K Value for the Fertility Organic Mesh", fontweight="bold", fontsize=16)
plt.grid(True)
plt.savefig("preprocessing_time_K_big.png", dpi=300)  # Save the plot with higher resolution
plt.close()

# Plot for Query Time
plt.figure(figsize=(8, 6))
plt.plot(K, query_time, marker='s', color='green', label="Query Time")
for i, txt in enumerate(query_time):
    plt.text(K[i] - 0.1, query_time[i], f"{txt:.5f}", fontsize=12, ha='right', va='bottom')
plt.xlabel("K (Number of Rings)", fontsize=14)
plt.ylabel("Query Time (seconds)", fontsize=14)
plt.title("Query Time vs. K Value for the Fertility Organic Mesh", fontweight="bold", fontsize=16)
plt.grid(True)
plt.savefig("query_time_K_big.png", dpi=300)  # Save the plot with higher resolution
plt.close()

print("Plots saved as 'preprocessing_time_plot.png' and 'query_time_plot.png' in the current directory.")
