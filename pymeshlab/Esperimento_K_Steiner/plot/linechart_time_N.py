# import matplotlib.pyplot as plt

# # Data from the new dataset
# N = [1, 3, 5, 7, 9]
# query_time = [0.0249958, 0.0628777, 0.101319, 0.149553, 0.190584]
# pre_time = [0.08322, 0.164728, 0.289012, 0.380927, 0.521494]

# # Plot for Preprocessing Time
# plt.figure(figsize=(8, 6))
# plt.plot(N, pre_time, marker='o', color='blue', label="Preprocessing Time")
# for i, txt in enumerate(pre_time):
#     plt.text(N[i] - 0.2, pre_time[i], f"{txt:.5f}", fontsize=8, ha='right', va='bottom')  # Move labels slightly to the left
# plt.xlabel("N: Number of Steiner Points per Edge")
# plt.ylabel("Preprocessing Time (seconds)")
# plt.title("Preprocessing Time vs. N for the Fertility Organic Mesh", fontweight="bold")
# plt.grid(True)
# plt.savefig("preprocessing_time_steiner_plot.png", dpi=300)  # Save the plot with higher resolution
# plt.close()

# # Plot for Query Time
# plt.figure(figsize=(8, 6))
# plt.plot(N, query_time, marker='s', color='green', label="Query Time")
# for i, txt in enumerate(query_time):
#     plt.text(N[i] - 0.2, query_time[i], f"{txt:.5f}", fontsize=8, ha='right', va='bottom')  # Move labels slightly to the left
# plt.xlabel("N: Number of Steiner Points per Edge")
# plt.ylabel("Query Time (seconds)")
# plt.title("Query Time vs. N for the Fertility Organic Mesh", fontweight="bold")
# plt.grid(True)
# plt.savefig("query_time_steiner_plot.png", dpi=300)  # Save the plot with higher resolution
# plt.close()

# print("Plots saved as 'preprocessing_time_steiner_plot.png' and 'query_time_steiner_plot.png' in the current directory.")


# WITH BIGGER TEXT
import matplotlib.pyplot as plt

# Data from the new dataset
N = [1, 3, 5, 7, 9]
query_time = [0.0249958, 0.0628777, 0.101319, 0.149553, 0.190584]
pre_time = [0.08322, 0.164728, 0.289012, 0.380927, 0.521494]

# Plot for Preprocessing Time
plt.figure(figsize=(8, 6))
plt.plot(N, pre_time, marker='o', color='blue', label="Preprocessing Time")
for i, txt in enumerate(pre_time):
    plt.text(N[i] + 0.02, pre_time[i] - 0.005, f"{txt:.5f}", fontsize=12, ha='left', va='top')  
plt.xlabel("N: Number of Steiner Points per Edge", fontsize=14)
plt.ylabel("Preprocessing Time (seconds)", fontsize=14)
plt.title("Preprocessing Time vs. N for the Fertility Organic Mesh", fontweight="bold", fontsize=16)
plt.grid(True)
plt.savefig("preprocessing_time_steiner_big.png", dpi=300)  # Save the plot with higher resolution
plt.close()

# Plot for Query Time
plt.figure(figsize=(8, 6))
plt.plot(N, query_time, marker='s', color='green', label="Query Time")
for i, txt in enumerate(query_time):
    plt.text(N[i] + 0.02, query_time[i] - 0.002, f"{txt:.5f}", fontsize=12, ha='left', va='top') 
plt.xlabel("N: Number of Steiner Points per Edge", fontsize=14)
plt.ylabel("Query Time (seconds)", fontsize=14)
plt.title("Query Time vs. N for the Fertility Organic Mesh", fontweight="bold", fontsize=16)
plt.grid(True)
plt.savefig("query_time_steiner_big.png", dpi=300)  # Save the plot with higher resolution
plt.close()

print("Plots saved as 'preprocessing_time_steiner_plot.png' and 'query_time_steiner_plot.png' in the current directory.")