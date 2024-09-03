# import trimesh
# import os

# # Directory containing the .obj files
# directory = os.path.dirname(os.path.realpath(__file__))

# # File to store the results in sorted order
# output_file = os.path.join(directory, 'sorted_mesh_complexity_scores.txt')

# def compute_metrics(mesh):
#     V = len(mesh.vertices)
#     F = len(mesh.faces)
#     E = len(mesh.edges)
#     A = mesh.area
#     D = 2 * E / V  # average degree calculation
#     return V, F, E, D, A

# def complexity_score(V, F, E, D, A):
#     # Preliminary normalization factors (adjust these based on your data)
#     max_V, max_F, max_E, max_D, max_A = 10000, 20000, 30000, 12, 5000

#     # Weights for each metric
#     w_v, w_f, w_e, w_d, w_a = 0.2, 0.2, 0.2, 0.2, 0.2

#     # Score calculation using normalized values
#     score = (w_v * (V / max_V) +
#              w_f * (F / max_F) +
#              w_e * (E / max_E) +
#              w_d * (D / max_D) +
#              w_a * (A / max_A))
#     return score

# def gather_and_sort_data():
#     data = []
#     print("Starting to process files...")
#     for filename in os.listdir(directory):
#         if filename.endswith('.obj'):
#             filepath = os.path.join(directory, filename)
#             print(f"Processing {filename}...")
#             try:
#                 mesh = trimesh.load_mesh(filepath)
#                 V, F, E, D, A = compute_metrics(mesh)
#                 score = complexity_score(V, F, E, D, A)
#                 data.append((filename, score, V, F, E, D, A))
#                 print(f"Completed processing {filename}.")
#             except Exception as e:
#                 print(f"Failed to process {filename}: {e}")
#     print("Sorting data by complexity score...")
#     # Sort by complexity score, descending
#     data.sort(key=lambda x: x[1], reverse=True)
#     return data

# def write_data_to_file(data):
#     print("Writing data to file...")
#     with open(output_file, 'w') as f:
#         for entry in data:
#             f.write(f"File: {entry[0]}, Score: {entry[1]:.4f}, Vertices: {entry[2]}, Faces: {entry[3]}, Edges: {entry[4]}, Avg Degree: {entry[5]:.2f}, Area: {entry[6]:.2f}\n")
#     print(f"Data has been written to {output_file}")

# def main():
#     data = gather_and_sort_data()
#     write_data_to_file(data)

# if __name__ == "__main__":
#     main()



import trimesh
import os

# Directory containing the .obj files
directory = os.path.dirname(os.path.realpath(__file__))

# File to store the results in sorted order
output_file = os.path.join(directory, 'sorted_mesh_complexity_scores_dynamic.txt')

def compute_metrics(mesh):
    V = len(mesh.vertices)
    F = len(mesh.faces)
    E = len(mesh.edges)
    A = mesh.area
    D = 2 * E / V  # average degree calculation
    return V, F, E, D, A

def compute_max_values(data):
    # Compute max values across all meshes
    print("---- Computing max values ----")
    max_V = max(data, key=lambda x: x[1])[1]
    max_F = max(data, key=lambda x: x[2])[2]
    max_E = max(data, key=lambda x: x[3])[3]
    max_D = max(data, key=lambda x: x[4])[4]
    max_A = max(data, key=lambda x: x[5])[5]
    print(f"Max values: V={max_V}, F={max_F}, E={max_E}, D={max_D}, A={max_A}")
    return max_V, max_F, max_E, max_D, max_A

def complexity_score(V, F, E, D, A, max_V, max_F, max_E, max_D, max_A):
    # Weights for each metric
    w_v, w_f, w_e, w_d, w_a = 0.2, 0.2, 0.2, 0.2, 0.2

    # Score calculation using normalized values
    score = (w_v * (V / max_V) +
             w_f * (F / max_F) +
             w_e * (E / max_E) +
             w_d * (D / max_D) +
             w_a * (A / max_A))
    return score

def gather_data_and_compute_max():
    data = []
    for filename in os.listdir(directory):
        if filename.endswith('.obj'):
            filepath = os.path.join(directory, filename)
            print(f"Processing {filename}...")
            try:
                mesh = trimesh.load_mesh(filepath)
                metrics = compute_metrics(mesh)
                data.append((filename, *metrics))
                print(f"Completed processing {filename}.")
            except Exception as e:
                print(f"Failed to process {filename}: {e}")
    return data, compute_max_values(data)

def sort_and_write_data(data, max_values):
    # Sort by complexity score, descending
    data.sort(key=lambda x: complexity_score(x[1], x[2], x[3], x[4], x[5], *max_values), reverse=True)
    with open(output_file, 'w') as f:
        for entry in data:
            score = complexity_score(entry[1], entry[2], entry[3], entry[4], entry[5], *max_values)
            f.write(f"File: {entry[0]}, Score: {score:.4f}, Vertices: {entry[1]}, Faces: {entry[2]}, Edges: {entry[3]}, Avg Degree: {entry[4]:.2f}, Area: {entry[5]:.2f}\n")

def main():
    data, max_values = gather_data_and_compute_max()
    sort_and_write_data(data, max_values)
    print(f"Data has been written to {output_file}")

if __name__ == "__main__":
    main()
