# import os
# import pymeshlab
# import numpy as np
# from scipy.stats import skew, kurtosis

# def compute_mesh_complexity(file_path):
#     ms = pymeshlab.MeshSet()
#     ms.load_new_mesh(file_path)
#     m = ms.current_mesh()
    
#     # Basic topology measures
#     num_vertices = m.vertex_number()
#     num_faces = m.face_number()
#     num_edges = m.edge_number()
    
#     # Compute Euler characteristic
#     euler_characteristic = num_vertices - num_edges + num_faces
    
#     # Compute genus (assuming a single connected component)
#     genus = (2 - euler_characteristic) // 2
    
#     # Compute average valence (number of edges per vertex)
#     avg_valence = 2 * num_edges / num_vertices
    
#     # Compute mesh surface area and volume
#     ms.apply_filter('get_geometric_measures')
#     measures = ms.get_geometric_measures()
#     area = measures['surface_area']
#     volume = measures.get('mesh_volume', 0)  # Use 0 if mesh is not closed
    
#     # Compute aspect ratio
#     ms.apply_filter('compute_scalar_by_aspect_ratio_per_face')
#     quality = m.face_quality()
#     aspect_ratio = np.mean(quality)
    
#     # Compute curvature
#     ms.apply_filter('compute_curvature_principal_directions_per_vertex')
#     curvature = m.vertex_quality()
#     mean_curvature = np.mean(curvature)
    
#     # Compute face angle statistics
#     ms.compute_scalar_by_aspect_ratio_per_face()
#     face_angles = ms.get_current_mesh().face_quality()
#     angle_min = np.min(face_angles)
#     angle_max = np.max(face_angles)
#     angle_avg = np.mean(face_angles)
#     angle_std = np.std(face_angles)
#     angle_skewness = skew(face_angles)
#     angle_kurtosis = kurtosis(face_angles)
    
#     # Compute edge length statistics
#     ms.compute_quality_edge_by_length()
#     edge_lengths = ms.get_current_mesh().edge_quality()
#     length_min = np.min(edge_lengths)
#     length_max = np.max(edge_lengths)
#     length_avg = np.mean(edge_lengths)
#     length_std = np.std(edge_lengths)
#     length_skewness = skew(edge_lengths)
#     length_kurtosis = kurtosis(edge_lengths)
    
#     # Compute non-manifold elements
#     ms.compute_topological_measures()
#     topo_measures = ms.get_topological_measures()
#     non_manifold_edges = topo_measures.get('edges_nmanifold', 0)
#     non_manifold_vertices = topo_measures.get('vert_nmanifold', 0)
    
#     # Compute connected components
#     connected_components = topo_measures.get('connected_components', 1)
    
#     # Compute a simple complexity score (you may want to adjust this)
#     complexity_score = (num_vertices * num_faces * avg_valence) / (area * (abs(euler_characteristic) + 1))
    
#     return {
#         'file_name': os.path.basename(file_path),
#         'num_vertices': num_vertices,
#         'num_faces': num_faces,
#         'num_edges': num_edges,
#         'euler_characteristic': euler_characteristic,
#         'genus': genus,
#         'avg_valence': avg_valence,
#         'surface_area': area,
#         'volume': volume,
#         'aspect_ratio': aspect_ratio,
#         'mean_curvature': mean_curvature,
#         'angle_statistics': {
#             'min': angle_min,
#             'max': angle_max,
#             'avg': angle_avg,
#             'std': angle_std,
#             'skewness': angle_skewness,
#             'kurtosis': angle_kurtosis
#         },
#         'edge_length_statistics': {
#             'min': length_min,
#             'max': length_max,
#             'avg': length_avg,
#             'std': length_std,
#             'skewness': length_skewness,
#             'kurtosis': length_kurtosis
#         },
#         'non_manifold_edges': non_manifold_edges,
#         'non_manifold_vertices': non_manifold_vertices,
#         'connected_components': connected_components,
#         'complexity_score': complexity_score
#     }

# def analyze_mesh_folder(folder_path):
#     mesh_complexities = []
    
#     for file_name in os.listdir(folder_path):
#         if file_name.endswith('.obj'):
#             file_path = os.path.join(folder_path, file_name)
#             complexity_info = compute_mesh_complexity(file_path)
#             mesh_complexities.append(complexity_info)
    
#     # Sort meshes by complexity score in descending order
#     sorted_meshes = sorted(mesh_complexities, key=lambda x: x['complexity_score'], reverse=True)
    
#     return sorted_meshes

# def print_mesh_analysis(mesh):
#     print(f"File: {mesh['file_name']}")
#     print(f"Complexity Score: {mesh['complexity_score']:.2f}")
#     print(f"Vertices: {mesh['num_vertices']}, Faces: {mesh['num_faces']}, Edges: {mesh['num_edges']}")
#     print(f"Euler Characteristic: {mesh['euler_characteristic']}, Genus: {mesh['genus']}")
#     print(f"Average Valence: {mesh['avg_valence']:.2f}")
#     print(f"Surface Area: {mesh['surface_area']:.2f}, Volume: {mesh['volume']:.2f}")
#     print(f"Aspect Ratio: {mesh['aspect_ratio']:.4f}")
#     print(f"Mean Curvature: {mesh['mean_curvature']:.4f}")
#     print("Face Angle Statistics:")
#     for key, value in mesh['angle_statistics'].items():
#         print(f"  {key}: {value:.4f}")
#     print("Edge Length Statistics:")
#     for key, value in mesh['edge_length_statistics'].items():
#         print(f"  {key}: {value:.4f}")
#     print(f"Non-manifold Edges: {mesh['non_manifold_edges']}, Non-manifold Vertices: {mesh['non_manifold_vertices']}")
#     print(f"Connected Components: {mesh['connected_components']}")
#     print("---")

# if __name__ == "__main__":
#     folder_path = os.getcwd()  # Use the current directory
#     results = analyze_mesh_folder(folder_path)

#     # Print results
#     for mesh in results:
#         print_mesh_analysis(mesh)

#     # Optionally, save results to a file
#     with open('mesh_analysis_results.txt', 'w') as f:
#         for mesh in results:
#             f.write(f"File: {mesh['file_name']}\n")
#             f.write(f"Complexity Score: {mesh['complexity_score']:.2f}\n")
#             f.write(f"Vertices: {mesh['num_vertices']}, Faces: {mesh['num_faces']}, Edges: {mesh['num_edges']}\n")
#             f.write("---\n")

#     print(f"Analysis complete. Results saved to 'mesh_analysis_results.txt'")



import os
import pymeshlab
import numpy as np
import trimesh

def compute_mesh_complexity(file_path):
    print(f"Loading mesh from {file_path}...")
    mesh = trimesh.load(file_path, force='mesh')

    num_vertices = len(mesh.vertices)
    num_faces = len(mesh.faces)
    num_edges = len(mesh.edges)

    euler_characteristic = num_vertices - num_edges + num_faces
    genus = max(0, (2 - euler_characteristic) // 2)  # Ensure non-negative

    avg_valence = 2 * num_edges / num_vertices if num_vertices > 0 else 0

    ms = pymeshlab.MeshSet()
    ms.load_new_mesh(file_path)
    ms.apply_filter('get_geometric_measures')
    measures = ms.get_geometric_measures()
    area = measures.get('surface_area', 0)
    volume = measures.get('mesh_volume', 0)

    complexity_score = ((num_vertices + num_faces + num_edges) / (area + 1)) if area > 0 else 0

    return {
        'file_name': os.path.basename(file_path),
        'num_vertices': num_vertices,
        'num_faces': num_faces,
        'num_edges': num_edges,
        'euler_characteristic': euler_characteristic,
        'genus': genus,
        'avg_valence': avg_valence,
        'surface_area': area,
        'volume': volume,
        'complexity_score': complexity_score
    }

def analyze_mesh_folder(folder_path):
    mesh_complexities = []
    
    for file_name in os.listdir(folder_path):
        if file_name.endswith('.obj'):
            file_path = os.path.join(folder_path, file_name)
            print(f"Analyzing {file_name}...")
            try:
                complexity_info = compute_mesh_complexity(file_path)
                mesh_complexities.append(complexity_info)
            except Exception as e:
                print(f"Error processing {file_name}: {str(e)}")
    
    # Sort meshes by complexity score in descending order
    sorted_meshes = sorted(mesh_complexities, key=lambda x: x['complexity_score'], reverse=True)
    
    print("Sorting completed. Preparing to display results...")
    return sorted_meshes

def print_mesh_analysis(mesh):
    print(f"File: {mesh['file_name']}")
    print(f"Complexity Score: {mesh['complexity_score']:.2f}")
    print(f"Vertices: {mesh['num_vertices']}, Faces: {mesh['num_faces']}, Edges: {mesh['num_edges']}")
    print(f"Euler Characteristic: {mesh['euler_characteristic']}, Genus: {mesh['genus']}")
    print(f"Average Valence: {mesh['avg_valence']:.2f}")
    print(f"Surface Area: {mesh['surface_area']:.2f}, Volume: {mesh['volume']:.2f}")
    print("---")

if __name__ == "__main__":
    folder_path = os.getcwd()  # Use the current directory
    print("Starting analysis...")
    results = analyze_mesh_folder(folder_path)

    # Print results
    for mesh in results:
        print_mesh_analysis(mesh)

    # Save results to a file
    with open('mesh_analysis_results.txt', 'w') as f:
        for mesh in results:
            f.write(f"File: {mesh['file_name']}\n")
            f.write(f"Complexity Score: {mesh['complexity_score']:.2f}\n")
            f.write(f"Vertices: {mesh['num_vertices']}, Faces: {mesh['num_faces']}, Edges: {mesh['num_edges']}\n")
            f.write("---\n")

    print(f"Analysis complete. Results saved to 'mesh_analysis_results.txt'")
