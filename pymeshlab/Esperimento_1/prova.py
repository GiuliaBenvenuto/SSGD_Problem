import pymeshlab as ml
import numpy as np
from scipy.spatial import KDTree

def generate_refined_meshes(mesh_base_name, num_subdivisions=3):
    # Load the original mesh
    mesh_path = f'{mesh_base_name}.ply'
    ms = ml.MeshSet()
    ms.load_new_mesh(mesh_path)

    # List to store mesh data (vertex and face matrices)
    mesh_data = []

    # Apply Loop subdivision multiple times
    for i in range(num_subdivisions):
        print(f'Applying Loop subdivision {i+1}/{num_subdivisions}')
        ms.apply_filter('meshing_surface_subdivision_loop', loopweight='Loop', iterations=1, threshold=ml.PercentageValue(0))
        current_mesh = ms.current_mesh()
        # Store vertices and faces
        mesh_data.append((current_mesh.vertex_matrix(), current_mesh.face_matrix()))

    # Re-index vertices to match the finest mesh
    fine_mesh_vertices = mesh_data[-1][0]  # Get the vertex matrix of the finest mesh
    fine_mesh_kdtree = KDTree(fine_mesh_vertices)  # Using KDTree for efficient nearest neighbor search

    # Create new meshes with re-indexed vertices
    refined_meshes = []
    for i, (vertices, faces) in enumerate(mesh_data):
        print('Re-indexing vertices for mesh', i + 1)
        if i < num_subdivisions - 1:  # Do not re-index the finest mesh
            _, correspondence = fine_mesh_kdtree.query(vertices)
            matched_vertices = fine_mesh_vertices[correspondence]
            refined_mesh = ml.Mesh(vertex_matrix=matched_vertices, face_matrix=faces)
        else:
            refined_mesh = ml.Mesh(vertex_matrix=vertices, face_matrix=faces)
        refined_meshes.append(refined_mesh)

    # Save each mesh
    for i, mesh in enumerate(refined_meshes):
        ms.clear()
        ms.add_mesh(mesh, set_as_current=True)
        save_path = f'{mesh_base_name}_subdiv_{i+1}.ply'
        print(f'Saving mesh to {save_path}')
        ms.save_current_mesh(save_path)

# Example usage
generate_refined_meshes('data/bunny', 3)
