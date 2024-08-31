import numpy as np
import trimesh
import os

class Mesh:
    def __init__(self, vertices, faces):
        self.vertices = np.array(vertices)
        self.faces = np.array(faces)

def load_mesh(file_path):
    print(f"Loading mesh from {file_path}...")
    mesh = trimesh.load(file_path, process=False)
    print(f"Loaded {file_path}")
    return Mesh(vertices=mesh.vertices.tolist(), faces=mesh.faces.tolist())

def save_mesh(mesh, file_path):
    print(f"Saving mesh to {file_path}...")
    trimesh.Trimesh(vertices=mesh.vertices, faces=mesh.faces).export(file_path)
    print(f"Saved {file_path}")

def create_mapping(previous_mesh, current_mesh):
    print("Creating vertex mapping...")
    previous_vertices = np.array(previous_mesh.vertices)
    current_vertices = np.array(current_mesh.vertices)

    # Create a dictionary for fast lookup
    vertex_to_index = {tuple(v): i for i, v in enumerate(previous_vertices)}

    # Mapping of vertex indices from current_mesh to previous_mesh
    map_vector = []
    new_vertices = []

    for vertex in current_vertices:
        vertex_tuple = tuple(vertex)
        if vertex_tuple in vertex_to_index:
            map_vector.append(vertex_to_index[vertex_tuple])
        else:
            map_vector.append(len(previous_vertices) + len(new_vertices))
            new_vertices.append(vertex)

    print("Vertex mapping completed.")
    return map_vector, new_vertices

def reorder_mesh(previous_mesh, current_mesh, map_vector, new_vertices):
    print("Reordering vertices and faces...")
    # Combine previous vertices with new vertices
    reordered_vertices = np.vstack((previous_mesh.vertices, new_vertices))

    # Reorder faces of current_mesh
    reordered_faces = [[map_vector[vertex_index] for vertex_index in face] for face in current_mesh.faces]

    print("Reordering completed.")
    return Mesh(reordered_vertices, reordered_faces)

def process_meshes(input_folder, output_folder):
    print("Processing meshes...")
    # Assume file names are M0.obj, M1.obj, ..., M8.obj
    file_paths = [os.path.join(input_folder, f'M{i}.obj') for i in range(8)]  # Changed to 9 to include M8.obj

    print("Loading all meshes...")
    # Load all meshes
    meshes = [load_mesh(file_path) for file_path in file_paths]

    # List for reordered meshes
    mapped_meshes = [meshes[0]]  # Initialize with M0

    # Process each pair of meshes
    for i in range(1, len(meshes)):
        print(f"Processing mesh {i}...")
        previous_mesh = mapped_meshes[-1]
        current_mesh = meshes[i]

        # Create mapping between previous_mesh and current_mesh
        map_vector, new_vertices = create_mapping(previous_mesh, current_mesh)

        # Reorder current_mesh to match previous_mesh and add new vertices
        reordered_mesh = reorder_mesh(previous_mesh, current_mesh, map_vector, new_vertices)

        # Add reordered mesh to the list
        mapped_meshes.append(reordered_mesh)

    print("Saving all reordered meshes...")
    # Save all reordered meshes
    os.makedirs(output_folder, exist_ok=True)
    for i, mesh in enumerate(mapped_meshes):
        save_mesh(mesh, os.path.join(output_folder, f'M{i}_reordered.obj'))

    print("All meshes processed and saved.")

if __name__ == "__main__":
    input_folder = "asian_dragon_decimated"
    output_folder = "dragon_ordered"
    process_meshes(input_folder, output_folder)