import pymeshlab as ml
import numpy as np
import os

def repair_mesh(mesh_input, mesh_output):
    # Load the mesh
    # ms = pymeshlab.MeshSet()
    # ms.load_new_mesh(mesh_input)

    # Function to parse the .obj file
    def parse_obj_file(file_path):
        vertices = []
        faces = []
        with open(file_path, 'r') as file:
            for line in file:
                parts = line.split()
                if line.startswith('v '):
                    vertices.append([float(parts[1]), float(parts[2]), float(parts[3])])
                elif line.startswith('f '):
                    # Faces in OBJ are 1-indexed, pymeshlab uses 0-indexed
                    face_indices = [int(idx.split('/')[0]) - 1 for idx in parts[1:]]
                    faces.append(face_indices)
        return np.array(vertices, dtype=np.float32), np.array(faces, dtype=np.int32)

    # Read vertices and faces from the original file
    vertices, faces = parse_obj_file(mesh_input)

    # Create a new Mesh with vertices and faces
    mesh = ml.Mesh(vertex_matrix=vertices, face_matrix=faces)

    # Add the mesh to the MeshSet
    ms = ml.MeshSet()
    ms.add_mesh(mesh, mesh_input)
    
    # Initial diagnostics
    print("Initial diagnostics:")
    print("Number of vertices:", ms.current_mesh().vertex_number())
    print("Number of faces:", ms.current_mesh().face_number())
    # print("Number of holes:", ms.compute_number_of_holes())
    # print("Number of connected components:", ms.compute_number_of_connected_components())

    # Close holes
    print("Attempting to close holes...")
    ms.meshing_close_holes(maxholesize=100)  # Adjust 'maxholesize' as needed
    # print("Holes after operation:", ms.compute_number_of_holes())

    # Remove non-manifold edges
    print("Removing non-manifold edges...")
    ms.meshing_repair_non_manifold_edges()

    # Remove non-manifold vertices
    print("Removing non-manifold vertices...")
    ms.meshing_repair_non_manifold_vertices()

    # Final diagnostics
    print("Final diagnostics:")
    print("Number of vertices:", ms.current_mesh().vertex_number())
    print("Number of faces:", ms.current_mesh().face_number())
    # print("Number of holes:", ms.compute_number_of_holes())
    # print("Number of connected components:", ms.compute_number_of_connected_components())

    # Save the repaired mesh
    # ms.save_current_mesh(mesh_output)
    # print(f"Repaired mesh saved as {mesh_output}")

    # Save the repaired mesh
    try:
        ms.save_current_mesh(mesh_output)
        print(f"Repaired mesh saved as {mesh_output}")
    except Exception as e:
        print(f"Failed to save the mesh: {e}")

# if __name__ == "__main__":
#     import sys
#     if len(sys.argv) != 3:
#         print("Usage: python repair_holes.py [input_mesh.ply] [output_mesh.ply]")
#     else:
#         input_mesh = sys.argv[1]
#         output_mesh = sys.argv[2]
#         repair_mesh(input_mesh, output_mesh)

# call the function with the input and output mesh paths
# repair_mesh('Esperimento_1/data/blub/blub_tri.obj', 'blub_tri_output_mesh.obj')

# Ensure the output folder exists
output_folder = 'repaired_blub'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Process all .obj files in the input folder
#input_folder = 'Esperimento_1/data/altro'
input_folder = '/Users/giuliabenvenuto/Desktop/organic_meshes/blub'
for file in os.listdir(input_folder):
    if file.endswith('.obj'):
        input_mesh = os.path.join(input_folder, file)
        output_mesh = os.path.join(output_folder, file)
        print(f"Processing {input_mesh}...")
        repair_mesh(input_mesh, output_mesh)
        print("\n")

print("All meshes repaired")
