# ----------------------------
import pymeshlab as ml
import numpy as np

def generate_refined_meshes(mesh_base_name, num_subdivisions):
    print(f"----- Generating refined meshes for {mesh_base_name} with {num_subdivisions} subdivisions -----")
    mesh_name = f'{mesh_base_name}.obj'

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

    # Read vertices and faces from file
    vertices, faces = parse_obj_file(mesh_name)

    # Create a new Mesh with vertices and faces
    mesh = ml.Mesh(vertex_matrix=vertices, face_matrix=faces)

    # Add the mesh to the MeshSet
    ms = ml.MeshSet()
    ms.add_mesh(mesh, mesh_name)

    # print 20 vertices coords of the mesh
    # for i in range(20):
    #     print(f"Vertex {i}: {ms.current_mesh().vertex_matrix()[i]}")
    
    # print 20 faces of the mesh
    # for i in range(20):
    #     print(f"Face {i}: {ms.current_mesh().face_matrix()[i]}")


    # Apply the subdivision filter multiple times
    for i in range(num_subdivisions):
        ms.apply_filter('meshing_surface_subdivision_loop', loopweight='Loop', iterations=1, threshold=ml.PercentageValue(0), selected=False)
        current_vertices = ms.current_mesh().vertex_matrix()
        print(f"After subdivision {i+1}: {current_vertices.shape[0]} vertices, sample vertex: {current_vertices[0]}")

        save_path = f'{mesh_base_name}_subdiv_{i+1}.obj'
        
        # Save without vertex normals and colors:
        # ms.save_current_mesh(save_path)
        ms.save_current_mesh(save_path, save_vertex_normal=False, save_vertex_color=False)
    
        print(f"Mesh saved to {save_path}")

        # Clear and reload the mesh
        ms.clear()
        ms.load_new_mesh(save_path)
        print(f"Reloaded mesh vertices: {ms.current_mesh().vertex_matrix()[0]}")

    print(f"Final mesh has {ms.current_mesh().vertex_number()} vertices")


# Use the function to generate refined meshes
generate_refined_meshes('data/altro/bunny_1k', 3)
# generate_refined_meshes('data/bob/bob_tri', 3)
# generate_refined_meshes('data/altro/Nefertiti', 3)







