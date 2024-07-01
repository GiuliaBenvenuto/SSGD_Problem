# # ----------------------------
# import pymeshlab as ml
# import numpy as np

# def generate_refined_meshes(mesh_base_name, num_subdivisions):
#     print(f"----- Generating refined meshes for {mesh_base_name} with {num_subdivisions} subdivisions -----")
#     mesh_name = f'{mesh_base_name}.obj'

#     # Function to parse the .obj file
#     def parse_obj_file(file_path):
#         vertices = []
#         faces = []
#         with open(file_path, 'r') as file:
#             for line in file:
#                 parts = line.split()
#                 if line.startswith('v '):
#                     vertices.append([float(parts[1]), float(parts[2]), float(parts[3])])
#                 elif line.startswith('f '):
#                     # Faces in OBJ are 1-indexed, pymeshlab uses 0-indexed
#                     face_indices = [int(idx.split('/')[0]) - 1 for idx in parts[1:]]
#                     faces.append(face_indices)
#         return np.array(vertices, dtype=np.float32), np.array(faces, dtype=np.int32)

#     # Read vertices and faces from file
#     vertices, faces = parse_obj_file(mesh_name)

#     # Create a new Mesh with vertices and faces
#     mesh = ml.Mesh(vertex_matrix=vertices, face_matrix=faces)

#     # Add the mesh to the MeshSet
#     ms = ml.MeshSet()
#     ms.add_mesh(mesh, mesh_name)

#     # print 20 vertices coords of the mesh
#     # for i in range(20):
#     #     print(f"Vertex {i}: {ms.current_mesh().vertex_matrix()[i]}")
    
#     # print 20 faces of the mesh
#     # for i in range(20):
#     #     print(f"Face {i}: {ms.current_mesh().face_matrix()[i]}")


#     # Apply the subdivision filter multiple times
#     for i in range(num_subdivisions):
#         ms.apply_filter('meshing_surface_subdivision_loop', loopweight='Loop', iterations=1, threshold=ml.PercentageValue(0), selected=False)
#         current_vertices = ms.current_mesh().vertex_matrix()
#         print(f"After subdivision {i+1}: {current_vertices.shape[0]} vertices, sample vertex: {current_vertices[0]}")

#         save_path = f'{mesh_base_name}_subdiv_{i+1}.obj'
        
#         # Save without vertex normals and colors:
#         # ms.save_current_mesh(save_path)
#         ms.save_current_mesh(save_path, save_vertex_normal=False, save_vertex_color=False)
    
#         print(f"Mesh saved to {save_path}")

#         # Clear and reload the mesh
#         ms.clear()
#         ms.load_new_mesh(save_path)
#         print(f"Reloaded mesh vertices: {ms.current_mesh().vertex_matrix()[0]}")

#     print(f"Final mesh has {ms.current_mesh().vertex_number()} vertices")
#     print("\n")


#     # Substitute the vertices of the original mesh vertices with the final mesh vertices
#     # take the original mesh from mesh_base_name.obj
#     vertices_original, _ = parse_obj_file(mesh_name)
#     n_vertices_original = len(vertices_original)
#     print(f"Original mesh has {n_vertices_original} vertices")

#     # Read vertices from data/altro/bunny_1k_subdiv_3.obj
#     # VERTICES THAT I USE TO REPLACE THE ORIGINAL VERTICES
#     vertices_final = ms.current_mesh().vertex_matrix()

#     # Print vertices with high precision
#     # for i, vertex in enumerate(vertices_final):
#     #     if i >= 20:
#     #         break
#     #     print(f"Vertex {i}: {vertex[0]:.16f}, {vertex[1]:.16f}, {vertex[2]:.16f}")


#     n_vertices_final = len(vertices_final)
#     print(f"Final mesh has {n_vertices_final} vertices")

#     # Substitute the vertices of the original mesh with the final mesh vertices
#     save_path = mesh_base_name + '_final.obj'
#     vertex_index = 0
#     with open(save_path, 'w') as file:
#         with open(mesh_name, 'r') as original_file:
#             for line in original_file:
#                 if line.startswith('v ') and vertex_index < n_vertices_original:
#                     # Replace original vertex coordinates with final mesh coordinates
#                     new_vertex = vertices_final[vertex_index]
#                     file.write(f"v {new_vertex[0]:.16f} {new_vertex[1]:.16f} {new_vertex[2]:.16f}\n")
#                     vertex_index += 1
#                 else:
#                     # Copy the line as is for the faces of the mesh
#                     file.write(line)



# # Use the function to generate refined meshes
# generate_refined_meshes('data/altro/bunny_1k', 3)
# # generate_refined_meshes('data/bob/bob_tri', 3)
# # generate_refined_meshes('data/altro/Nefertiti', 3)



import numpy as np
import pymeshlab as ml
import os

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

    # Read vertices and faces from the original file
    vertices, faces = parse_obj_file(mesh_name)

    # Create a new Mesh with vertices and faces
    mesh = ml.Mesh(vertex_matrix=vertices, face_matrix=faces)

    # Add the mesh to the MeshSet
    ms = ml.MeshSet()
    ms.add_mesh(mesh, mesh_name)

    # Apply the subdivision filter multiple times and save each step
    for i in range(num_subdivisions):
        ms.apply_filter('meshing_surface_subdivision_loop', loopweight='Loop', iterations=1, threshold=ml.PercentageValue(0), selected=False)
        current_vertices = ms.current_mesh().vertex_matrix()
        print(f"After subdivision {i+1}: {current_vertices.shape[0]} vertices.")

        save_path = f'{mesh_base_name}_subdiv_{i+1}.obj'
        ms.save_current_mesh(save_path, save_vertex_normal=False, save_vertex_color=False)
        print(f"Mesh saved to {save_path}")

        # Clear and reload the mesh
        ms.clear()
        ms.load_new_mesh(save_path)
        # print(f"Reloaded mesh vertices: {ms.current_mesh().vertex_matrix()[0]}")


    # Substitute the vertices of the original and intermediate meshes with the final refined mesh vertices
    ms.clear()
    ms.load_new_mesh(f'{mesh_base_name}_subdiv_{num_subdivisions}.obj')
    vertices_final = ms.current_mesh().vertex_matrix()

    n_vertices_original = len(vertices)
    print(f"Original mesh has {n_vertices_original} vertices")

    save_path = mesh_base_name + '_final.obj'
    vertex_index = 0
    with open(save_path, 'w') as file:
        with open(mesh_name, 'r') as original_file:
            for line in original_file:
                if line.startswith('v ') and vertex_index < n_vertices_original:
                    # Replace original vertex coordinates with final mesh coordinates
                    new_vertex = vertices_final[vertex_index]
                    file.write(f"v {new_vertex[0]:.16f} {new_vertex[1]:.16f} {new_vertex[2]:.16f}\n")
                    vertex_index += 1
                else:
                    # Copy the line as is for the faces of the mesh
                    file.write(line)

    # Perform substitution on all files except the last subdivided file
    for i in range(num_subdivisions):
        current_subdiv_path = f'{mesh_base_name}_subdiv_{i+1}.obj'
        if i == num_subdivisions - 1:
            continue  # Skip the last subdivision file

        vertices_current, _ = parse_obj_file(current_subdiv_path)
        n_vertices_current = len(vertices_current)
        updated_save_path = f'{mesh_base_name}_subdiv_{i+1}_final.obj'
        vertex_index = 0

        with open(updated_save_path, 'w') as file:
            with open(current_subdiv_path, 'r') as subdiv_file:
                for line in subdiv_file:
                    if line.startswith('v ') and vertex_index < n_vertices_current:
                        new_vertex = vertices_final[vertex_index]
                        file.write(f"v {new_vertex[0]:.16f} {new_vertex[1]:.16f} {new_vertex[2]:.16f}\n")
                        vertex_index += 1
                    else:
                        file.write(line)
        # print(f"Updated {current_subdiv_path} saved to {updated_save_path}")


    # remove intermediate files
    for i in range(num_subdivisions - 1):
        current_subdiv_path = f'{mesh_base_name}_subdiv_{i+1}.obj'
        os.remove(current_subdiv_path)

    # rename the last subdivided file adding the _final suffix
    os.rename(f'{mesh_base_name}_subdiv_{num_subdivisions}.obj', f'{mesh_base_name}_subdiv_{num_subdivisions}_final.obj')

# Use the function to generate refined meshes
generate_refined_meshes('data/altro/bunny_1k', 6)
