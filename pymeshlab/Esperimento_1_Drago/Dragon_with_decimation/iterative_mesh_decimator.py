# import pymeshlab
# import time

# def simplify_mesh(input_mesh_path, output_mesh_basepath):
#     ms = pymeshlab.MeshSet()
#     ms.load_new_mesh(input_mesh_path)

#     print(f"Mesh loaded successfully. Initial number of faces: {ms.current_mesh().face_number()}")

#     mesh_sequence = [ms.current_mesh()]
#     target_faces = 200  # Target number of faces
#     current_faces = ms.current_mesh().face_number()

#     i = 0
#     while current_faces > target_faces:
#         new_face_count = max(target_faces, current_faces // 4)
#         #new_face_count = max(target_faces, current_faces - int(current_faces * 0.25))  # Reduce by 25%
#         print("face target: ", new_face_count)

#         print(f"##### Iteration {i}: reducing to {new_face_count} faces #####")
        
#         start_time = time.time()
#         #ms.meshing_decimation_quadric_edge_collapse(targetfacenum=new_face_count, optimalplacement=False, preserveboundary=True, preservetopology=True)
#         ms.meshing_decimation_quadric_edge_collapse(
#             targetfacenum=new_face_count,
#             qualitythr=0.3,  # Increase this value to preserve mesh quality
#             preserveboundary=True,
#             boundaryweight = 1.0,
#             preservenormal = False,
#             preservetopology=True,
#             optimalplacement=False,
#             planarquadric=False,
#             planarweight=0.001,
#             qualityweight=False,
#             autoclean=True,
#             selected=False
#         )
#         end_time = time.time()
        
#         current_faces = ms.current_mesh().face_number()
#         print(f"Number of faces after iteration {i}: {current_faces}")
#         print(f"Time taken for iteration {i}: {end_time - start_time:.2f} seconds")
        
#         ms.save_current_mesh(f"{output_mesh_basepath}_M{i}.obj", save_vertex_normal=False, save_vertex_color=False)
#         print(f"Mesh saved: {output_mesh_basepath}_M{i}.obj")
#         print("")
        
#         i += 1
#         mesh_sequence.append(ms.current_mesh())

#     print("Simplification completed.")
#     return mesh_sequence

# # Usage:
# # input_mesh_path = 'data/dragon.ply'
# # output_mesh_basepath = 'data/dragon_simplified'

# # input_mesh_path = 'dode/dode_tri.obj'
# # output_mesh_basepath = 'dode/dode_simplified'

# # input_mesh_path = 'dragon.ply'
# # output_mesh_basepath = 'prova_decimator/dragon_simplified'
# input_mesh_path = 'thai_d.obj'
# output_mesh_basepath = 'prova_decimator/thai_simpl'

# mesh_sequence = simplify_mesh(input_mesh_path, output_mesh_basepath)

# print(f"Total number of meshes generated: {len(mesh_sequence)}")


import pymeshlab
import time

def truncate_vertex_coords(obj_content, precision=7):
    """
    Truncate vertex coordinates in the .obj content to the specified precision.
    """
    truncated_lines = []
    for line in obj_content.splitlines():
        if line.startswith('v '):
            parts = line.split()
            truncated_parts = [f'{float(coord):.{precision}f}' for coord in parts[1:]]
            truncated_line = f'v {" ".join(truncated_parts)}'
            truncated_lines.append(truncated_line)
        else:
            truncated_lines.append(line)
    return "\n".join(truncated_lines)

def simplify_mesh(input_mesh_path, output_mesh_basepath):
    ms = pymeshlab.MeshSet()
    ms.load_new_mesh(input_mesh_path)

    print(f"Mesh loaded successfully. Initial number of faces: {ms.current_mesh().face_number()}")

    mesh_sequence = [ms.current_mesh()]
    target_faces = 500  # Target number of faces
    current_faces = ms.current_mesh().face_number()

    i = 0
    while current_faces > target_faces:
        new_face_count = max(target_faces, current_faces // 2)
        print("face target: ", new_face_count)

        print(f"##### Iteration {i}: reducing to {new_face_count} faces #####")

        start_time = time.time()
        ms.meshing_decimation_quadric_edge_collapse(
            targetfacenum=new_face_count,
            #qualitythr=0.3,
            preserveboundary=True,
            boundaryweight=1.0,
            preservenormal=False,
            preservetopology=True,
            optimalplacement=False,
            planarquadric=False,
            planarweight=0.001,
            qualityweight=False,
            autoclean=True,
            selected=False
        )
        end_time = time.time()

        current_faces = ms.current_mesh().face_number()
        print(f"Number of faces after iteration {i}: {current_faces}")
        print(f"Time taken for iteration {i}: {end_time - start_time:.2f} seconds")

        # Save the mesh temporarily
        temp_output_path = f"{output_mesh_basepath}_M{i}_temp.obj"
        ms.save_current_mesh(temp_output_path, save_vertex_normal=False, save_vertex_color=False)
        print(f"Mesh saved temporarily: {temp_output_path}")

        # # Read, truncate vertex coordinates, and rewrite the .obj file
        # with open(temp_output_path, 'r') as file:
        #     obj_content = file.read()
        # truncated_content = truncate_vertex_coords(obj_content, precision=7)

        # final_output_path = f"{output_mesh_basepath}_M{i}.obj"
        # with open(final_output_path, 'w') as file:
        #     file.write(truncated_content)
        # print(f"Mesh saved with truncated coordinates: {final_output_path}")
        # print("")

        i += 1
        mesh_sequence.append(ms.current_mesh())

    print("Simplification completed.")
    return mesh_sequence

# Usage:
input_mesh_path = 'bunny.ply'
output_mesh_basepath = 'prova_decimator/bunny_simpl'

mesh_sequence = simplify_mesh(input_mesh_path, output_mesh_basepath)

print(f"Total number of meshes generated: {len(mesh_sequence)}")

