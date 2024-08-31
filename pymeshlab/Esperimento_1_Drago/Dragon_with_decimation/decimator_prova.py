import pymeshlab
import time

def simplify_mesh(input_mesh_path, output_mesh_basepath):
    ms = pymeshlab.MeshSet()
    ms.load_new_mesh(input_mesh_path)

    print(f"Mesh loaded successfully. Initial number of faces: {ms.current_mesh().face_number()}")

    mesh_sequence = [ms.current_mesh()]
    target_faces = 9  # Target number of faces
    current_faces = ms.current_mesh().face_number()

    i = 0
    while current_faces > target_faces:
        new_face_count = max(target_faces, current_faces // 4)
        print(f"##### Iteration {i}: reducing to {new_face_count} faces #####")
        
        start_time = time.time()
        # ms.meshing_decimation_quadric_edge_collapse(targetfacenum=new_face_count, optimalplacement=False)
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
        
        ms.save_current_mesh(f"{output_mesh_basepath}_M{i}.obj", save_vertex_normal=False, save_vertex_color=False)
        print(f"Mesh saved: {output_mesh_basepath}_M{i}.obj")
        print("")
        
        i += 1
        mesh_sequence.append(ms.current_mesh())

    print("Simplification completed.")
    return mesh_sequence

# Usage:
# input_mesh_path = 'data/dragon.ply'
# output_mesh_basepath = 'data/dragon_simplified'

input_mesh_path = 'bunny-saved-ply.ply'
output_mesh_basepath = 'decimator_prova/bunny_simplified'


mesh_sequence = simplify_mesh(input_mesh_path, output_mesh_basepath)

print(f"Total number of meshes generated: {len(mesh_sequence)}")
