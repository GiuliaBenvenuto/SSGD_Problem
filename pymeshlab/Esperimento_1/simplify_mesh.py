import pymeshlab as pml
import numpy as np

def decimate_mesh(verts, faces, target, optimalplacement=True):
    _ori_vert_shape = verts.shape
    _ori_face_shape = faces.shape
    
    m = pml.Mesh(verts, faces)
    ms = pml.MeshSet()
    ms.add_mesh(m, 'mesh')
    
    # Gradual decimation
    current_faces = faces.shape[0]
    target_faces = int(target)
    
    while current_faces > target_faces:
        try:
            # Calculate intermediate target (reduce by 20% each iteration)
            intermediate_target = max(int(current_faces * 0.8), target_faces)
            
            ms.meshing_decimation_quadric_edge_collapse(
                targetfacenum=intermediate_target,
                optimalplacement=optimalplacement,
                preserveboundary=True,
                preservenormal=True,
                qualitythr=0.3
            )
            
            # Check the result
            m = ms.current_mesh()
            current_faces = m.face_number()
            
            if current_faces == 0:
                print("[WARNING] Decimation resulted in zero faces. Reverting to previous state.")
                break
            
        except Exception as e:
            print(f"[ERROR] An error occurred during decimation: {str(e)}")
            break
    
    
    # Extract mesh
    m = ms.current_mesh()
    verts = m.vertex_matrix()
    faces = m.face_matrix()

    print(f'mesh decimation: {_ori_vert_shape} --> {verts.shape}, {_ori_face_shape} --> {faces.shape}')
    return verts, faces


# Call the function
import trimesh

input_file = 'data/spot/spot_tri.obj'
output_file = 'data/simplified/spot_500f.obj'

# Load the mesh
mesh = trimesh.load_mesh(input_file)

# Get vertices and faces
verts = mesh.vertices
faces = mesh.faces

# Decimate the mesh
new_verts, new_faces = decimate_mesh(verts, faces, target=500)

# Create a new mesh with the decimated data
new_mesh = trimesh.Trimesh(vertices=new_verts, faces=new_faces)

# Save the new mesh
new_mesh.export(output_file)

