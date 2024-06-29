import pymeshlab
from multiprocessing import Process
import numpy as np

def generate_subdivisions(input_mesh, num_subdivisions):
    # Load the original mesh
    ms = pymeshlab.MeshSet()
    ms.load_new_mesh(input_mesh)
    
    # Apply subdivision repeatedly, saving each result
    for i in range(num_subdivisions):
        # Apply the Loop subdivision filter with correct parameters
        ms.apply_filter('meshing_surface_subdivision_loop', loopweight='Loop', iterations=1)
        # Save the current mesh in a new file with .obj extension
        ms.save_current_mesh(f'data/subdivided_mesh_level_{i+1}.ply')

# Call the function with the path to your mesh and the desired number of subdivisions
generate_subdivisions('data/blub_triangulated.ply', 7)


