# NOT WORKING
import pymeshlab as ml

def convert_obj_to_ply(input_file, output_file):
    # Load mesh from .obj file
    ms = ml.MeshSet()
    ms.load_new_mesh(input_file)

    # Save mesh to .ply file
    ms.save_current_mesh(output_file, save_vertex_normal=False, save_vertex_color=False)

if __name__ == "__main__":
    #input_file = "data/3holes.obj"  # Replace with the path to your .obj file
    input_file = "Esperimento_1/data/blub_triangulated.obj"
    #output_file = "3holes.ply"  # Desired path for the output .ply file
    output_file = "Esperimento_1/data/blub_triangulated.ply"
    convert_obj_to_ply(input_file, output_file)
