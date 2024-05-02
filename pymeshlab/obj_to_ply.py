import pymeshlab as ml

def convert_obj_to_ply(input_file, output_file):
    # Load mesh from .obj file
    ms = ml.MeshSet()
    ms.load_new_mesh(input_file)

    # Save mesh to .ply file
    ms.save_current_mesh(output_file, save_vertex_normal=False, save_vertex_color=False)

if __name__ == "__main__":
    input_file = "data/cube.obj"  # Replace with the path to your .obj file
    output_file = "cube.ply"  # Desired path for the output .ply file
    convert_obj_to_ply(input_file, output_file)
