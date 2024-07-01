import pymeshlab as ml

def convert_ply_to_obj(input_file, output_file):
    # Load mesh from .ply file
    ms = ml.MeshSet()
    ms.load_new_mesh(input_file)

    # Save mesh to .obj file
    ms.save_current_mesh(output_file, save_vertex_normal=True, save_vertex_color=True)

if __name__ == "__main__":
    # Replace with the path to your .ply file
    input_file = "Esperimento_1/data/bob_tri_subdiv_2.ply"
    # Desired path for the output .obj file
    output_file = "Esperimento_1/data/bob_tri_subdiv_2.obj"
    convert_ply_to_obj(input_file, output_file)
