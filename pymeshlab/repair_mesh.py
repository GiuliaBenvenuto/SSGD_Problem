import pymeshlab

def repair_mesh(mesh_input, mesh_output):
    # Load the mesh
    ms = pymeshlab.MeshSet()
    ms.load_new_mesh(mesh_input)
    
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
    ms.save_current_mesh(mesh_output)
    print(f"Repaired mesh saved as {mesh_output}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python repair_holes.py [input_mesh.ply] [output_mesh.ply]")
    else:
        input_mesh = sys.argv[1]
        output_mesh = sys.argv[2]
        repair_mesh(input_mesh, output_mesh)
