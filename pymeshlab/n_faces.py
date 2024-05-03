import pymeshlab

def count_faces(file_path):
    ms = pymeshlab.MeshSet()
    ms.load_new_mesh(file_path)
    mesh = ms.current_mesh()
    num_faces = mesh.face_number()
    print("Number of faces:", num_faces)

if __name__ == "__main__":
    file_path = "bunny.ply"  # Replace with the path to your .ply file
    count_faces(file_path)
