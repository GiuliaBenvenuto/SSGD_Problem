from fcmatch import fcmatch

# https://www.cim.mcgill.ca/~shape/benchMark/#download
# bunny number of faces: 28576
# fcmatch(mesh_base_name='cube', fine_faces_min=50000, coarse_faces=5000)
fcmatch(mesh_base_name='bunny', fine_faces_min=100000, coarse_faces=5000)
