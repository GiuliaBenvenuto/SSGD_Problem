import pymeshlab as ml
from multiprocessing import Process
import numpy as np
from scipy.spatial import distance_matrix

### reFine-Coarsen-Match ###

# PyMeshLab script to refine and then coarsen a mesh, while maintaining vertex
# correspondence between the refined and coarsened meshes.
# The user requests a minimum number of faces NF for the fine mesh, and a target
# number NC of faces for the coarse mesh. The first NC vertices of the fine mesh
# are exactly the vertices of the coarse mesh, in the same order.

### Constants ###

def fcmatch(mesh_base_name, fine_faces_min = 100000, coarse_faces = 5000, quiet=False):
	mesh_name = f'{mesh_base_name}.ply'
	mesh_name_f = f'{mesh_base_name}_fine.ply'
	mesh_name_c = f'{mesh_base_name}_coarse.ply'

	# mesh_name = f'{mesh_base_name}.obj'
	# mesh_name_f = f'{mesh_base_name}_fine.obj'
	# mesh_name_c = f'{mesh_base_name}_coarse.obj' 

	### Load mesh ###
	ms = ml.MeshSet()
	ms.load_new_mesh(mesh_name)
	mf = ms.current_mesh()
	idf = ms.current_mesh_id()
	
	print("Initial mesh loaded with", mf.vertex_number(), "vertices.")
	print("Sample initial vertices:", mf.vertex_matrix()[:5])

	### Refine to target ###
	fnum = ms.current_mesh().face_number()
	if not quiet: print(f'Refining mesh {mesh_base_name}...')
	while fnum < fine_faces_min:
		if not quiet: print('\t', fnum, '<', fine_faces_min)
		ms.meshing_surface_subdivision_ls3_loop(
			loopweight = 'Loop',
			iterations = 1,
			threshold = ml.PercentageValue(0)
		)
		fnum = ms.current_mesh().face_number()

	### Geometry-preserving decimation ###
	ms.add_mesh(ms.current_mesh(), mesh_name='coarse', set_as_current=True)
	mc = ms.current_mesh()
	idc = ms.current_mesh_id()
	if not quiet: print(f'Decimating mesh {mesh_base_name}...')
	ms.meshing_decimation_quadric_edge_collapse(
		targetfacenum = 5000,
		qualitythr = .5,
		preservenormal = True,
		optimalplacement = False
	)
	if not quiet: print(f'Flipping edges of mesh {mesh_base_name}...')
	# for i in range(3):
	ms.meshing_edge_flip_by_planar_optimization(
		planartype = 'delaunay',
		iterations = 0,
		pthreshold = 45
	)

	### Match points ###
	if not quiet: print(f'Matching vertices for mesh {mesh_base_name}...')
	gf = mf.vertex_matrix()
	gc = mc.vertex_matrix()
	vfNum = gf.shape[0]
	vcNum = gc.shape[0]
	dm = distance_matrix(gc, gf, p=1)	# shape (nc, nf)
	cor = np.argmin(dm, axis=1)	# shape (nc,)

	### Re-index fine mesh ###
	if not quiet: print(f'Re-indexing mesh {mesh_base_name}...')
	allv = np.arange(vfNum, dtype='int32')
	diff = np.setdiff1d(allv, cor)
	newIndices = np.concatenate((cor, diff), axis=0)
	newIndicesInv = np.zeros_like(newIndices)
	for i,x in enumerate(newIndices):
		newIndicesInv[x] = i
	newGf = gf[newIndices]

	tf = mf.face_matrix().flatten()
	ffNum = mf.face_number()
	for i in range(tf.size):
		tf[i] = newIndicesInv[tf[i]]
	newTf = tf.reshape(ffNum,3)

	newMesh = ml.Mesh(newGf, newTf)

	### Save ###
	if not quiet: print(f'Saving {mesh_base_name}...')
	# ms.set_current_mesh(idf)
	ms.add_mesh(newMesh, set_as_current=True)
	ms.save_current_mesh(
		mesh_name_f,
		save_vertex_normal=False,
		save_vertex_color=False
	)
	ms.set_current_mesh(idc)
	ms.save_current_mesh(
		mesh_name_c,
		save_vertex_normal=False,
		save_vertex_color=False
	)
	if not quiet: print(f'Finished {mesh_base_name}.')