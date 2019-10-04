import numpy as np
import numpy_indexed as npi
import h5py

H5_MESHVERTS = "coordinates"
H5_MESHELEMS = "topology"
H5_MESHFACETS = "facets"

H5_VERTMARKERS = "vertex_markers"
H5_FACETMARKERS = "facet_markers"
H5_ELEMMARKERS = "element_markers"

H5_MESHFIBRES = "fibres"
H5_MESH_FIBRE_SHEETS = "fibre_sheets"
H5_MESH_FIBRE_SHEETNORMALS = "fibre_sheetnormals"

#From h5.py File.create_dataset.
H5_COMPRESSION = "gzip"
H5_COMPRESSION_LEVEL = 6

class MeshData(object):
	def __init__(self, 
				 verts,
				 elems,
				 facets = None,
				 vert_markers = None,
				 facet_markers = None,
				 elem_markers = None,
				 fibres = None,
				 fibre_sheets = None,
				 fibre_sheet_normals = None):

		self.verts = verts
		self.facets = facets
		self.elems = elems
		
		self.vert_markers = vert_markers
		self.facet_markers = facet_markers
		self.elem_markers = elem_markers
		self.fibres = fibres
		self.fibre_sheets = fibre_sheets
		self.fibre_sheet_normals = fibre_sheet_normals

		self.edim = elems.shape[1]
		self.vdim = verts.shape[1]

	def get_facets(self):
		facets = np.vstack([np.roll(self.elems, i, axis = 1).flatten() for i in range(self.edim - 1)]).transpose()
		return npi.unique(np.sort(facets, axis = 1))

	#def save(self, path, mode = "a"):
	# 	meshformat_str = path.split(".")[-1]
	# 	from formats import FORMATS
	# 	FORMATS[meshformat_str]().write(path, self, mode = mode)

	# def get_boundary(self):
	# 	import dolfin as df
	# 	from meshtoolz.formats.dolfin_format import DolfinFormat
	# 	dolfin_mesh = DolfinFormat().make_mesh(self)

	# 	bmesh = df.BoundaryMesh(dolfin_mesh, "exterior")
	# 	vb_v = np.array(bmesh.entity_map(0), dtype = int) #boundary vertex to interior vertex
	# 	return MeshData(bmesh.coordinates(),
	# 					bmesh.cells(),
	# 					boundary_map = vb_v)

	def __str__(self):
		return "{}D MeshData with {} verticies and {} elements".format(self.vdim, len(self.verts), len(self.elems))


class H5Format(object):
	def parse_path(self, path):
		if ":" in path:
			path_in_db = path.split(":")[-1]
			db_path = path.rstrip(":" + path_in_db)
		else:
			db_path, path_in_db = path, ""

		if len(db_path.split(".")) == 1: db_path += ".hdf5"
		return db_path, path_in_db

	def read(self, path):
		db_path, path_in_db = self.parse_path(path)

		with h5py.File(db_path, "r") as f:
			if path_in_db:
				f = f[path_in_db]

			verts = np.array(f[H5_MESHVERTS])
			elems = np.array(f[H5_MESHELEMS])

			if H5_MESHFIBRES in f.keys():
				fibres = np.array(f[H5_MESHFIBRES])
			else: fibres = None

			if H5_MESH_FIBRE_SHEETS in f.keys():
				fibre_sheets = np.array(f[H5_MESH_FIBRE_SHEETS])
			else: fibre_sheets = None

			if H5_MESH_FIBRE_SHEETNORMALS in f.keys():
				fibre_sheet_normals = np.array(f[H5_MESH_FIBRE_SHEETNORMALS])
			else: fibre_sheet_normals = None

			if H5_ELEMMARKERS in f.keys():
				elem_markers = np.array(f[H5_ELEMMARKERS], dtype = int)
			else: elem_markers = None

			if H5_FACETMARKERS in f.keys():
				facet_markers = np.array(f[H5_FACETMARKERS], dtype = int)
			else: facet_markers = None

			if H5_MESHFACETS in f.keys():
				facets = np.array(f[H5_MESHFACETS])
			else: facets = None

		return MeshData(verts,
					    elems, 
						elem_markers = elem_markers,
						facets = facets,
						facet_markers = facet_markers,
						fibres = fibres,
						fibre_sheets = fibre_sheets,
						fibre_sheet_normals = fibre_sheet_normals)

	def h5_save_dataset(self, f, path, dataset):
		if f.get(path): del f[path]
		f.create_dataset(path,
						 data = dataset,
						 compression = H5_COMPRESSION,
						 compression_opts = H5_COMPRESSION_LEVEL)

	def write(self, path, mesh_data, mode = "a"):
		db_path, path_in_db = self.parse_path(path)
		with h5py.File(db_path, mode) as f_db:
			#Mesh
			self.h5_save_dataset(f_db,
							path_in_db + "/" + H5_MESHVERTS,
							mesh_data.verts)

			self.h5_save_dataset(f_db,
							path_in_db + "/" + H5_MESHELEMS,
							mesh_data.elems)

			if mesh_data.facets is not None:
				self.h5_save_dataset(f_db,
								path_in_db + "/" + H5_MESHFACETS,
								mesh_data.facets)
			#Markers
			if mesh_data.vert_markers is not None:
				self.h5_save_dataset(f_db,
								path_in_db + "/" + H5_VERTMARKERS,
								mesh_data.vert_markers)
			
			if mesh_data.facet_markers is not None:
				self.h5_save_dataset(f_db,
								path_in_db + "/" + H5_FACETMARKERS,
								mesh_data.facet_markers)
			
			if mesh_data.elem_markers is not None:
				self.h5_save_dataset(f_db,
								path_in_db + "/" + H5_ELEMMARKERS,
								mesh_data.elem_markers)

			#Other data
			if mesh_data.fibres is not None:
				#Pad 2-D fibres for paraview.
				if mesh_data.verts.shape[1] == 2 and mesh_data.fibres.shape[1] == 2:
					mesh_data.fibres = np.hstack([mesh_data.fibres, np.zeros((mesh_data.fibres.shape[0], 1))])
				self.h5_save_dataset(f_db,
								path_in_db + "/" + H5_MESHFIBRES,
								mesh_data.fibres)

			if mesh_data.fibre_sheets is not None:
				self.h5_save_dataset(f_db,
								path_in_db + "/" + H5_MESH_FIBRE_SHEETS,
								mesh_data.fibre_sheets)

			if mesh_data.fibre_sheet_normals is not None:
				self.h5_save_dataset(f_db,
								path_in_db + "/" + H5_MESH_FIBRE_SHEETNORMALS,
								mesh_data.fibre_sheet_normals)


		print "created " + db_path