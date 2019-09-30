import numpy as np
import numpy_indexed as npi

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
