from vertex_splitter import split_mesh_along_facets
from scipy.ndimage.filters import convolve
import numpy as np
import numpy_indexed as npi
import time
import os
import pandas as pd

MARK_BACKGROUND = 0

def facetsplit_mesh(mesh,
				    imdata,
				    seg,
				    scar_marker,
				    maxden,
				    anisotropy,
				    transpose = False,
				    unit_conversion = 1.0):
	"""
	Create splits in the mesh 
	according to the image intensities
	in imdata where the scars are marked
	by the scar_markers in seg.
	"""

	#Make sure elements are int
	mesh.elems = mesh.elems.astype(int)

	regist = ImageMeshLgeRegistration(mesh,
									  imdata,
									  seg,
									  scar_marker,
									  transpose,
									  unit_conversion)

	facets = mesh.get_facets()
	facet_midpoints = np.mean(mesh.verts[facets], axis = 1)
	elem_midpoints = np.mean(mesh.verts[mesh.elems], axis = 1)
	facet_intensities = regist.register(facet_midpoints)
	elem_intensities = regist.register(elem_midpoints)

	assert (facet_intensities !=0).any()
	assert (elem_intensities !=0).any()

	lge_facets, lge_facet_elems = get_lge_facet_elements(mesh, scar_marker)

	local_global_facetmap = npi.indices(facets, lge_facets)

	#Map intesity data into facets
	lge_facet_intense = facet_intensities[local_global_facetmap]
	lge_facet_intense = relative_intensity(lge_facet_intense, lge_facet_intense.min())
	lge_facet_fibres = mean_direction(mesh.fibres[lge_facet_elems])

	if mesh.edim == 3:
		#Cosine between fibre and edge
		costheta = get_fibre_edge_cosine(mesh.verts[lge_facets], 
										 lge_facet_fibres[:,:2])

	elif mesh.edim == 4:
		#Cosine between sheet normal and facet normal
		lge_facet_fibre_sheets = mean_direction(mesh.fibre_sheets[lge_facet_elems])
		lge_facet_sheet_normals = np.cross(lge_facet_fibres, lge_facet_fibre_sheets)

		facet_normals = np.cross(mesh.verts[lge_facets][:,0] - mesh.verts[lge_facets][:,1],
								 mesh.verts[lge_facets][:,0] - mesh.verts[lge_facets][:,2])

		facet_normals /= np.linalg.norm(facet_normals, axis = 1)[:,np.newaxis]
		costheta = np.abs(np.sum(lge_facet_sheet_normals*facet_normals, axis = 1))

	#--------------------------------------------
	#The probability formula for the edge splitting
	p = maxden*(costheta**anisotropy)*lge_facet_intense
	#--------------------------------------------

	assert not np.isnan(p).any()

	#Determine which facets are to be split
	is_split_facet = p >= np.random.rand(len(p))
	#from IPython import embed; embed()
	split_facets = lge_facets[is_split_facet]
	return split_mesh_along_facets(mesh, split_facets)

def relative_intensity(intensities, baseline):
	imax = max(intensities)
	intense = (intensities - baseline)/float(imax - baseline)
	return intense

def mean_direction(x):
	mean = np.sum(x, axis = 1)
	mean /= np.linalg.norm(mean, axis = 1)[:, np.newaxis]
	return mean

def get_fibre_edge_cosine(edges, fibres):
	d = edges[:,0,:] - edges[:,1,:]
	d = d/np.expand_dims(np.linalg.norm(d, axis = 1), axis = 1)
	
	#f = edges/np.expand_dims(np.linalg.norm(edges, axis = 1), axis = 1)

	ctheta = np.abs(np.sum(d*fibres, axis = 1)) #np.abs(f[:,0]*d[:,0] + f[:,1]*d[:,1])
	ctheta[ctheta > 1] = 1.0
	return ctheta

def get_lge_facet_elements(mesh, scar_marker):
	edim = mesh.elems.shape[1]

 	facet_data = np.array([np.roll(mesh.elems, i, axis = 1).flatten() for i in range(edim - 1)]).transpose().astype(int)
 
	enums = np.arange(facet_data.shape[0]) / edim
	#elem_markers = mesh.elem_markers[enums]

	facet_data = np.hstack((facet_data, enums[:, np.newaxis]))

	#Sort each element internally so lowest vnum is first
	facet_data[:,:edim -1] =  np.sort(facet_data[:,:edim -1], axis = 1) 
	pkeys = ["p" + str(i) for i in range(edim -1)]

	facet_data = pd.DataFrame(facet_data, columns = pkeys + ["elemnum"])
	
	#Get rid of facets outside of LGE
	facet_data = facet_data.loc[mesh.elem_markers[facet_data["elemnum"]] == scar_marker]

	#Sort so that duplicate facets are beside one another
	facet_data = facet_data.sort_values(by = pkeys).reset_index(drop = True)

	#Get the interior facets
	facets_1 = facet_data[facet_data.duplicated(subset = pkeys, keep = "last")]
	facets_2 = facet_data.iloc[facets_1.index + 1]

	lge_facets = np.array(facets_1[pkeys])
	lge_elems = np.vstack((np.array(facets_1["elemnum"]), facets_2["elemnum"])).transpose()
	return lge_facets, lge_elems

class ImageMeshLgeRegistration(object):
	def __init__(self, mesh, image, segmentation, scar_marker, transpose, unit_conversion):
		self.d = mesh.vdim
	
		self.imdata = np.squeeze(image.get_data())

		#Get rid of negative intensities
		self.imdata += min(0, self.imdata.min()) 
		self.segdata = np.squeeze(segmentation.get_data()) #squeeze gets rid of uncessary dimensions

		if transpose:
			self.imdata = self.imdata.transpose()
			self.segdata = self.segdata.transpose()

		#Check that image and mesh dimensions match
	
		assert len(self.segdata.shape) == self.d
		assert len(self.imdata.shape) == self.d

		self.lge_intensities = np.zeros(self.imdata.shape)
		self.lge_intensities[self.segdata == scar_marker] = self.imdata[self.segdata == scar_marker]
		self.lge_intensities = self.extrapolate_image(self.lge_intensities)
	
		assert scar_marker in np.unique(mesh.elem_markers)
		self.lge_elemnums = np.where(mesh.elem_markers == scar_marker)[0]
		self.lge_centroids = np.average(mesh.verts[mesh.elems[self.lge_elemnums]], axis = 1)

		if self.d == 2:
			self.im2coord = Im2MeshCoord2d(segmentation,
										   unit_conversion = unit_conversion)
		elif self.d == 3:
			self.im2coord = Im2MeshCoord3d(segmentation,
			 							   unit_conversion = unit_conversion)

	def register(self, points):
		points_pixel = self.im2coord.inverse(points, round_to_pixel = True)
		voxnums = np.array([points_pixel[:, i] for i in range(self.d)])
		
		if self.d ==2:
			return self.lge_intensities[voxnums[0], voxnums[1]]
		elif self.d == 3:
			return self.lge_intensities[voxnums[0], voxnums[1], voxnums[2]]

	def extrapolate_image(self, image):
		KERNEL_SIZE = 7

		bin_image = np.zeros(image.shape)
		bin_image[image  > MARK_BACKGROUND] = 1

		#This determines the extent of the extrapolation
		kernel = np.ones([KERNEL_SIZE]*len(bin_image.shape))
	
		num_neighbours = convolve(bin_image, kernel, mode = "constant", cval = 0)

		num_neighbours[num_neighbours == 0] = 1
		extrapolated = convolve(image, kernel)/num_neighbours #Mean of nonzero neighbour values

		#Put original intensity values back into extrapolated image
		extrapolated[image != MARK_BACKGROUND] = image[image != MARK_BACKGROUND]
		return extrapolated

class Im2MeshCoord2d(object):
	"""
	Maps a 2-D image onto physical mesh coordinates 
	so that the mesh and image look the same when
	viewed in their standard way.
	"""

	def __init__ (self, nii_image, unit_conversion = 1.0):
		self.sx, self.sy = np.array(nii_image.get_header().get_zooms()[:2])*unit_conversion
		imdata = nii_image.get_data()

		self.xmax = np.where(imdata >0)[1].max()
		self.xmin = np.where(imdata >0)[1].min()
		self.ymin = np.where(imdata >0)[0].min()
		self.ymax = np.where(imdata >0)[0].max()

	def inverse(self, coords, round_to_pixel = True):
		X,Y = coords[:,0], coords[:,1]
		x = X/self.sx + self.xmin
		y = self.ymax - Y/self.sy
		#x = self.xmax - X/self.sx
		#y = Y/self.sy + self.ymin
		if round_to_pixel:
			x,y = x.astype(int), y.astype(int)
		return np.vstack((x,y)).transpose()

class Im2MeshCoord3d(object):
	def __init__ (self, nii_image, unit_conversion = 1.0):
		self.sx, self.sy, self.sz = np.array(nii_image.get_header().get_zooms()[:3])*unit_conversion

	def inverse(self, coords, round_to_pixel = True):
		scale = np.array([self.sx, self.sy, self.sz])

		imcoords = coords /scale
		if round_to_pixel:
			imcoords = np.around(imcoords).astype(int)
		return imcoords