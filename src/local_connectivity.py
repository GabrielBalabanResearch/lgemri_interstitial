import numpy as np
import numpy_indexed as npi
import networkx as nx
import pandas as pd

class ConnectivityCalculator2d(object):
	def __init__(self, mesh, global_elem_connectivity, split_facets):
		self.mesh = mesh
 		self.orig_elems = mesh.elems.copy()
 		self.global_elem_connectivity = global_elem_connectivity
 		self.split_facets = split_facets

	def get_local_components(self, vgroup):
		local_connectivity = nx.subgraph(self.global_elem_connectivity,
										 list(vgroup["elemnum"]))

		local_components = [list(eg) for eg in nx.connected_components(local_connectivity)]
		if len(local_components) == 1:
			local_components = self.extend_splitedge(vgroup, local_connectivity)
 		return local_components

	def extend_splitedge(self, vgroup, local_connectivity):
		#The local elements are not disconnected by split edges
		#This means there is only one split edge. Extend this edge to create a disconnected graph
		
		#Step 1: Locate the split edge
		verts = self.mesh.verts
		local_elemconnect = np.sort(np.array(local_connectivity.edges()), axis = 1)
		split_edge = self.split_facets[np.where(np.sum(self.split_facets == vgroup["vnum"].iloc[0], axis = 1))[0]][0]

		#Step 2: Locate the edge that is closest to 180 degrees from the current edge and remove it 
		#from the local connectivity graph

		local_edges = np.array([np.intersect1d(e1, e2) for e1, e2 in zip(self.orig_elems[local_elemconnect[:,0]],
																		 self.orig_elems[local_elemconnect[:,1]])])
		
		local_vecs = verts[local_edges[:,1]] - verts[local_edges[:,0]]
		split_vec = verts[split_edge[0]] - verts[split_edge[1]]
		
		thetas = np.array([self.calc_theta(localvec, split_vec) for localvec in local_vecs])

		remove_edge = local_elemconnect[np.argmin((np.array(thetas) - np.pi)**2)].tolist()
		local_connectivity.remove_edge(*remove_edge)
		return [list(eg) for eg in nx.connected_components(local_connectivity)]

	def calc_theta(self, a, b):
		theta = np.arctan2(np.linalg.det(np.array([a, b])), a.dot(b))
		if theta < 0: theta += 2*np.pi
		return theta

class ConnectivityCalculator3d(object):
	def __init__(self, mesh, global_elem_connectivity, split_facets):
		self.mesh = mesh
 		self.orig_elems = mesh.elems.copy()
 		self.global_elem_connectivity = global_elem_connectivity
 		self.split_facets = split_facets

	def get_local_components(self, vgroup):
		local_connectivity = nx.subgraph(self.global_elem_connectivity,
										 list(vgroup["elemnum"]))

		#If already disconnected then nothing need be done
		if len(local_connectivity.edges()) == 0:
			return [[n] for n in local_connectivity.nodes()]

		#If connections present make sure elements on opposite side of a split edge
		#are in separate components
		local_components = [list(eg) for eg in nx.connected_components(local_connectivity)]
		local_splitplanes = self.get_local_splitplanes(vgroup, local_connectivity)

		#Loop through every split plane, if the two elements are not disconnected then create a disconnection
		for split_plane in local_splitplanes:
			if self.are_elems_connected(split_plane[0], split_plane[1], local_components):
				local_connectivity, local_components = self.extend_splitplane(vgroup, 
																			  local_connectivity,
																			  split_plane)	
		#Safety check
		assert True not in [self.are_elems_connected(split_plane[0],
													  split_plane[1],
													  local_components) for split_plane in local_splitplanes]
 		return local_components

	def extend_splitplane(self,
						  vgroup,
						  local_connectivity,
						  split_plane):

		verts = self.mesh.verts

		#Todo: store all the normals beforehand
		splitplane_normal = self.get_plane_normal(split_plane[0],
												  split_plane[1])
		
		candidate_normals = np.array([self.get_plane_normal(e[0], e[1]) for e in local_connectivity.edges()])
		dot_products = np.abs(np.dot(splitplane_normal, candidate_normals.transpose()))

	 	connected = True
	 	
	 	current_planes = np.arange(len(dot_products))
	 	while connected:
	 		i_remove = dot_products[current_planes].argmax()
			current_planes = np.delete(current_planes, i_remove)
			
			proposed_connectivity = nx.Graph(np.array(local_connectivity.edges())[current_planes].tolist())
			proposed_connectivity.add_nodes_from(local_connectivity)

			proposed_components = [list(eg) for eg in nx.connected_components(proposed_connectivity)]
			
			connected = self.are_elems_connected(split_plane[0],
												 split_plane[1],
												 proposed_components)

		return proposed_connectivity, proposed_components

 	def are_elems_connected(self, enum0, enum1, local_components):
 		loc0 = int(np.where([enum0 in component for component in local_components])[0])
 		loc1 = int(np.where([enum1 in component for component in local_components])[0])
 		return loc0 == loc1

 	def get_local_splitplanes(self, vgroup, local_connectivity):
 		all_local_connectivity = get_elem_connectivity(self.orig_elems[vgroup["elemnum"]], None)
		all_elemconnect = np.sort(np.array(vgroup["elemnum"])[np.array(all_local_connectivity.edges())], axis =1)

		local_elemconnect = np.sort(np.array(local_connectivity.edges()), axis = 1)
		return npi.difference(all_elemconnect, local_elemconnect)

	def get_plane_normal(self, e0, e1):
		plane = np.intersect1d(self.orig_elems[e0],
							   self.orig_elems[e1])
		N = np.cross(self.mesh.verts[plane][0] - self.mesh.verts[plane][1],
			            self.mesh.verts[plane][0] - self.mesh.verts[plane][2])
		return N/np.linalg.norm(N)

def get_elem_connectivity(elems, split_facets):
	edim = elems.shape[1]

 	flatelems = np.array([np.roll(elems, i, axis = 1).flatten() for i in range(edim - 1)]).transpose()

 	enums = np.arange(flatelems.shape[0]) / edim
	facet_data = np.hstack((flatelems, enums[:, np.newaxis]))

	#Sort each element internally so lowest vnum is first
	facet_data[:,:edim -1] =  np.sort(facet_data[:,:edim -1], axis = 1) 
	pkeys = ["p" + str(i) for i in range(edim -1)]

	facet_data = pd.DataFrame(facet_data, columns = pkeys + ["elemnum"])

	#Sort so that duplicate facets are beside one another
	facet_data = facet_data.sort_values(by = pkeys).reset_index(drop = True)

	facets_1 = facet_data[facet_data.duplicated(subset = pkeys, keep = "last")]

	#Remove facets that are to be split
	if split_facets is not None and split_facets.any():
		facets_1 = facets_1[np.logical_not(npi.in_(np.array(facets_1[pkeys]),
												   split_facets,
												   axis = 0))]

	#Get matching elements, due to sorting they should be next to one another
	facets_2 = facet_data.iloc[facets_1.index + 1]

	C_elems = nx.Graph(np.array([facets_1["elemnum"].tolist(),
								 facets_2["elemnum"].tolist()]).transpose().tolist())

	#Make sure all elements are included
	C_elems.add_nodes_from(np.arange(len(elems)))
	return C_elems