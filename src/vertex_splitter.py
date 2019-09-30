import numpy as np
import numpy_indexed as npi
import networkx as nx
import pandas as pd
import os
from meshtoolz.meshdata import MeshData
from meshtoolz.percolation import renumber_mesh
from local_connectivity import *

import time

def split_mesh_along_facets(mesh,
							split_facets,
							remove_disconnected_elems = True):
	"""
	A caveat for the disconnected elems. The algorithm removes elements
	that are edge-disconnected from the network. Some elements may be accessible
	via a single node and should be included.
	"""

	t0_total = time.time()

	edim = mesh.elems.shape[1]
 	split_facets = np.sort(split_facets, axis = 1)

 	C_elems = get_elem_connectivity(mesh.elems, split_facets)

 	if remove_disconnected_elems:
 		mesh, split_facets, C_elems = remove_disconnected_elems_from_mesh(mesh,
 																		  split_facets,
 																		  C_elems)
 	#Create edge markers to indicate
	facets = mesh.get_facets()
	facet_markers = np.zeros(len(facets))
	facet_markers[npi.indices(facets, split_facets)] = 1
	mesh.facets = facets
	mesh.facet_markers = facet_markers

 	#Make vertex-element data structure
 	split_verticies = np.unique(split_facets)

 	v_elem = pd.DataFrame(np.vstack((mesh.elems.flatten(),
 						  np.arange(len(mesh.elems)*edim) / edim)).transpose(), 
 						  columns = ["vnum", "elemnum"])

	v_elem = v_elem[npi.in_(v_elem["vnum"], split_verticies)]
 	v_elem.loc[v_elem["vnum"].isin(split_verticies)]

 	if edim ==3:
 		connectivity_calc = ConnectivityCalculator2d(mesh,
 													 C_elems,
 													 split_facets)
 	elif edim == 4:
 		connectivity_calc = ConnectivityCalculator3d(mesh,
 													 C_elems,
 													 split_facets)
 		
 	else:
 		raise Exception("Cannot handle {}D mesh".format(edim -1))

 	#Loop over every vertex that is to be split
 	t0_facets = time.time()

 	for vnum, vgroup in v_elem.groupby("vnum"):
 		local_components = connectivity_calc.get_local_components(vgroup)

 		for group in local_components[1:]:
 		 	#Add a new vertex for every connected element group after the first
 		 	mesh.verts = np.vstack((mesh.verts, [mesh.verts[vnum]]))
 		 	new_vnum = len(mesh.verts) - 1

 		 	#Update the group elements with the new vertex
 		 	new_elems = mesh.elems[group]
 		 	new_elems[new_elems == vnum] = new_vnum

 			mesh.elems[group] = new_elems

 	print "Finished edgesplitting mesh"
 	print "Total time = {:.3f}".format(time.time() - t0_total)
 	print "Facet growing time ={:.3f}".format((time.time() - t0_facets))
 	return mesh

def remove_disconnected_elems_from_mesh(mesh, split_facets, C_elems):
	"""Elements which are completely disconnected from the rest of the mesh can safely be removed"""

	components = list(nx.connected_components(C_elems))
	print "Number of components after facet disconnection = ", len(components)

	if len(components) == 1:
		#Everything connected so nothing to be done
		return mesh, split_facets, C_elems

	#Some elems are disconnected by split facets so remove them from the mesh
	lc_num = np.argmax([len(c) for c in components]) #Largest component

	connected_elems = np.array(list(components[lc_num]))
	isolated_elems = np.setdiff1d(range(len(mesh.elems)), connected_elems)
	print "Number elements removed = ", len(isolated_elems)

	assert isolated_elems.any()

	#update split facets
	removed_facets = npi.unique(np.sort(np.sort(np.array([np.roll(mesh.elems[isolated_elems], i, axis = 1).flatten() for i in range(mesh.edim -1)]).transpose(), axis = 1), axis =1))
	split_facets = split_facets[np.logical_not(npi.in_(split_facets, np.array(removed_facets)))]

	#Update the connectivity graph
	C_elems.remove_nodes_from(isolated_elems)
	ediff = np.zeros(len(mesh.elems), dtype = int)
	for enum in isolated_elems:
		ediff[enum:] += 1
	elem_connections = np.array(C_elems.edges())
	elem_connections -= ediff[elem_connections]
	C_elems = nx.Graph(elem_connections.tolist())

	# #Renumber and update mesh
	is_elem_removed = np.zeros(len(mesh.elems))
	is_elem_removed[isolated_elems] = 1

	mesh, vdiff = renumber_mesh(mesh,
	  						    is_elem_removed, 
	  						    return_vdiff  = True)

	#Renumber the split vertex numbers 
	split_facets -= vdiff[split_facets]

	# debugmesh = MeshData(mesh.verts, mesh.elems, facets = split_facets)
	# debugmesh.save("debugmesh.h5")
	# os.system("pviewh5 debugmesh.h5: -elems_db_path debugmesh.h5:facets -xdmf_only -xdmf_file facets.xdmf")
	# os.system("pviewh5 debugmesh.h5:")
	# exit()
	return mesh, split_facets, C_elems