import sys
from os import path
#Make the 'src' folder visible to this script 
sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) ) )

from src.vertex_splitter import split_mesh_along_facets
from src.meshdata import MeshData
import networkx as nx
import os
import numpy as np

MESH_2D = MeshData(np.array([[0,0],
						     [1,0],
						     [2,0],
						     [3,0],
						     [0,1],
						     [1,1],
						     [2,1],
					   	     [3,1],
					  	     [0,2],
					  	     [1,2],
					  	     [2,2],
					  	     [3,2]], dtype = float),
				   np.array([[0,1,5],
					  		 [0,5,4],
					  		 [1,2,6],
					  		 [1,6,5],
					  		 [2,3,7],
					  		 [2,7,6],
					  		 [4,5,9],
					  		 [4,9,8],
					  		 [5,6,10],
					  	     [5,10,9],
					  		 [6,7,11],
					  		 [6,11,10]])
					)

MESH_3D = MeshData(np.array([[0, 0, 0],
       						 [1, 0, 0],
       						 [0, 1, 0],
       						 [1, 1, 0],
       						 [0, 0, 1],
       						 [1, 0, 1],
       						 [0, 1, 1],
       						 [1, 1, 1],
       						 [0, 0, 2],
       						 [1, 0, 2],
       						 [0, 1, 2],
       						 [1, 1, 2]], dtype = float),
					np.array([[ 0,  1,  3,  7],
       						  [ 0,  1,  5,  7],
       						  [ 0,  4,  5,  7],
       						  [ 0,  2,  3,  7],
       					      [ 0,  4,  6,  7],
       						  [ 0,  2,  6,  7],
       						  [ 4,  5,  7, 11],
       						  [ 4,  5,  9, 11],
       						  [ 4,  8,  9, 11],
       						  [ 4,  6,  7, 11],
       						  [ 4,  8, 10, 11],
       						  [ 4,  6, 10, 11]]))

class TestCase(object):
	def __init__(self, meshdata, split_facets, num_expected_components):
		self.meshdata = MeshData(meshdata.verts.copy(), 
								 meshdata.elems.copy(),
								 elem_markers = np.arange(len(meshdata.elems)))

		self.meshdata.facets = split_facets
		self.num_expected_components = num_expected_components

	def calculate_splitmesh(self):
		meshcopy = MeshData(self.meshdata.verts.copy(), self.meshdata.elems.copy())
		return split_mesh_along_facets(meshcopy, self.meshdata.facets, remove_disconnected_elems = False)

def get_vertgraph(mesh):
	return nx.Graph(flatlist([zip(e, np.roll(e, 1)) for e in mesh.elems]))

##################################################################
#3D tests
##################################################################

def test_split_mesh_along_facets_3DCheckCorrectNumberOfComponents():
	split_in_two = [[4,5,7], [4,6,7]]

	testcase_in_two = TestCase(MESH_3D,
					  		   split_in_two,
							   2)

	split_mesh = testcase_in_two.calculate_splitmesh()

	vertgraph_split = get_vertgraph(split_mesh)
	assert len(list(nx.connected_components(vertgraph_split))) == testcase_in_two.num_expected_components

def test_split_mesh_along_facets_3DoneFacetShouldProduce3Verts():
	split_triangles = [[0,2,7]]
	testcase = TestCase(MESH_3D,
					  	split_triangles,
						1)

	split_mesh = testcase.calculate_splitmesh()
	assert len(MESH_3D.verts) + 3 == len(split_mesh.verts)

def test_split_mesh_along_facets_3D2intersectingtriangles_shouldNotDisconnectMesh():
	split_triangles = [[0,2,7], [4,7,6]]
	testcase = TestCase(MESH_3D,
					  	split_triangles,
						1)

	split_mesh = testcase.calculate_splitmesh()
	vertgraph_split = get_vertgraph(split_mesh)

	assert len(list(nx.connected_components(vertgraph_split))) == 1

##################################################################
#2D tests
##################################################################
def test_split_mesh_along_facets_2DcheckCorrectNumberOfComponents():
	split_bound_3way = [[1, 5],
						[5, 9],
						[0, 5],
						[5, 6]]
	
	split_bound_H = [[1,5],
					 [5,9],
					 [5,6],
					 [2,6],
					 [6,10]]

	testcase_3way = TestCase(MESH_2D,
							 split_bound_3way,
							 3)
	
	testcase_H = TestCase(MESH_2D,
	 					  split_bound_H,
	 					  4)

	vertgraph_orig = get_vertgraph(MESH_2D)
	assert len(list(nx.connected_components(vertgraph_orig))) == 1

	for testcase in [testcase_3way, testcase_H]:
		splitmesh = testcase.calculate_splitmesh()

		vertgraph_split = get_vertgraph(splitmesh)
		assert len(list(nx.connected_components(vertgraph_split))) == testcase.num_expected_components

def test_split_mesh_along_facets_2DOneEdgeShouldAddTwoVerticies():
	split = [[5, 6]]

	testcase = TestCase(MESH_2D, split, 1)
	
	splitmesh = testcase.calculate_splitmesh()  #split_mesh_along_facets(testcase.copy_mesh(), split)

	#2 extra verticies
	assert len(testcase.meshdata.verts) + 2 == len(splitmesh.verts)

	#Only one connected component
	vertgraph_split = get_vertgraph(splitmesh)
	assert len(list(nx.connected_components(vertgraph_split))) == testcase.num_expected_components

def flatlist(listtoflatten):
	return [item for sublist in listtoflatten for item in sublist]

if __name__ == "__main__":	
	test_split_mesh_along_facets_2DOneEdgeShouldAddTwoVerticies()
	test_split_mesh_along_facets_2DcheckCorrectNumberOfComponents()
	test_split_mesh_along_facets_3DCheckCorrectNumberOfComponents()
	test_split_mesh_along_facets_3DoneFacetShouldProduce3Verts()
	test_split_mesh_along_facets_3D2intersectingtriangles_shouldNotDisconnectMesh()
	print "Tests passed!"