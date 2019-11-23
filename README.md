# What this script does
This repository provides a script that creates a series of cracks in a triangular or tetrahedral mesh that serve as a model for interstitial fibrosis. Each crack is a mesh edge (2D) or facet (3D) whose nodes are doubled.
If an electrical PDE such as the monodomain equation is solved over the mesh then the cracks are effectively
local no-flux boundaries which disrupt the flow of electricity (reference upcoming paper).

The cracks are created with the following probability

p = &rho; I cos<sup>&alpha;</sup>(&theta;)

where &rho; is the global crack density, I the local normalized image intensity, and &alpha;, the 
anisotropy parameter. You can change the parameteres &rho; and &alpha; yourself to get different patterns of cracks. The quantity I is determined from the input images, and is the max-min normalized intensity in the areas marked as scars. 

# Tests
A few unit-tests are present that check the facet and edge splitters against some simple cases. To
run the tests use the command 

`python test/test_vertex_splitter.py`

# How to use it
You will need a medical image and it's segmentation in nifti format, as well as a mesh whose coordinates 
match the world coordinates in the images (e.g. pixels locations transformed by the affine transformation matrix.) 

![Example images](/images/example_images.png)

The segmentation should mark the myocardium and the scar with seperate markers. All areas that are not myocardium or scar should be marked 0. The mesh should be formatted as in 'example_inputs/imageslice_mesh.h5': it should have coordinates, element markers, fibres, and topology. The elements markers should be 0 for nonscar and scar should be marked as in the segmentation image. The fibres are a vector field giving the direction of anisotropy. Cracks will be more likely to occur if the mesh entity is aligned with this direction. The topology is how the coordinates are connected (triangles or tetrahedra).

The script can then be run with the command

`python create_split_network.py`

To see a list of the script command line options write 

`python create_split_network.py -h`

Running the script without any parameters will analyze the data in the 'example_inputs' folder. This should 
create a mesh similar to the one in 'example_output'. Also in the 'example_output' are .xdmf files 
which can be used to visualize the example mesh in paraview. These files can be freely modified if you want to use
them to visualize your own meshes.  

If you want to modify the script 'create_split_network.py' to work with other mesh formats, please see the class 'H5Format' in 'src/meshdata.py', which can serve as a template for other mesh input/output classes.

# Known compatible dependencies

* python 2.7.12
* scipy 1.1.0
* numpy 1.15.2
* nibabel 2.1.0
* cv2 2.4.9.1
* pandas 0.24.1
* h5py 2.6.0
* networkx 1.10
* paraview 5.3.0 (for xdmf file visualization)

# Citations
The original crack generation method was developed in 

Costa, Caroline Mendonca, et al. "An efficient finite element approach for modeling fibrotic clefts in the heart." IEEE Transactions on Biomedical Engineering 61.3 (2013): 900-910.

For the extension with topological analysis please cite

(upcoming paper)


# Lisence 
CC-BY 4.0 or later version (https://creativecommons.org/licenses/)