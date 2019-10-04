# What this script does
This repository provides a script that creates a series of cracks in a triangular or tetrahedral mesh that serve as a model for interstitial fibrosis. Each crack is a mesh edge (2D) or facet (3D) whose nodes are doubled.
If an electrical PDE such as the monodomain equation is solved over the mesh than the cracks are effectively
local no-flux boundaries which disrupt the flow of electricity (reference upcoming paper).

The cracks are created with the following probability

h<sub>&theta;</sub>(x) = &theta;<sub>o</sub> x + &theta;<sub>1</sub>x

# Tests
A few unit-tests are present that check the facet and edge splitters against some simple cases. To
run the tests use the command 

`python test/test_vertex_splitter.py`

# How to use it
You will need a medical image and it's segmentation in nifti format, as well as a mesh whose coordinates 
match the world coordinates in the images (e.g. pixels transformed by the affine transformation matrix.) 

![Example images](/images/example_images.png)

One image should contain the raw pixel values, while the other should contain a segmentation of the myocardium and the 
scar with seperate markers. The myocardium should be a complete ring with a blood pool in the middle. All non-myocardial and scar areas should be marked 0. The mesh should be formatted as in 'example_inputs/imageslice_mesh.h5', it should have coordinates, element markers, fibres, and topology. The elements markers should be 0 for
nonscar and scar should be marked as in the segmentation image. The fibres are a vector field giving the direction
of anisotropy. Cracks will be more likely to occur if the mesh entity is aligned with this direction. The topology
is how the coordinates are connected (triangles or tetrahedra).

The script can then be run with the command

`python create_split_network.py`

To see a list of the script command line options write 

`python create_split_network.py -h`

Running the script without any parameters will analyze the data in the 'example_inputs' folder. This should 
create a mesh similar to the one in 'example_output'. Also in the 'example_output' are .xdmf files 
which can be used to visualize the example mesh in paraview. These files can be freely modified if you want to used 
them to visualize your own meshes.  

If you want to modify the create_split_network.py script to work with some other mesh formats. Please see the class H5Format in 'src/meshdata.py', which can serve as template for another mesh format input/output class.

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


# Lisence 
CC-BY 4.0 or later version (https://creativecommons.org/licenses/)