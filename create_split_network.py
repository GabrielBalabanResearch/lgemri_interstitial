import argparse
import nibabel as nib
from src.meshdata import H5Format
from src.split_facet_network import facetsplit_mesh

MESHPATH = "example_data/imageslice_mesh.h5"
IMAGEPATH = "example_data/rawimage.nii"
SEGPATH = "example_data/segmentation.nii"


def main(args):
	mesh = H5Format().read(args.meshpath)
	
	imdata = nib.load(args.imagepath)
	seg = nib.load(args.segpath)

	mesh = facetsplit_mesh(mesh,
						   imdata,
						   seg,
						   args.mark_scar,
						   args.maxden,
						   args.anisotropy,
						   transpose = True)
	
	H5Format().write(args.output,
					 mesh)

if __name__ == "__main__":
	parser = argparse.ArgumentParser("Create a network of cracks (interstitial fibrosis) in a mesh based on LGE-MRI image intensity")
	meshgroup = parser.add_argument_group("Meshes")

	meshgroup.add_argument("-meshpath", default = MESHPATH, help = "Path to input mesh in .h5 format")
	meshgroup.add_argument("-output", default = "output_split_mesh.h5", help = "Path to output mesh")
	
	imgroup = parser.add_argument_group("Images")
	imgroup.add_argument("-imagepath", default = IMAGEPATH, help = "Path to image intensity data in .nii format")
	imgroup.add_argument("-segpath", default = SEGPATH, help = "Path to segmentation data in .nii format")	
	imgroup.add_argument("-mark_scar", default = 3, help = "Marker of enhanced myocardium (scar) in the segmentation (default 3)")

	fibrosisgroup = parser.add_argument_group('Fibrosis Probability Parameters')

	fibrosisgroup.add_argument("-anisotropy",
								 default = 4.0,
								 help = "Exponent in cos(theta) relation for probability of interstitial fibrosis, where theta is fiber-mesh entity angle (default 4.0).")
	fibrosisgroup.add_argument("-maxden", default = 1.0, help = "Maximum global scar density parameters (default 1.0).")
	args = parser.parse_args()
	main(args)