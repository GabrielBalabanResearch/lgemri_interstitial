import argparse

MESHPATH = "example_data/imageslice_mesh.h5"
IMAGEPATH = "example_data/rawimage.nii"
SEGPATH = "example_data/segmentation.nii"


def main(args):
	#TODO load the data and put it back into mm to match the image.
	mesh = H5Format().read(args.input)
	
	imdata = nib.load(args.imagepath)
	seg = nib.load(args.segpath)


	#mesh = facetsplit_mesh(mesh,
	#					   imdata,
#						   seg,
	#					   args.scar_marker,
	#					   args.maxden,
	#					   args.anisotropy,
	#					   args.transpose)
	#mesh.save(args.output)




if __name__ == "__main__":
	parser = argparse.ArgumentParser("Create a network of cracks in a mesh based on LGE-MRI image intensity")
	parser.add_argument("meshpath", default = MESHPATH, help = "Path to input mesh in .h5 format")
	parser.add_argument("imagepath", default = IMAGEPATH, help = "Path to image intensity data in .nii format")
	parser.add_argument("segpath", default = SEGPATH, help = "Path to segmentation data in .nii format")
	parser.add_argument("-output", default = "output_split_mesh.h5", help = "Path to output mesh")
	parser.add_argument("-mark_myo", default = 1, help = "Marker of non-enhanced myocardium in the segmentation (default 1)")
	parser.add_argument("-mark_scar", default = 3, help = "Marker of enhanced myocardium (scar) in the segmentation (default 3)")
	args = parser.parse_args()
	main(args)