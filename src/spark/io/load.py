from astropy.io import fits
import pathlib


from spark.io.readers.jwst import NIRSPEC_IFU_Cube, MIRI_MRS_Cube


CUBE_FMTS = {
	'NIRSPEC_IFU': NIRSPEC_IFU_Cube,
	'MIRI_MRS': MIRI_MRS_Cube,
}


def get_fits_fmt(hdu):
	# TODO
	print('get_fits_fmt unimplemented!')
	return None


def load_cube(cube_file, fmt=None, **kwargs):
	hdu = fits.open(cube_file)

	if fmt is None:
		fmt = get_fits_fmt(hdu)
	if fmt is None:
		raise Exception('Could not determine fits format')
	if not fmt in CUBE_FMTS:
		raise Exception('Unsupported format: %s'%fmt)

	cube = CUBE_FMTS[fmt].from_fits(hdu, **kwargs)
	if cube.file is None:
		cube.file = pathlib.Path(cube_file)

	hdu.close()
	return cube
