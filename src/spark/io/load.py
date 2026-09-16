from astropy.io import fits
from astropy.wcs import FITSFixedWarning
import pathlib
import warnings

warnings.simplefilter('ignore', FITSFixedWarning)

from spark.io.readers.jwst import NIRSPEC_IFU_Cube, MIRI_MRS_Cube, DJASpec
from spark.io.readers.sdss import SDSSSpectrum

SPEC_FMTS = {
	'SDSS': SDSSSpectrum,
	'NIRSPEC_IFU': NIRSPEC_IFU_Cube,
	'MIRI_MRS': MIRI_MRS_Cube,
	'DJA': DJASpec,
}

def get_supported_instruments():
	return list(SPEC_FMTS.keys())


def get_fits_fmt(hdu):
	# TODO
	print('get_fits_fmt unimplemented!')
	return None


def load_spectral_product(spec_file, fmt=None, **kwargs):
	if fmt is None:
		hdu = fits.open(spec_file)
		fmt = get_fits_fmt(hdu)
		hdu.close()
	if fmt is None:
		raise Exception('Could not determine fits format')
	if not fmt in SPEC_FMTS:
		raise Exception('Unsupported format: %s'%fmt)

	spec = SPEC_FMTS[fmt].from_fits(spec_file, **kwargs)
	if spec.file is None and isinstance(spec_file, (str,pathlib.Path)):
		spec.file = pathlib.Path(spec_file)

	return spec
