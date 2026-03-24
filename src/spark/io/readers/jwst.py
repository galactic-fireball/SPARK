import astropy.units as u
from dataclasses import dataclass
import numpy as np
import pathlib

from spark.io.models import SparkCube
import spark.constants as sc

JWST_WAVE_UNIT = u.um

@dataclass
class JWSTCube(SparkCube):

	@classmethod
	def from_fits(cls, hdu, **kwargs):

		hdu_close = False
		if isinstance(hdu, (str,pathlib.Path)):
			hdu = fits.open(hdu)
			hdu_close = True

		header = hdu['SCI'].header
		cunit = u.Unit(header['CUNIT3'])
		bunit = u.Unit(header['BUNIT'])

		nwave = hdu['SCI'].data.shape[0]
		wave0 = header['CRVAL3'] - (header['CRPIX3'] - 1) * header['CDELT3']
		obs_wave = (wave0 + np.arange(nwave)*header['CDELT3']) * cunit

		cube_spec = hdu['SCI'].data.T * bunit
		cube_err = hdu['ERR'].data.T * bunit
		if ('/sr' in str(bunit)) or ('/ sr' in str(bunit)):
		    pxar = header['PIXAR_SR'] * u.sr
		    cube_spec *= pxar
		    cube_err *= pxar

		cube_spec = cube_spec.to(sc.FLUX_UNIT_UM, equivalencies=u.spectral_density(obs_wave))
		cube_err = cube_err.to(sc.FLUX_UNIT_UM, equivalencies=u.spectral_density(obs_wave))
		cube_err[np.isnan(cube_err)] = np.nanmedian(cube_err)

		inst_data = {}
		inst_data['wave'] = obs_wave
		inst_data['wave_unit'] = cunit
		inst_data['flux'] = cube_spec
		inst_data['err'] = cube_err
		inst_data['flux_unit'] = sc.FLUX_UNIT_UM
		print('WARNING: JWST fits headers do not store redshift data, please supply your own if needed')
		inst_data['z'] = 0.0
		inst_data['ra'] = hdu[0].header['TARG_RA']
		inst_data['dec'] = hdu[0].header['TARG_DEC']

		inst_data.update(kwargs)
		if hdu_close: hdu.close()
		return cls(**inst_data)


	def __post_init__(self):
		super().__post_init__()
		self.wave = self.wave.to(JWST_WAVE_UNIT)
		self.wave_unit = JWST_WAVE_UNIT


@dataclass
class MIRI_MRS_Cube(JWSTCube):
	pass


@dataclass
class NIRSPEC_IFU_Cube(JWSTCube):
	pass
