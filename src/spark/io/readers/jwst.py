from astropy.io import fits
import astropy.units as u
from astropy.wcs import WCS
from dataclasses import dataclass
import numpy as np
import pathlib

from spark.io.models import SparkInstrument, SparkCube, ObservationTarget
import spark.constants as sc

JWST_WAVE_UNIT = u.um

@dataclass
class JWSTCube(SparkCube):

    pix_scale: u.Unit = 0.0

    @classmethod
    def from_fits(cls, hdu, **kwargs):

        name = None
        inst_data = {}
        hdu_close = False
        if isinstance(hdu, (str,pathlib.Path)):
            hdu = pathlib.Path(hdu)
            inst_data['file'] = hdu
            hdul = fits.open(hdu)
            hdu_close = True
            name = hdu.stem
            hdu = hdul

        header = hdu['SCI'].header
        cunit = u.Unit(header['CUNIT3'])
        bunit = u.Unit(header['BUNIT'])

        nwave = hdu['SCI'].data.shape[0]
        wave0 = header['CRVAL3'] - (header['CRPIX3'] - 1) * header['CDELT3']
        obs_wave = (wave0 + np.arange(nwave)*header['CDELT3']) * cunit

        cube_flux = hdu['SCI'].data * bunit
        cube_err = hdu['ERR'].data * bunit
        if ('/sr' in str(bunit)) or ('/ sr' in str(bunit)):
            pxar = header['PIXAR_SR'] * u.sr
            cube_flux *= pxar
            cube_err *= pxar

        cube_flux = cube_flux.T.to(sc.FLUX_UNIT_UM, equivalencies=u.spectral_density(obs_wave)).T
        cube_err = cube_err.T.to(sc.FLUX_UNIT_UM, equivalencies=u.spectral_density(obs_wave)).T
        cube_err[np.isnan(cube_err)] = np.nanmedian(cube_err)

        inst_data['wcs'] = WCS(hdu['SCI'])
        inst_data['wave'] = obs_wave
        inst_data['wave_unit'] = cunit
        inst_data['flux'] = cube_flux
        inst_data['err'] = cube_err
        inst_data['flux_unit'] = sc.FLUX_UNIT_UM
        # print('WARNING: JWST fits headers do not store redshift data, please supply your own if needed')

        inst_data['target'] = ObservationTarget(
            name = kwargs.pop('name', hdu[0].header['TARGNAME']),
            ra = hdu[0].header['TARG_RA'],
            dec = hdu[0].header['TARG_DEC'],
            z = kwargs.pop('z', 0.0),
            instrument = SparkInstrument(name = hdu[0].header['INSTRUME']),
        )

        inst_data['pix_scale'] = (header['CDELT1'] * u.Unit(header['CUNIT1'])).to(u.arcsec)

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
