from astropy.io import fits
import astropy.units as u
from astropy.wcs import WCS
from dataclasses import dataclass
import numpy as np
import pathlib

from spark.io.models import SparkInstrument, SparkCube, ObservationTarget

@dataclass
class MUSECube(SparkCube):

    pix_scale: u.Unit = 0.0
    spec_res: float = 0.0

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

        primary = hdu['PRIMARY'].header
        header = hdu['DATA'].header

        cunit = u.Unit(header['CUNIT3'])
        bunit = u.Unit(header['BUNIT'])

        obs_wave = np.array(header['CRVAL3'] + header['CD3_3']*np.arange(header['NAXIS3'])) * cunit

        cube_flux = hdu['DATA'].data * bunit
        cube_err = np.sqrt(hdu['STAT'].data) * bunit

        inst_data['wcs'] = WCS(hdu['DATA'])
        inst_data['wave'] = obs_wave
        inst_data['wave_unit'] = cunit
        inst_data['flux'] = cube_flux
        inst_data['err'] = cube_err
        inst_data['flux_unit'] = bunit

        inst_data['target'] = ObservationTarget(
            name = kwargs.pop('name', primary['ESO OBS TARG NAME']),
            ra = primary['RA'],
            dec = primary['DEC'],
            z = kwargs.pop('z', 0.0),
            instrument = SparkInstrument(name = primary['INSTRUME']),
        )

        inst_data['pix_scale'] = (np.abs(header['CD1_1'])*u.Unit(header['CUNIT1'])).to(u.arcsec)
        inst_data['spec_res'] = primary['SPEC_RES']

        inst_data.update(kwargs)
        if hdu_close: hdu.close()
        return cls(**inst_data)

