import astropy.units as u
from dataclasses import dataclass
import numpy as np
import pathlib

from spark.io.models import SparkInstrument, SparkSpec, ObservationTarget
import spark.constants as sc

SDSS_WAVE_UNIT = u.AA

@dataclass
class SDSSSpectrum(SparkSpec):

    @classmethod
    def from_fits(cls, hdu, **kwargs):

        name = None
        hdu_close = False
        if isinstance(hdu, (str,pathlib.Path)):
            hdu = pathlib.Path(hdu)
            hdu = fits.open(hdu)
            hdu_close = True
            name = hdu.stem

        if name is None:
            name = kwargs.pop('name', str(hdu[2].data['SPECOBJID'][0]))

        data = hdu[1].data

        inst_data = {}
        inst_data['wave'] = np.power(10, data['loglam'])
        inst_data['wave_unit'] = SDSS_WAVE_UNIT
        inst_data['flux'] = data['flux']
        err = np.sqrt(1 / data['ivar'])
        err[np.isnan(err)] = np.nanmedian(err)
        err[np.isinf(err)] = np.nanmedian(err)
        inst_data['err'] = err
        inst_data['flux_unit'] = sc.FLUX_UNIT_AA

        inst_data['target'] = ObservationTarget(
            name = name,
            ra = hdu[0].header['RA'],
            dec = hdu[0].header['DEC'],
            z = kwargs.pop('z', hdu[2].data['z'][0]),
            instrument = SparkInstrument(name='SDSS'),
        )

        inst_data.update(kwargs)
        if hdu_close: hdu.close()
        return cls(**inst_data)
