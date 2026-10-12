from astropy.io import fits
import astropy.units as u
from astropy.wcs import WCS
from dataclasses import dataclass
import numpy as np
import pathlib

from spark.io.models import SparkInstrument, SparkCube, SparkSpec, SparkSlit, ObservationTarget
import spark.constants as sc

JWST_WAVE_UNIT = u.um

@dataclass
class JWSTCube(SparkCube):
    inst_class = SparkInstrument

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

        inst_name = hdu[0].header['INSTRUME']
        # TODO: fix
        if inst_name == 'NIRSPEC':
            filt = hdu[0].header['FILTER']
            grating = hdu[0].header['GRATING']
            instrument = NIRSpecInstrument(name=inst_name, filter=filt, grating=grating)
        else:
            instrument = SparkInstrument(name=inst_name)

        inst_data['target'] = ObservationTarget(
            name = kwargs.pop('name', hdu[0].header['TARGNAME']),
            ra = hdu[0].header['TARG_RA'],
            dec = hdu[0].header['TARG_DEC'],
            z = kwargs.pop('z', 0.0),
            instrument = instrument,
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
class NIRSpecInstrument(SparkInstrument):
    filter: str = None
    grating: str = None


@dataclass
class MIRI_MRS_Cube(JWSTCube):
    pass


@dataclass
class NIRSPEC_IFU_Cube(JWSTCube):
    inst_class = NIRSpecInstrument


@dataclass
class NIRSPEC_MSA_Slit(SparkSlit):
    inst_class = NIRSpecInstrument
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

        header = hdu['PRIMARY'].header
        sci_header = hdu['SCI'].header

        obs_wave_2d = hdu['WAVELENGTH'].data
        obs_wave = obs_wave_2d[0,:]
        is_equal = np.array_equal(obs_wave_2d, np.broadcast_to(obs_wave, obs_wave_2d.shape))
        diff = np.abs(obs_wave_2d - obs_wave[None, :])

        # check if the first wavelength row can be used for the whole slit
        # if not, obs_wave needs to be 2d
        if not is_equal and np.nanmax(diff) > 1e-5:
            obs_wave = obs_wave_2d

        cunit = u.Unit(header.get('CUNIT', 'um'))
        obs_wave = obs_wave * cunit

        bunit = u.Unit(sci_header['BUNIT'])
        slit_flux = hdu['SCI'].data * bunit
        slit_err = hdu['ERR'].data * bunit
        if ('/sr' in str(bunit)) or ('/ sr' in str(bunit)):
            pxar = sci_header['PIXAR_SR'] * u.sr
            slit_flux *= pxar
            slit_err *= pxar

        slit_flux = slit_flux.to(sc.FLUX_UNIT_UM, equivalencies=u.spectral_density(obs_wave))
        slit_err = slit_err.to(sc.FLUX_UNIT_UM, equivalencies=u.spectral_density(obs_wave))
        slit_err[np.isnan(slit_err)] = np.nanmedian(slit_err)


        inst_data['wcs'] = WCS(hdu['SCI'])
        inst_data['wave'] = obs_wave
        inst_data['wave_unit'] = cunit
        inst_data['flux'] = slit_flux
        inst_data['err'] = slit_err
        inst_data['flux_unit'] = sc.FLUX_UNIT_UM

        instrument = NIRSpecInstrument(name=header['INSTRUME'], filter=header['FILTER'], grating=header['GRATING'])
        inst_data['target'] = ObservationTarget(
            name = kwargs.pop('name', sci_header['SRCNAME']),
            ra = sci_header['SRCRA'],
            dec = sci_header['SRCDEC'],
            z = kwargs.pop('z', 0.0),
            instrument = instrument,
        )

        inst_data['pix_scale'] = np.sqrt(sci_header['PIXAR_A2']) * u.arcsec

        inst_data.update(kwargs)
        if hdu_close: hdu.close()
        return cls(**inst_data)


@dataclass
class DJASpec(SparkSpec):
    # Spectrum from the DAWN JWST Archive (https://dawn-cph.github.io/dja/)

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

        header = hdu['SPEC1D'].header
        data = hdu['SPEC1D'].data
        sci_header = hdu['SCI'].header

        wave_unit = u.Unit(header['TUNIT1'])
        flux_unit = u.Unit(header['TUNIT2'])
        wcs = WCS(hdu['SCI'])

        obs_wave = data['wave']*wave_unit
        flux = data['flux']*flux_unit
        err = data['err']*flux_unit

        flux = flux.to(sc.FLUX_UNIT_UM, equivalencies=u.spectral_density(obs_wave))
        err = err.to(sc.FLUX_UNIT_UM, equivalencies=u.spectral_density(obs_wave))
        err[np.isnan(err)] = np.nanmedian(err)

        inst_data['wave'] = obs_wave
        inst_data['wave_unit'] = wave_unit
        inst_data['flux'] = flux
        inst_data['err'] = err
        inst_data['flux_unit'] = sc.FLUX_UNIT_UM

        if name is None:
            name = hdu['SCI'].header['SRCNAME']

        ra, dec = 0.0, 0.0
        for label in ['TARG_RA', 'SRCRA']:
            if label in sci_header:
                ra, dec = sci_header[label], sci_header[label.replace('RA','DEC')]
                break

        if 'INSTRUME' in sci_header:
            instrument = sci_header['INSTRUME']
        else:
            instrument = 'NS_%s_%s'%(sci_header['GRATING'], sci_header['FILTER'])

        inst_data['target'] = ObservationTarget(
            name = name,
            ra = ra,
            dec = dec,
            z = kwargs.pop('z', 0.0),
            instrument = SparkInstrument(name = instrument),
        )

        inst_data.update(kwargs)
        if hdu_close: hdu.close()
        return cls(**inst_data)

