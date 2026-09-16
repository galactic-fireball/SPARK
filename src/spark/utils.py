from astropy import constants as const
from astropy.cosmology import LambdaCDM
import astropy.units as u
import numpy as np

def deredden(wave, z=0.0):
    return wave / (1 + z)


def redden(wave, z=0.0):
    return wave * (1 + z)


def flux_to_lum(flux, z, cosmo=None):
    if cosmo is None:
        cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
    dL = cosmo.luminosity_distance(z).to(u.cm)
    if not isinstance(flux, u.Quantity):
        dL = dL.value
    return 4 * np.pi * dL ** 2 * (flux)


def get_voff(lam, lam2):
    dlam = lam2 - lam
    return (dlam / lam) * const.c.to(u.km/u.s)


def calc_redshift(line_center, feature_center):
    return feature_center/line_center - 1


def unlog_flux(flux, err):
    unlog_flux = 10**flux
    unlog_err = np.log(10) * unlog_flux * err

    return unlog_flux, unlog_err
