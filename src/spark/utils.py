from astropy import constants as const
from astropy.coordinates import SkyCoord
from astropy.cosmology import LambdaCDM
from astroquery.irsa_dust import IrsaDust
import astropy.units as u
import numpy as np

def deredden(wave, z=0.0):
    return wave / (1 + z)


def redden(wave, z=0.0):
    return wave * (1 + z)


def to_rest(wave, z=0.0):
    return deredden(wave, z=z)


def to_obs(wave, z=0.0):
    return redden(wave, z=z)


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


def to_velocity(wave, lam0):
    return const.c.to(u.km/u.s) * (wave - lam0) / lam0


def calc_redshift(line_center, feature_center):
    return feature_center/line_center - 1


def unlog_flux(flux, err):
    unlog_flux = 10**flux
    unlog_err = np.log(10) * unlog_flux * err

    return unlog_flux, unlog_err


def get_target_ebv(ra, dec):
    co = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='fk5')
    try:
        table = IrsaDust.get_query_table(co, section='ebv')
        ebv = table['ext SandF mean'][0]
    except:
        ebv = 0.04

    if (ebv >= 1.0):
        ebv = 0.04

    return ebv


def ccm_unred(wave, flux, ebv, r_v=""):

    wave = np.array(wave, float)
    flux = np.array(flux, float)

    # if wave.size != flux.size: raise TypeError('ERROR - wave and flux vectors must be the same size')

    if not bool(r_v): r_v = 3.1

    x = 10000.0 / wave
    # Correction invalid for x>11:
    if np.any(x > 11):
        return flux

    npts = wave.size
    a = np.zeros(npts, float)
    b = np.zeros(npts, float)

    ###############################
    # Infrared

    good = np.where((x > 0.3) & (x < 1.1))
    a[good] = 0.574 * x[good] ** (1.61)
    b[good] = -0.527 * x[good] ** (1.61)

    ###############################
    # Optical & Near IR

    good = np.where((x >= 1.1) & (x < 3.3))
    y = x[good] - 1.82

    c1 = np.array([1.0, 0.104, -0.609, 0.701, 1.137, \
                   -1.718, -0.827, 1.647, -0.505])
    c2 = np.array([0.0, 1.952, 2.908, -3.989, -7.985, \
                   11.102, 5.491, -10.805, 3.347])

    a[good] = np.polyval(c1[::-1], y)
    b[good] = np.polyval(c2[::-1], y)

    ###############################
    # Mid-UV

    good = np.where((x >= 3.3) & (x < 8))
    y = x[good]
    F_a = np.zeros(np.size(good), float)
    F_b = np.zeros(np.size(good), float)
    good1 = np.where(y > 5.9)

    if np.size(good1) > 0:
        y1 = y[good1] - 5.9
        F_a[good1] = -0.04473 * y1 ** 2 - 0.009779 * y1 ** 3
        F_b[good1] = 0.2130 * y1 ** 2 + 0.1207 * y1 ** 3

    a[good] = 1.752 - 0.316 * y - (0.104 / ((y - 4.67) ** 2 + 0.341)) + F_a
    b[good] = -3.090 + 1.825 * y + (1.206 / ((y - 4.62) ** 2 + 0.263)) + F_b

    ###############################
    # Far-UV

    good = np.where((x >= 8) & (x <= 11))
    y = x[good] - 8.0
    c1 = [-1.073, -0.628, 0.137, -0.070]
    c2 = [13.670, 4.257, -0.420, 0.374]
    a[good] = np.polyval(c1[::-1], y)
    b[good] = np.polyval(c2[::-1], y)

    # Applying Extinction Correction

    a_v = r_v * ebv
    a_lambda = a_v * (a + b / r_v)

    funred = flux * 10.0 ** (0.4 * a_lambda)

    return funred