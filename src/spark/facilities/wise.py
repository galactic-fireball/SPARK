from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.irsa import Irsa
import numpy as np

ALLWISE_CATALOG = 'allwise_p3as_psd'
DEFAULT_SEARCH_RAD_ARCSEC = 3*u.arcsec
WISE_COLS = ['w1mpro', 'w1sigmpro', 'w2mpro', 'w2sigmpro',
             'w3mpro', 'w3sigmpro', 'w4mpro', 'w4sigmpro',]

def search_wise(coord, catalog=ALLWISE_CATALOG, search_radius=DEFAULT_SEARCH_RAD_ARCSEC):
    results = []

    if not isinstance(coord, SkyCoord):
        # TODO
        print('non-SkyCoord unimplemented')
        return results

    query_results = Irsa.query_region(coord, catalog=catalog, radius=search_radius)
    for res in query_results:
        wise_dict = {}

        for band in range(1,5):
            mag_attr = 'w%dmpro' % band
            mag_s_attr = 'w%dsigmpro' % band
            wise_dict[mag_attr.upper()] = res[mag_attr]
            wise_dict[mag_s_attr.upper()] = res[mag_s_attr]

            flux_attr = 'W%d_FLUX' % band
            flux_s_attr = 'W%d_FLUX_sigma' % band
            flux, flux_sigma = mag_to_flux(res[mag_attr], res[mag_s_attr], 'W%d'%band)
            wise_dict[flux_attr] = flux
            wise_dict[flux_s_attr] = flux_sigma

        for lband in range(1,4):
            hband = lband+1
            wise_dict['W%d%d'%(lband,hband)] = res['w%dmpro'%lband] - res['w%dmpro'%hband]

        results.append(wise_dict)

    return results


def mag_to_flux(wise, wise_err, band):
    Fnu0 = {'W1': 309.540*u.Jy, 'W2': 171.787*u.Jy, 'W3': 31.674*u.Jy, 'W4': 8.363*u.Jy}
    Fnu = Fnu0[band] * 10**(-wise/2.5)
    Fnu_sigma = np.abs(Fnu0[band] * -np.log(10)/2.5 * 10 ** (-wise/2.5)) * wise_err

    flux_unit = u.erg/u.s/(u.cm**2)/u.Hz
    Fnu = Fnu.to(flux_unit)
    Fnu_sigma = Fnu_sigma.to(flux_unit)

    # "Integrate" passband
    F = Fnu.value * 1.4653e13
    F_sigma = np.sqrt((1.4653e13 * Fnu_sigma.value)**2 + (Fnu.value * 1.1759e10)**2)

    return F, F_sigma