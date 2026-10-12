from astroquery.sdss import SDSS
import numpy as np
import re

LATEST_DR = 20

# TODO: coord option
# def download_sdss_spec(spec_name=None, plate=None, mjd=None, fiberid=None, dest=None):
def download_sdss_spec(spec_name, dest=None):
    if not spec_name is None:
        _, plate, mjd, fiberid = spec_name.split('-')

    sp = SDSS.get_spectra(plate=int(plate), mjd=int(mjd), fiberID=int(fiberid), data_release=LATEST_DR)
    if sp is None or len(sp) == 0:
        print('Could not find spec')
        return None

    if not dest is None:
        if len(sp) > 1:
            if not dest.is_dir():
                print('Multiple results found, using dest parent')
                dest = dest.parent
            for s in sp:
                s.writeto(dest.joinpath(spec_name+'.fits'), overwrite=True)
        else:
            if dest.is_dir():
                dest = dest.joinpath(spec_name+'.fits')
            sp[0].writeto(dest, overwrite=True)

    return sp


# See: https://github.com/astropy/astroquery/blob/main/astroquery/sdss/core.py
def get_spec_url(spec_name):
    if not spec_name is None:
        _, plate, mjd, fiberid = spec_name.split('-')

    base_url = 'https://data.sdss.org/sas'
    SPECTRA_URL_SUFFIX = ('{base}/dr{dr}/spectro/sdss/redux/{run2d}/{spectra_path}/{plate:0>4d}/spec-{plate:0>4d}-{mjd}-{fiber:04d}.fits')
    parse_run2d = re.compile(r'v(?P<major>[0-9]+)_(?P<minor>[0-9]+)_(?P<bugfix>[0-9]+)')


    linkstr = SPECTRA_URL_SUFFIX

    sp = SDSS.query_specobj(plate=int(plate), mjd=int(mjd), fiberID=int(fiberid), data_release=LATEST_DR)
    if sp is None or len(sp) == 0:
        print('Could not find spec')
        return []

    res = []
    for s in sp:
        if isinstance(s['run2d'], bytes):
            run2d = s['run2d'].decode()
        elif isinstance(s['run2d'], (np.integer, int)):
            run2d = str(s['run2d'])
        else:
            run2d = s['run2d']

        format_args = {
            'base': base_url,
            'dr': LATEST_DR,
            'run2d': run2d,
            'spectra_path': 'spectra/full' if run2d in ('26', '103', '104') else 'spectra',
            'mjd': s['mjd'],
        }

        try:
            format_args['plate'] = s['plate']
            format_args['fiber'] = s['fiberID']
        except:
            format_args['fieldid'] = s['fieldID']
            format_args['catalogid'] = s['catalogID']

        match_run2d = parse_run2d.match(run2d)
        if not match_run2d is None:
            major = int(match_run2d.group('major'))
            if major > 5:
                linkstr = linkstr.replace('/{plate:0>4d}/', '/{fieldid:0>4d}p/{mjd:5d}/')
                linkstr = linkstr.replace('spec-{plate:0>4d}-{mjd}-{fiber:04d}.fits',
                                              'spec-{fieldid:0>4d}-{mjd:5d}-{catalogid:0>11d}.fits')

        link = linkstr.format(**format_args)
        res.append(link)

    return res

