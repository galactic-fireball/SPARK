from astropy.coordinates import SkyCoord
import astropy.units as u
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np

import spark.constants as sc

def add_ax_labels(ax, unit, yscale=0, **kwargs):
    if unit == 'um':
        ax.set_xlabel(sc.WAVE_UM_LABEL, **kwargs)
        ylabel = sc.FLUX_UM_LABEL
        if yscale != 0:
            ylabel = sc.FLUX_UM_SCALE_LABEL%yscale
        ax.set_ylabel(ylabel, **kwargs)
    elif unit == 'AA':
        ax.set_xlabel(sc.WAVE_AA_LABEL, **kwargs)
        ylabel = sc.FLUX_AA_LABEL
        if yscale != 0:
            ylabel = sc.FLUX_AA_SCALE_LABEL%yscale
        ax.set_ylabel(ylabel, **kwargs)
    else:
        print('Unsupported unit: %s'%unit)


def add_wcs_grid_labels(ax):
    lon_edge_frac = 0.12
    lat_edge_frac = 0.31
    offset = 0.8*u.arcsec

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    xmin, xmax = min(xlim), max(xlim)
    ymin, ymax = min(ylim), max(ylim)
    dx = xmax - xmin
    dy = ymax - ymin

    text_kw = dict(color='white',transform=ax.get_transform('world'),rotation_mode='anchor',ha='center',va='center')

    px = xmin + lon_edge_frac * dx
    py = ymin + 0.5 * dy
    c = ax.wcs.pixel_to_world(px, py)

    c1 = SkyCoord(ra=c.ra, dec=c.dec-offset, frame=c.frame)
    c2 = SkyCoord(ra=c1.ra, dec=c1.dec+offset, frame=c.frame)
    x1, y1 = ax.wcs.world_to_pixel(c1)
    x2, y2 = ax.wcs.world_to_pixel(c2)
    lon_angle = np.degrees(np.arctan2(y2-y1,x2-x1))

    # ax.text(c.ra.deg, c.dec.deg, 'TEST LON', rotation=lon_angle, **text_kw)


    px = xmin + 0.5 * dx
    py = ymin + lat_edge_frac * dy
    c = ax.wcs.pixel_to_world(px, py)

    dra = offset / np.cos(c.dec)
    c1 = SkyCoord(ra=c.ra-dra, dec=c.dec, frame=c.frame)
    c2 = SkyCoord(ra=c.ra+dra, dec=c.dec, frame=c.frame)
    x1, y1 = ax.wcs.world_to_pixel(c1)
    x2, y2 = ax.wcs.world_to_pixel(c2)
    lat_angle = np.degrees(np.arctan2(y2-y1,x2-x1))

    # ax.text(c1.ra.value, c1.dec.value, 'TEST LAT', rotation=lat_angle, **text_kw)

    ax.coords.grid(color='0.8', linestyle=':', linewidth=0.5, alpha=0.7)


def show_2d(image, log=False):
    plt.figure()
    if log:
        plt.imshow(image, origin='lower', norm=LogNorm())
    else:
        plt.imshow(image, origin='lower')
    plt.show()
    plt.close()


def show_median_map(cube_flux, log=False):
    medcube = np.nanmedian(cube_flux, axis=0)
    medcube[np.isnan(medcube)] = 0.0
    plt.figure()
    if log:
        plt.imshow(medcube, origin='lower', norm=LogNorm())
    else:
        plt.imshow(medcube, origin='lower')
    plt.show()
    plt.close()
