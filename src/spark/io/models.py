from __future__ import annotations
import astropy.units as u
from astropy.wcs import WCS
from dataclasses import asdict, dataclass, field, fields
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
import pathlib
from photutils.aperture import CircularAperture
from photutils.detection import DAOStarFinder
from scipy import ndimage
from scipy.signal import fftconvolve
from typing import NamedTuple, Tuple

from spark.utils import deredden


@dataclass(kw_only=True)
class SparkInstrument:
    name: str = None


@dataclass(kw_only=True)
class Target:
    name: str = None
    ra: float = 0.0
    dec: float = 0.0
    z: float = 0.0


@dataclass(kw_only=True)
class ObservationTarget(Target):
    instrument: SparkInstrument = None


class Coord(NamedTuple):
    x: int
    y: int

    def __str__(self):
        return f'({self.x}, {self.y})'


def as_coord(c):
    return c if isinstance(c, Coord) else Coord(*c)


@dataclass(kw_only=True)
class SparkSpec():
    file: pathlib.Path = None
    name: str = None
    wave: np.ndarray = None
    wave_unit: u.Unit = None
    wave_is_rest: bool = False
    flux: np.ndarray = None
    err: np.ndarray = None
    flux_unit: u.Unit = None
    target: Target = None


    def __post_init__(self):
        print('SparkSpec __post_init__')
        if self.target.z != 0.0:
            self.wave = deredden(self.wave, z=self.target.z)
            self.wave_is_rest = True


    @classmethod
    def from_fits(cls, hdu, **kwargs):
        raise Exception('Unimplemented')


    def to_dict(self):
        return asdict(self)


    def __str__(self):
        s = 'Spectral Data from: %s' % str(self.file) if not self.file is None else '--'
        s += '\n\tRA: %.04f ; DEC: %.04f ; z = %.04f' % (self.target.ra, self.target.dec, self.target.z)
        s += '\n\tWave: [%.04f - %.04f] %s (%s)' % (self.wave[0], self.wave[-1], str(self.wave_unit), 'Rest' if self.wave_is_rest else 'Obs')
        s += '\n\tFlux: [%.04f - $.04f] %s' % (self.flux[0], self.flux[-1], str(self.flux_unit))
        return s


    def quick_view(self):
        if len(self.flux.shape) > 1:
            return
        fix, ax = plt.subplots()
        ax.plot(self.wave, self.flux, linewidth=1.0, color='black')
        plt.show()



@dataclass(kw_only=True)
class SparkSpaxel(SparkSpec):
    cube: SparkCube = None
    coord: Coord = None
    name: str = None

    def __post_init__(self):
        self.wave_is_rest = self.cube.wave_is_rest

        if self.name is None:
            self.name = 'spaxel_%d_%d'%(self.coord.x,self.coord.y)

        for key, val in self.__dict__.items():
            if val is None:
                setattr(self, key, getattr(self.cube, key))


    def __str__(self):
        s = 'Spaxel (%d, %d) Data from: %s' % (self.coord.x, self.coord.y, self.file.stem if not self.file is None else '--')
        s += '\n\twave array: sizeof %d' % (len(self.wave))
        s += '\n\tflux array: sizeof %s' % (self.flux.shape)
        return s


@dataclass(kw_only=True)
class SparkAperture(SparkSpec):
    cube: SparkCube = None
    center: Coord | tuple = None
    weights: np.ndarray = None
    method: str = 'sum'
    name: str = None

    def __post_init__(self):
        self.wave_is_rest = self.cube.wave_is_rest
        self.center = as_coord(self.center)
        self.set_weights()

        flux_funcs = {
            'sum': self.sum_flux,
            'mean': self.mean_flux,
        }

        if not self.method in flux_funcs:
            raise Exception('Invalid flux method: %s'%self.method)
        flux_funcs[self.method]()

        for key, val in self.__dict__.items():
            if val is None:
                setattr(self, key, getattr(self.cube, key))


    def set_weights(self):
        ny, nx = self.cube.flux.shape[1:]
        # weights of all spaxels = 0.0
        self.weights = np.zeros((ny,nx), dtype=float)


    def sum_flux(self):
        w = self.weights[None, :, :]
        self.flux = np.nansum(self.cube.flux*w, axis=(1,2))
        self.err = np.sqrt(np.nansum((self.cube.err*w)**2, axis=(1,2)))


    def mean_flux(self):
        w = self.weights[None, :, :]
        norm = np.nansum(self.weights)
        self.flux = np.nansum(self.cube.flux*w, axis=(1,2)) / norm
        self.err = np.sqrt(np.nansum((self.cube.err*w)**2, axis=(1,2))) / norm


    def show_weights(self):
        plt.figure()
        plt.imshow(self.weights, origin='lower', cmap='grey')
        plt.show()


@dataclass(kw_only=True)
class SparkCircularAperture(SparkAperture):
    radius: float

    def __post_init__(self):
        super().__post_init__()
        if self.name is None:
            self.name = 'circap_x%d_y%d_r%d'%(self.center.x,self.center.y,self.radius)


    def set_weights(self):
        ny, nx = self.cube.flux.shape[1:]
        ap = CircularAperture((self.center.x, self.center.y), self.radius)
        self.weights = ap.to_mask(method='exact').to_image((ny, nx))


@dataclass(kw_only=True)
class SparkRectangularAperture(SparkAperture):
    width: float
    height: float = -1.0

    def __post_init__(self):
        if self.height < 0.0: self.height = self.width
        super().__post_init__()

        if self.name is None:
            self.name = 'rectap_x%d_y%d_w%d_h%d'%(self.center.x,self.center.y,self.width,self.height)


    def set_weights(self):
        ny, nx = self.cube.flux.shape[1:]
        x, y = np.meshgrid(np.arange(nx), np.arange(ny))
        self.weights = ((np.abs(x - self.center.x) <= self.width / 2) & (np.abs(y - self.center.y) <= self.height / 2)).astype(float)


@dataclass(kw_only=True)
class SparkCube(SparkSpec):
    spaxel_class = SparkSpaxel
    ap_shapes = {
        'circular': SparkCircularAperture,
        'rectangular': SparkRectangularAperture,
    }

    shape: Tuple[int, int] = (0,0) # (nx,ny)
    wcs: WCS = None

    def __post_init__(self):
        super().__post_init__()
        if (self.shape == (0,0)) and (not self.flux is None):
            self.shape = self.flux.shape[2], self.flux.shape[1] # flux.shape = (nz,ny,nx)

        # self.flux[np.isnan(self.flux)] = 0.0
        # self.err[np.isnan(self.err)] = np.nanmedian(self.err)


    def get_median_map(self):
        medcube = np.nanmedian(self.flux, axis=0)
        medcube[np.isnan(medcube)] = 0.0
        return medcube


    def show_median_map(self):
        medcube = self.get_median_map()
        plt.figure()
        plt.imshow(medcube.value, origin='lower', norm=LogNorm())
        plt.show()


    def quick_view(self):
        self.show_median_map()


    def spax(self, x, y=None):
        if not y is None:
            coord = Coord(x=x,y=y)
        elif isinstance(x, tuple) and len(x) > 1:
            coord = Coord(x=x[0],y=x[1])
        elif isinstance(x, Coord):
            coord = x
        else:
            return None

        spax_data = {
            'cube': self,
            'coord': coord,
            'wave': self.wave,
            'flux': self.flux[:,coord.y,coord.x], # NumPy indexes (row, column) -> (y,x)
            'err': self.err[:,coord.y,coord.x]
        }

        return self.spaxel_class(**spax_data)


    def aperture(self, shape, center, **kwargs):
        if not shape in self.ap_shapes:
            print('Invalid aperture shape: %s'%shape)
            return None
        shape_class = self.ap_shapes[shape]
        vfields = {f.name for f in fields(shape_class)}

        ap_data = {
            'cube': self,
            'center': as_coord(center),
            'wave': self.wave,
        }
        ap_data.update({k:v for k,v in kwargs.items() if k in vfields})
        return self.ap_shapes[shape](**ap_data)


    def get_brightest_spaxel(self, mode='median'):
        dao = DAOStarFinder(fwhm=3.0, threshold=1e-17*self.flux_unit)
        sources = dao(self.get_median_map())
        source = sources[np.argmax(sources['flux'])]
        return self.spax(round(source['xcentroid']), round(source['ycentroid']))


    def circular_smooth(self, radius):
        cube = self.flux.copy().value

        ny, nx = cube.shape[1:]

        y, x = np.ogrid[-radius:radius+1, -radius:radius+1]
        mask = x**2 + y**2 <= radius**2
        kernel = mask.astype(float)
        kernel /= kernel.sum()

        smoothed_cube = np.empty_like(cube)
        for i in range(cube.shape[0]):
            smoothed_cube[i] = ndimage.convolve(cube[i], kernel, mode='reflect')
            # smoothed_cube[i] = fftconvolve(cube[i], kernel, mode='same')
        self.flux = smoothed_cube*self.flux.unit


@dataclass(kw_only=True)
class SparkSlit(SparkSpec):
    nx: int = 0

