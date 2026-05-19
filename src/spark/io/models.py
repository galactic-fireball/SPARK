from abc import ABC, abstractmethod
import astropy.units as u
from astropy.wcs import WCS
from dataclasses import asdict, dataclass, field
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
import pathlib
from photutils.aperture import CircularAperture, ApertureStats
from photutils.detection import DAOStarFinder
from scipy import ndimage
from scipy.signal import fftconvolve
from typing import NamedTuple, Tuple

from spark.utils import deredden


@dataclass
class SparkInstrument:
    name: str = None


@dataclass
class Target:
    name: str = None
    ra: float = 0.0
    dec: float = 0.0
    z: float = 0.0


@dataclass
class ObservationTarget(Target):
    instrument: SparkInstrument = None


class Coord(NamedTuple):
    x: int
    y: int

    def __str__(self):
        return f'({self.x}, {self.y})'


@dataclass
class SparkSpec(ABC):
    file: pathlib.Path = None
    wave: np.ndarray = None
    wave_unit: u.Unit = None
    wave_is_rest: bool = False
    flux: np.ndarray = None
    err: np.ndarray = None
    flux_unit: u.Unit = None
    target: Target = None


    def __post_init__(self):
        if self.target.z != 0.0:
            self.wave = deredden(self.wave, z=self.target.z)
            self.wave_is_rest = True


    @classmethod
    @abstractmethod
    def from_fits(cls, hdu, **kwargs):
        pass


    def to_dict(self):
        return asdict(self)


    def __str__(self):
        s = 'Spectral Data from: %s' % str(self.file) if not self.file is None else '--'
        s += '\n\tRA: %.04f ; DEC: %.04f ; z = %.04f' % (self.target.ra, self.target.dec, self.target.z)
        s += '\n\tWave: [%.04f - %.04f] %s (%s)' % (self.wave[0], self.wave[-1], str(self.wave_unit), 'Rest' if self.wave_is_rest else 'Obs')
        s += '\n\tFlux: [%.04f - $.04f] %s' % (self.flux[0], self.flux[-1], str(self.flux_unit))
        return s


@dataclass
class SparkCube(SparkSpec):
    # flux.shape = (nx, ny, nwave)
    shape: Tuple[int, int] = (0,0)
    wcs: WCS = None

    def __post_init__(self):
        super().__post_init__()
        if (self.shape == (0,0)) and (not self.flux is None):
            self.shape = self.flux.shape[0], self.flux.shape[1]


    def get_median_map(self):
        medcube = np.nanmedian(self.flux, axis=2)
        medcube[np.isnan(medcube)] = 0.0
        return medcube


    def show_median_map(self):
        medcube = self.get_median_map()
        plt.figure()
        plt.imshow(medcube.value, origin='lower', norm=LogNorm())
        plt.show()


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
            'flux': self.flux[coord.x,coord.y,:],
            'err': self.err[coord.x,coord.y,:]
        }

        return SparkSpaxel(**spax_data)


    def get_brightest_spaxel(self, mode='median'):
        dao = DAOStarFinder(fwhm=3.0, threshold=1e-17*self.flux_unit)
        sources = dao(self.get_median_map())
        source = sources[np.argmax(sources['flux'])]
        return self.spax(round(source['ycentroid']), round(source['xcentroid']))


    def get_circular_aperture(self, ap_coord, radius, sum_method='exact'):
        aperture = CircularAperture((ap_coord.x,ap_coord.y), r=radius)
        apspec = np.zeros(len(self.wave))*self.flux_unit
        aperr = np.zeros(len(self.wave))*self.flux_unit
        for i in range(len(self.wave)):
            apstat = ApertureStats(self.flux[:,:,i], aperture, error=self.err[:,:,i], sum_method=sum_method, subpixels=100)
            apspec[i] = apstat.sum
            aperr[i] = apstat.sum_err
        return apspec, aperr


    def circular_smooth(self, radius):
        cube = self.flux.copy().T.value

        nx, ny = cube.shape[1:]

        y, x = np.ogrid[-radius:radius+1, -radius:radius+1]
        mask = x**2 + y**2 <= radius**2
        kernel = mask.astype(float)
        kernel /= kernel.sum()

        # ky, kx = kernel.shape
        # pad_y = ky // 2
        # pad_x = kx // 2
        # padded_cube = np.pad(cube, ((0,0), (pad_y,pad_y), (pad_x,pad_x)), mode='reflect')

        smoothed_cube = np.empty_like(cube)
        for i in range(cube.shape[0]):
            smoothed_cube[i] = ndimage.convolve(cube[i], kernel, mode='reflect')
            # smoothed_cube[i] = fftconvolve(cube[i], kernel, mode='same')
        # return smoothed_cube.T
        self.flux = smoothed_cube.T*self.flux.unit


@dataclass
class SparkSlit(SparkSpec):
    nx: int = 0


@dataclass
class SparkSpaxel(SparkSpec):
    cube: SparkCube = None
    coord: Coord = None


    def __post_init__(self):
        for attr in ['file', 'wave_unit', 'wave_is_rest', 'flux_unit', 'z', 'ra', 'dec']:
            if (getattr(self, attr, None) is None) or (getattr(self, attr, None) == 0.0):
                setattr(self, attr, getattr(self.cube, attr, None))


    def __str__(self):
        s = 'Spaxel (%d, %d) Data from: %s' % (self.coord.x, self.coord.y, self.file.stem if not self.file is None else '--')
        s += '\n\twave array: sizeof %d' % (len(self.wave))
        s += '\n\tflux array: sizeof %s' % (self.flux.shape)
        return s


    def from_fits(cls, hdu, **kwargs):
        raise Exception('Unimplemented')


