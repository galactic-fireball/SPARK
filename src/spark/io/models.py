from abc import ABC, abstractmethod
import astropy.units as u
from dataclasses import dataclass, field
import numpy as np
import pathlib
from typing import Tuple

from spark.utils import deredden

@dataclass
class SparkSpec(ABC):
	file: pathlib.Path = None
	wave: np.ndarray = None
	wave_unit: u.Unit = None
	wave_is_rest: bool = False
	flux: np.ndarray = None
	err: np.ndarray = None
	flux_unit: u.Unit = None
	z: float = 0.0
	ra: float = 0.0
	dec: float = 0.0


	def __post_init__(self):
		if self.z != 0.0:
			self.wave = deredden(self.wave, z=self.z)
			self.wave_is_rest = True


	@classmethod
	@abstractmethod
	def from_fits(cls, hdu, **kwargs):
		pass


	def __str__(self):
		s = 'Spectral Data from: %s' % str(self.file) if not self.file is None else '--'
		s += '\n\tRA: %.04f ; DEC: %.04f ; z = %.04f' % (self.ra, self.dec, self.z)
		s += '\n\twave array: sizeof %d [%s] (%s)' % (len(self.wave), str(self.wave_unit), 'Rest' if self.wave_is_rest else 'Obs')
		s += '\n\tflux array: sizeof %s [%s]' % (self.flux.shape, str(self.flux_unit))
		return s


@dataclass
class SparkCube(SparkSpec):
	shape: Tuple[int, int] = (0,0)

	def __post_init__(self):
		super().__post_init__()
		if (self.shape == (0,0)) and (not self.flux is None):
			self.shape = self.flux.shape[:-1]
