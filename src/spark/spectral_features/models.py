from dataclasses import asdict, dataclass, field

@dataclass
class SpectralLine:
	name: str
	center: float