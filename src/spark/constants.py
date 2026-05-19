import astropy.units as u

WAVE_UNIT_AA = u.AA
WAVE_UNIT_UM = u.um

INT_FLUX = u.erg / u.s / (u.cm**2)
FLUX_UNIT_UM = u.erg / u.s / (u.cm**2) / u.um
FLUX_UNIT_AA = u.erg / u.s / (u.cm**2) / u.AA


FLUX_AA_LABEL = r'$f_\lambda (~\rm{erg}~\rm{cm}^{-2}~\rm{s}^{-1}\rm{\AA}^{-1})$'
WAVE_AA_LABEL = r'$\lambda_{rest} (\AA)$'
FLUX_UM_LABEL = r'$f_\lambda (~\rm{erg}~\rm{cm}^{-2}~\rm{s}^{-1}\rm{\mu m}^{-1})$'
WAVE_UM_LABEL = r'$\lambda_{rest} (\mu m)$'
INT_FLUX_LABEL = r'$F (~\rm{erg}~\rm{cm}^{-2}~\rm{s}^{-1})$'
