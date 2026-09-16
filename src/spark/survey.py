from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import pandas as pd


def find_duplicates(catalog, ra='ra', dec='dec', name_col=None, cunit=u.deg, search_rad=1*u.deg, ounit=u.deg):

	ar_to_df = None
	if isinstance(catalog, pd.DataFrame):
		search_catalog = SkyCoord(ra=np.array(catalog[ra])*cunit, dec=np.array(catalog[dec])*cunit)
		ar_to_df = dict(zip(list(range(len(catalog))), catalog.index))

	all_matches = {k:([],[]) for k in range(len(search_catalog))}

	nn = 2
	while True:
		matches = search_catalog.match_to_catalog_sky(search_catalog, nthneighbor=nn)
		matched = np.where(matches.angular_separation <= search_rad)[0]
		if len(matched) == 0:
			break

		idxs = matches.indices_to_catalog
		dists = matches.angular_separation

		for i, m in enumerate(matched):
			all_matches[m][0].append(idxs[m])
			all_matches[m][1].append(dists[m].to(ounit))

		nn += 1

	if not ar_to_df is None:
		df_matches = {}
		for i, m in enumerate(all_matches.values()):
			k = ar_to_df[i]
			idxs = [ar_to_df[idx] for idx in m[0]]
			if not name_col is None:
				k = catalog.loc[k][name_col]
				idxs = [catalog.loc[idx][name_col] for idx in idxs]

			df_matches[k] = (
				idxs,
				m[1].copy()
			)

		all_matches = df_matches

	return all_matches
