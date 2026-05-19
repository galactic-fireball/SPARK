import spark.constants as sc

def add_ax_labels(ax, unit):
	if unit == 'um':
		ax.set_xlabel(sc.WAVE_UM_LABEL)
		ax.set_ylabel(sc.FLUX_UM_LABEL)
	elif unit == 'AA':
		ax.set_xlabel(sc.WAVE_AA_LABEL)
		ax.set_ylabel(sc.FLUX_AA_LABEL)
	else:
		print('Unsupported unit: %s'%unit)
