import json
import numpy as np
import pathlib
import shutil


def generate_cube_viewer(cube, out_dir):
	bin_dir = out_dir.joinpath('bin')
	bin_dir.mkdir(parents=True, exist_ok=True)

	wave_data_path = bin_dir.joinpath('wavelength.bin')
	wave_data = np.asarray(cube.wave.value, dtype=np.float32)
	wave_data = np.ascontiguousarray(wave_data)
	wave_data.tofile(wave_data_path)

	flux_data_path = bin_dir.joinpath('flux.bin')
	flux_data = np.asarray(cube.flux.value, dtype=np.float32)
	flux_data = np.ascontiguousarray(flux_data)
	flux_data.tofile(flux_data_path)

	metadata = {'shape':cube.flux.shape,}
	with open(out_dir.joinpath('metadata.json'), 'w') as f:
		json.dump(metadata, f)

	cwd = pathlib.Path(__file__).resolve().parent
	for path in ['cube_viewer.html', 'static']:
		src = cwd.joinpath(path)
		dest = out_dir.joinpath(path)
		dest.parent.mkdir(parents=True, exist_ok=True)
		if src.is_dir():
			shutil.copytree(src, dest, dirs_exist_ok=True)
		else:
			shutil.copy(src, dest)



def generate_cube_viewer_test(cube, out_dir):
	bin_dir = out_dir.joinpath('bin')
	bin_dir.mkdir(parents=True, exist_ok=True)

	rng = np.random.default_rng()
	test_data = rng.random((20,5,5)) # (nz,ny,nx)
	test_data = np.asarray(test_data, dtype=np.float32)
	test_data = np.ascontiguousarray(test_data)

	flux_data_path = bin_dir.joinpath('flux.bin')
	test_data.tofile(flux_data_path)

	test_data = rng.random((20,))
	test_data = np.asarray(test_data, dtype=np.float32)
	test_data = np.ascontiguousarray(test_data)

	wave_data_path = bin_dir.joinpath('wavelength.bin')
	test_data.tofile(wave_data_path)

	metadata = {'shape':(20,5,5),}
	with open(out_dir.joinpath('metadata.json'), 'w') as f:
		json.dump(metadata, f)

	cwd = pathlib.Path(__file__).resolve().parent
	for path in ['cube_viewer.html', 'static']:
		src = cwd.joinpath(path)
		dest = out_dir.joinpath(path)
		dest.parent.mkdir(parents=True, exist_ok=True)
		if src.is_dir():
			shutil.copytree(src, dest, dirs_exist_ok=True)
		else:
			shutil.copy(src, dest)
