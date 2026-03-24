import pathlib

from spark.io import load_cube

TEST_DATA_DIR = pathlib.Path(__file__).resolve().parent.joinpath('test_data')

def nirspec_basic_test():
	j1201_z = 0.003542211
	test_file = TEST_DATA_DIR.joinpath('j1201_nirspec_170.fits')
	cube = load_cube(test_file, fmt='NIRSPEC_IFU', z=j1201_z)
	print(cube)


def main():
	nirspec_basic_test()
	# miri_basic_test()


if __name__ == '__main__':
	main()
