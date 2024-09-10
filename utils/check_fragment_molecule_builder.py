import argparse
import glob

parser = argparse.ArgumentParser(description='Check output from fragment_molecule_builder when run in parallel with split inputs')

parser.add_argument('-i','--input', nargs='+', help='Input number for split files', required=True)

args = parser.parse_args()

all_depths = list(range(2,13))

print('split depths')

for i in args.input:

	depth_list = []

	files = glob.glob(f'phenylisoxazole-20-21-systematic-depth1_{i}/phenylisoxazole-20-21-systematic-bp10-depth1_{i}-depth*.inchi')

	for file in files:

		name = file[:-6]

		depth = name.rfind('depth')
		
		depth_list.append(int(name[depth+len('depth'):]))

	print(i, end=' ')

	for depth in all_depths:
		if depth not in depth_list:
			print('0', end=' ')
		else:
			print('1', end=' ')

	print()

print('Normal termination')
