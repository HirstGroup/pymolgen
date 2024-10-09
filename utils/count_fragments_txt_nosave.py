import argparse
import os

def count_depth(line):
	depth = len(line.strip().split(':')[0].split('-')) - 1
	return depth


def get_mw(line):
	mw = float(line.strip().split(':')[-1])
	return mw


def count_fragments_txt(input, mw_threshold=None):

	depth_count = {}

	for i in range(3,13):
		depth_count[i] = 0

	with open(input) as infile:

		for line in infile:
			depth = count_depth(line)

			if mw_threshold is not None:

				mw = get_mw(line)

				if mw > mw_threshold:

					continue

			depth_count[depth] += 1

	return depth_count


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Count number of fragments in txt file')

	# required arguments
	parser.add_argument('-i','--input', help='Input txt file', required=True)

	# optional arguments
	parser.add_argument('--mw_threshold', type=float, help='Molecular weight threshold', required=False)

	args = parser.parse_args()

	depth_count = count_fragments_txt(args.input, mw=args.mw_threshold)

	print(args.input, depth_count)







