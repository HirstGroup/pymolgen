import argparse
import os

def count_depth(line):
	depth = len(line.strip().split(':')[0].split('-')) - 1
	return depth


def count_fragments_txt(input):

	depth_count = {}

	for i in range(3,13):
		depth_count[i] = 0

	with open(input) as infile:

		for line in infile:
			depth = count_depth(line)

			depth_count[depth] += 1

	return depth_count


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Count number of fragments in inchi file by taking information from txt file')

	# required arguments
	parser.add_argument('-i','--input', help='Input txt file', required=True)

	args = parser.parse_args()

	depth_count = count_fragments_txt(args.input)

	print(args.input, depth_count)







