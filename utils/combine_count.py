import argparse
import ast

def combine_count(input):

	total_count = {}

	for i in range(3,13):
		total_count[i] = 0

	print(total_count)

	with open(input) as infile:

		for line in infile:

			part = ' '.join(line.strip().split()[1:])

			d = ast.literal_eval(part)

			for depth in d:

				total_count[depth] += d[depth]

			print(d)

		return total_count


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Combine count of number of molecules per depth')

	# required arguments
	parser.add_argument('-i','--input', help='Input txt file', required=True)

	args = parser.parse_args()

	total_count = combine_count(args.input)

	print(args.input, total_count)

	for key, val in total_count.items():
		print(key, val)