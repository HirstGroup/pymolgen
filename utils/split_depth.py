import argparse
import os

from functools import partial

print = partial(print, flush=True)

parser = argparse.ArgumentParser(description='Split inchi file into separate files according to depth')

# required arguments
parser.add_argument('-i','--input', help='Input inchi file', required=True)
parser.add_argument('-o','--output', help='Output inchi file base name', required=True)

# optional arguments
parser.add_argument('--delete', action='store_true', help='Delete input file after succesful run', required=False)

args = parser.parse_args()

print(args.input)

infile = open(args.input)

depth_list = [f'd{x}' for x in list(range(3,13))]

depth_count = {x:0 for x in depth_list}

outfile_list = []


for depth in depth_list:

	outfile = open(f'{args.output}_{depth}.inchi', 'w')

	outfile_list.append(outfile)


for n, line in enumerate(infile):

	inchi, count, line_depth = line.strip().split()

	for x, depth in enumerate(depth_list):

		if depth == line_depth:

			outfile_list[x].write(line)

			depth_count[depth] += 1

total_count = 0

for key, val in depth_count.items():

	total_count += val

assert n+1 == total_count

if args.delete and n+1 == total_count:

	os.remove(args.input)





