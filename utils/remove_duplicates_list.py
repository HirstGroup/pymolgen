import argparse
import os
import sys

from functools import partial

print = partial(print, flush=True)

parser = argparse.ArgumentParser(description='Remove Duplicates from List of inchi files')

# required arguments
parser.add_argument('-i', '--input', nargs='+', help='Space-separated list of inchi input files', required=True)
parser.add_argument('-o', '--output', help='Output inchi file', required=True)

# optional arguments
parser.add_argument('--delete', action='store_true', help='Delete input file after succesful run', required=False)
parser.add_argument('--read_count', action='store_true', help='Read count from input file', required=False)

args = parser.parse_args()

d = {}

for input in args.input:

    print(input)

    if input == args.output:
        sys.exit('ERROR: Same input and output')

    with open(input) as infile:

        for line in infile:

            inchi = line.strip().split()[0]

            if args.read_count:
                count = int(line.strip().split()[1])
            else:
                count = 1

            if inchi in d.keys():
                d[inchi] += count
            else:
                d[inchi] = count

    d = dict(sorted(d.items(), key=lambda item: item[1], reverse=True))

    with open(args.output, 'w') as outfile:
        for key, val in d.items():
            outfile.write('%s %s\n' %(key, val))

    if args.delete is True:
        os.remove(input)   


