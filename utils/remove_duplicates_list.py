import sys,os
import pandas as pd
import numpy as np
import argparse

from functools import partial

print = partial(print, flush=True)

parser = argparse.ArgumentParser(description='Remove Duplicates from List of inchi files')
parser.add_argument('-i', '--input', nargs='+', help='Space-separated list of inchi input files',required=True)
parser.add_argument('-o', '--output', help='Output file',required=True)

args = parser.parse_args()

d = {}

for input in args.input:

    print(input)

    infile = open(input)

    for line in infile:

        if len(line.split()) == 1:

            a = line.strip('\n')
            d[a] = 1

        elif len(line.split()) == 2:

            inchi = line.split()[0]
            probability = float(line.strip().split()[1])

            if inchi in d.keys():
                d[inchi] += probability
            else:
                d[inchi] = probability

        else:
            sys.error('Len line split > 2')

d = dict(sorted(d.items(), key=lambda item: item[1], reverse=True))

with open(args.output, 'w') as outfile:
    for key, val in d.items():
        outfile.write('%s %s\n' %(key, val))