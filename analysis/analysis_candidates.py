import sys,os
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Analyse candidates inchi file')
parser.add_argument('-i', '--input', nargs='+', help='Space-separated list of inchi input files',required=True)
parser.add_argument('-o', '--output', help='Output file',required=True)

args = parser.parse_args()

d = {}

for input in args.input:

    print(input)

    infile = open(input)

    a = []

    for line in infile:
        a = line.strip('\n')
        if a not in d:
            d[a] = 1
        else:
            d[a] += 1


d = {k: v for k, v in sorted(d.items(), key=lambda item: item[1], reverse=True)}

with open(args.output, 'w') as outfile:
    for key, val in d.items():
        outfile.write('%s %s\n' %(key, val))