import sys,os
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Remove Duplicates from List of inchi files')
parser.add_argument('-i', '--input', nargs='+', help='Space-separated list of inchi input files',required=True)
parser.add_argument('-o', '--output', help='Output file',required=True)

args = parser.parse_args()

d = set()

for input in args.input:

    print(input)

    infile = open(input)

    a = []

    for line in infile:
        a = line.strip('\n')
        if a not in d:
            d.add(a)


with open(args.output, 'w') as outfile:
    for i in d:
        outfile.write('%s\n' %i)