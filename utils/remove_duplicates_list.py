import sys,os
import pandas as pd
import numpy as np
import argsparse

parser = argparse.ArgumentParser(description='Remove Duplicates from List of inchi files')
parser.add_argument('-i', nargs='+', help='Space-separated list of inchi input files',required=True)
parser.add_argument('-o', help='Output file',required=True)

d = set()

for input in args.input:

    infile = open(input)

    a = []

    for line in infile:
        a = line.strip('\n')
        if a not in d:
            d.add(a)


with open(args.output, 'w') as outfile:
    for i in d:
        outfile.write('%s\n' %i)