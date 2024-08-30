import argparse
import sys,os

from functools import partial

print = partial(print, flush=True)

parser = argparse.ArgumentParser(description='Check difference between two inchi files')
parser.add_argument('-i1', '--input1', help='First inchi file',required=True)
parser.add_argument('-i2', '--input2', help='Second inchi file',required=True)
parser.add_argument('-o1', '--output1', help='First output file with inchis from input1 not in input2',required=True)
parser.add_argument('-o2', '--output2', help='First output file with inchis from input2 not in input1',required=True)

args = parser.parse_args()

if len(set([args.input1, args.output1, args.output2])) != 3:
    raise Exception('Same input and/or output files')

if len(set([args.input2, args.output1, args.output2])) != 3:
    raise Exception('Same input and/or output files')

infile1 = open(args.input1)

lines1 = infile1.readlines()

set1 = set([i.strip().split()[0] for i in lines1])

infile2 = open(args.input2)

lines2 = infile2.readlines()

set2 = set([i.strip().split()[0] for i in lines2])

outfile1 = open(args.output1, 'w')

outfile2 = open(args.output2, 'w')

for i in set1:
    if i not in set2:
        outfile1.write(f'{i}\n')

for i in set2:
    if i not in set1:
        outfile2.write(f'{i}\n')