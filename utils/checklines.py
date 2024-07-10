#!/usr/bin/env python

import sys,os 
import argparse

parser = argparse.ArgumentParser(description='Check that all lines have same number of fields in csv file')
parser.add_argument('-i','--input', help='Input file name',required=True)
parser.add_argument('-o1','--output1', help='Output file with correct csv',required=True)
parser.add_argument('-o2','--output2', help='Output file with incorrect csv',required=True)

args = parser.parse_args()

print(args.input, args.output1, args.output2)

if args.input == args.output1:
	sys.exit('Same input and output1')

if args.input == args.output2:
	sys.exit('Same input and output2')

infile = open(args.input)

outfile1 = open(args.output1, 'w')
outfile2 = open(args.output2, 'w')

first = True

for line in infile:

	if first:
		outfile1.write(line)
		outfile2.write(line)
		first = False
		n_lines = len(line.split(';'))

	else:

		if len(line.split(';')) != n_lines:
			outfile2.write(line)
		else:
			outfile1.write(line)



