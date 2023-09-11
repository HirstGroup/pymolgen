#!/usr/bin/python

import sys,os
import argparse

 
parser = argparse.ArgumentParser(description='This script splits SDF file')
parser.add_argument('-i','--input', help='Input file name', required=True)
parser.add_argument('-n','--n_mol', help='Number of Molecules per SDF file', required=True)
parser.add_argument('-o','--output', help='Output file name', required=True)
args = parser.parse_args()

infile = open(args.input)

n_mol = int(args.n_mol)

n = 0
n_file = 0
new = True

for line in infile:
	if new:
		outfile = open('%s_%s.sdf' %(args.output, n_file), 'w')
		new = False
	outfile.write(line)
	if '$$$$' in line:
		n += 1
		if n % n_mol == 0:
			n_file += 1
			new = True