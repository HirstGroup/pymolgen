#!/usr/bin/python

import sys,os
import argparse

 
parser = argparse.ArgumentParser(description='This script extracts molecules from SDF file')
parser.add_argument('-i','--input', help='Input file name', required=True)
parser.add_argument('-n','--n_mol_list', nargs='+', type=int, help='List of molecules to extract', required=True)
parser.add_argument('-o','--output', help='Output file name', required=True)
args = parser.parse_args()

infile = open(args.input)

outfile = open(args.output, 'w')

n = 1

new = True

for line in infile:
	if n in args.n_mol_list:
		outfile.write(line)
	if '$$$$' in line:
		n += 1
