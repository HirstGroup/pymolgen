#!/usr/bin/python

import sys,os
import argparse

 
parser = argparse.ArgumentParser(description='This script splits SDF file')
parser.add_argument('-i','--input', help='Input file name', required=True)
parser.add_argument('--n_in', help='Number of Input Files', required=True)
parser.add_argument('-n','--nfiles', help='Number of Files to group', required=True)
parser.add_argument('-o','--output', help='Output file name', required=True)
args = parser.parse_args()

n_in = int(args.n_in)
nfiles = int(args.nfiles)

n = 0
n_file = 0
n_group = 0

for i in range(n_in):

	if n % nfiles == 0:
		outfile = open('%s_%s.sdf' %(args.output, n_group), 'w' )
		n_group += 1

	with open('%s_%s.sdf' %(args.input, i)) as infile:
		for line in infile:
			outfile.write(line)

	n += 1
	