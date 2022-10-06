#!/usr/bin/python

import sys,os
import argparse

parser = argparse.ArgumentParser(description='This script splits SDF file')
parser.add_argument('-i','--input', help='Input file name', required=True)
parser.add_argument('-n','--n_mol', help='Number of Molecules', required=True)
parser.add_argument('-o','--output', help='Output file name', required=True)
args = parser.parse_args()

infile = open(args.input)

n_mol = int(args.n_mol)

name = args.output

n = 0
n_file = 0
new = True

for line in infile:
	if new:
		outfile = open('%s_%s.sdf' %(name, n_file), 'w')
		new = False
		n_file += 1
	outfile.write(line)
	if '$$$$' in line:
		n += 1
		if n % n_mol == 0:
			new = True

outfile.close()

cmd = 'pdftk '

for i in range(n_file):
	print(i)
	cmd += f' {name}_{i}.pdf'
	os.system(f'obabel -xd -isdf {name}_{i}.sdf -O {name}_{i}.svg')
	os.system(f'rsvg-convert -f pdf {name}_{i}.svg -o {name}_{i}.pdf')

cmd += f' cat output {name}_all.pdf'

os.system(cmd)
