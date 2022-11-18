#!/usr/bin/python

import sys,os
import argparse
import random

from pymolgen.generate import SDFDatasetLargeRAM
from pymolgen.molecule_formats import *
 
parser = argparse.ArgumentParser(description='This script splits SDF file')
parser.add_argument('-i','--input', help='Input file name', required=True)
parser.add_argument('-o','--output', help='Output file name', required=True)
parser.add_argument('-n','--n_batch', type=int, help='N batch number', required=True)
parser.add_argument('--batch_size', type=int, help='Batch size', required=True)
args = parser.parse_args()

mol_database = SDFDatasetLargeRAM(args.input)

n_mol = len(mol_database)

mol_size = int(n_mol / (args.batch_size))

row = args.n_batch

n_list = list(range(n_mol))

random.shuffle(n_list)

for col in range(args.batch_size):

	suffix = row * args.n_batch + col

	first_mol = mol_size * col
	last_mol = mol_size * (col + 1)

	print(row, col, suffix, first_mol, last_mol)

	outfile_name = args.output + str(suffix) + '.sdf'

	with open(outfile_name, 'w') as outfile:
		print('Writing to', outfile_name )

	for n in n_list[first_mol:last_mol]:

		mol = mol_database[n]

		lines = molecule_to_sdf(mol)

		with open(outfile_name, 'a') as outfile:
			for line in lines:
				outfile.write(line)
			outfile.write('$$$$\n')

