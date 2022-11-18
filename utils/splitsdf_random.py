#!/usr/bin/python

import sys,os
import argparse
import random

from pymolgen.generate import SDFDatasetLargeRAM
from pymolgen.molecule_formats import *
 
parser = argparse.ArgumentParser(description='This script splits SDF file')
parser.add_argument('-i','--input', help='Input file name', required=True)
parser.add_argument('-n','--n_mol', type=int, help='Number of Molecules', required=True)
parser.add_argument('-o','--output', help='Output file name', required=True)
parser.add_argument('--batch_size', type=int, help='Batch size, i.e. number of sections input split', required=True)
parser.add_argument('--num_batches', type=int, help='Number of batches, i.e. number of output sets with different random order of molecules', required=True)
args = parser.parse_args()

n_list = list(range(args.n_mol))
print(n_list)

mol_database = SDFDatasetLargeRAM(args.input)

if args.n_mol != len(mol_database):
	sys.exit('N_mol not equal to number of molecules in input file')

if args.batch_size * args.num_batches != args.n_mol:
	sys.exit('batch_size * num_batches != n_mol')

mol_size = int(args.n_mol / (args.batch_size))

for row in range(args.num_batches):

	random.shuffle(n_list)

	for col in range(args.batch_size):

		suffix = row * args.num_batches + col

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

