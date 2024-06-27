#!/usr/bin/env python

import argparse

from rdkit import Chem

parser = argparse.ArgumentParser(description='Classify molecules as meta, para or meta/para')
parser.add_argument('-i','--input', help='List of molecules as inchis',required=True)
parser.add_argument('-o','--output', help='Output file',required=True)

args = parser.parse_args()

# para
para = Chem.MolFromSmarts("[cH]1[cH]c(-*)[cH][cH]c1-c2c(-[#6])onc2-[#6]")
meta = Chem.MolFromSmarts("[cH]1c(-*)[cH][cH][cH]c1-c2c(-[#6])onc2-[#6]")
parameta = Chem.MolFromSmarts("[cH]1c(-*)c(-*)[cH][cH]c1-c2c(-[#6])onc2-[#6]")

infile = open(args.input)
outfile = open(args.output, 'w')

for line in infile:
	inchi = line.split()[0]

	m = Chem.inchi.MolFromInchi(inchi)

	mol_type = 'NONE'

	if m.HasSubstructMatch(para):
		mol_type = 'para'
	elif m.HasSubstructMatch(meta):
		mol_type = 'meta'
	elif m.HasSubstructMatch(parameta):
		mol_type = 'parameta'

	outfile.write(f'{inchi} {mol_type}\n')


