#!/usr/bin/env python

import argparse

from rdkit import Chem

para = Chem.MolFromSmarts("[cH]1[cH]c(-*)[cH][cH]c1-c2c(-[#6])onc2-[#6]")
meta = Chem.MolFromSmarts("[cH]1c(-*)[cH][cH][cH]c1-c2c(-[#6])onc2-[#6]")
parameta = Chem.MolFromSmarts("[cH]1c(-*)c(-*)[cH][cH]c1-c2c(-[#6])onc2-[#6]")


def classify(mol):
	"""
	Classify RDKit molecule into meta, para or meta/para

	Parameters
	----------
	mol : RDKit molecule object

	Returns
	-------
	mol_type : str
		classification of molecule into meta, para or parameta

	"""

	mol_type = 'NONE'

	if mol.HasSubstructMatch(para):
		mol_type = 'para'
	elif mol.HasSubstructMatch(meta):
		mol_type = 'meta'
	elif mol.HasSubstructMatch(parameta):
		mol_type = 'parameta'

	return mol_type


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Classify molecules as meta, para or meta/para')
	
	parser.add_argument('-i','--input', help='List of molecules as inchis',required=True)
	parser.add_argument('-o','--output', help='Output file',required=True)

	args = parser.parse_args()

	infile = open(args.input)
	outfile = open(args.output, 'w')

	for line in infile:
		inchi = line.split()[0]

		mol = Chem.inchi.MolFromInchi(inchi)

		mol_type = classify(mol)

		outfile.write(f'{inchi} {mol_type}\n')