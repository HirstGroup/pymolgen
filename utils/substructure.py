import argparse
import os
import sys

from rdkit import Chem


def load_patterns(pattern_list, pattern_format):
	"""
	Load pattern lists into RDKit objects

	Parameters
	----------
	pattern_list : list of str
		list of strings containing patterns in smarts format
	pattern_format : str
		format of pattern_list, inchi, smi or smarts

	Return
	------
	pattern_list_objects : list of RDKit molecules
	"""

	if pattern_format == 'smarts':
		pattern_list_objects = [Chem.MolFromSmarts(x) for x in pattern_list]
	elif pattern_format == 'smi':
		pattern_list_objects = [Chem.MolFromSmiles(x) for x in pattern_list]
	elif pattern_format == 'inchi':
		pattern_list_objects = [Chem.MolFromInchi(x) for x in pattern_list]
	else:
		raise Exception('Only inchi, smi or smarts formats accepted, not', pattern_format)

	return pattern_list_objects


def read_input(input, format_input):
	"""
	Read structures from input file

	Parameters
	----------
	input : str
		name of input file
	format_input : str
		format of input file, smi or inchi

	Returns
	-------
	mol_list : list of RDKit mol objects
	"""

	mol_list = []

	with open(input) as f:
		for line in f:
			x = line.strip('\n').split()[0]

			if format_input == 'smi':
				m = Chem.MolFromSmiles(x)
			elif format_input == 'inchi':
				m = Chem.MolFromInchi(x)
			else:
				raise Exception('Format of input can only be smi or inchi, not', format_input)

			mol_list.append(m)

	return mol_list


def find_substructure(mol_list, pattern_list_objects):
	"""
	Find list of substructures in input molecules

	Parameters
	----------
	mol_list : list of RDKit molecules to search for
	pattern_list_objects : list of RDKit molecules with patterns

	Returns
	-------
	output_false : list of str
		list of inchis without matching pattern
	output_true : list of str
		list of inchis with matching pattern
	"""

	output_false = []
	output_true = []

	for mol in mol_list:

		for pattern in pattern_list_objects:

			if mol.HasSubstructMatch(pattern):

				output_true.append(Chem.MolToInchi(mol))

			else:

				output_false.append(Chem.MolToInchi(mol))

	return output_false, output_true


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Find substructure in input file')

	parser.add_argument('-i','--input', help='Input file name', required=True)
	parser.add_argument('-fi','--format_input', help='Format of input file, inchi or smi', required=True)
	parser.add_argument('-p','--pattern_list', nargs='+', help='Smarts pattern list space separated, use individual quotation marks for each pattern', required=True)
	parser.add_argument('-o1','--output_true', help='Output file name of molecules that contain pattern', required=True)
	parser.add_argument('-o2','--output_false', help='Output file name of molecules that do not contain pattern', required=True)

	args = parser.parse_args()