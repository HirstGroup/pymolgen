import argparse
import copy
import os
import sys

from pymolgen.fragment_molecule import *
from pymolgen.generate import SDFDatasetLargeRAM
from pymolgen.molecule_formats import *
from pymolgen.fragment_builder import bond_frequencies_to_np, get_bond_frequencies, get_fragment_database, get_fragment_bond_frequencies_np

def extend_molecule(fragment_id, bond_frequencies, fragment_database):

	output_mol_list = []

	mol = fragment_database[fragment_id]

	free_valence_list = mol.free_valence_list
	for atom in free_valence_list:

		fragment_bonds, fragment_bond_frequencies = get_fragment_bond_frequencies_np(fragment_id, atom, bond_frequencies)

		for bond in fragment_bonds:

			i = bond[0]
			j = bond[1]
			k = bond[2]
			l = bond[3]

			f = FragmentMolecule()
			f.add_fragment(fragment_id, mol.free_valence_list)

			if i == fragment_id:

				f.add_fragment(j, fragment_database[j].free_valence_list)
				f.add_bond(0, 1, k, l)

			elif j == fragment_id:

				f.add_fragment(j, fragment_database[i].free_valence_list)
				f.add_bond(0, 1, l, k)

			else:
				sys.error('fragmend_id not in bond', bond)


			output_mol_list.append(f)

	return output_mol_list

def extend_molecule_list1(FragmentMolecule_list, bond_frequencies, fragment_database):

	output_mol_list = []

	for f in FragmentMolecule_list:

		free_valence_list = f.list_free_valence_points()

		for x in range(len(free_valence_list)):

			fragment_id = f.get_frag_id(x)

			for atom in free_valence_list[x]:

				fragment_bonds, fragment_bond_frequencies = get_fragment_bond_frequencies_np(fragment_id, atom, bond_frequencies)

				for bond in fragment_bonds:

					i = bond[0]
					j = bond[1]
					k = bond[2]
					l = bond[3]

					f2 = copy.deepcopy(f)

					if i == fragment_id:

						node_id = f2.add_fragment(j, fragment_database[j].free_valence_list)
						f2.add_bond(x, node_id, k, l)

					elif j == fragment_id:

						node_id = f2.add_fragment(i, fragment_database[i].free_valence_list)
						f2.add_bond(x, node_id, l, k)

					else:
						sys.error('fragmend_id not in bond', bond)

					output_mol_list.append(f2)

	return output_mol_list

def extend_molecule_list(FragmentMolecule_list, bond_frequencies, fragment_database, depth=False):

	output_mol_list = []

	for f in FragmentMolecule_list:

		free_valence_list = f.list_free_valence_points()

		for x in range(len(free_valence_list)):

			fragment_id = f.get_frag_id(x)

			for atom in free_valence_list[x]:

				atom_can = f.get_canonical_mapping(x, fragment_database)[atom]

				fragment_bonds, fragment_bond_frequencies = get_fragment_bond_frequencies_np(fragment_id, atom_can, bond_frequencies)

				for bond in fragment_bonds:

					i = bond[0]
					j = bond[1]
					k = bond[2]
					l = bond[3]

					f2 = copy.deepcopy(f)

					# if i corresponds to left fragment j is right fragment
					if i == fragment_id:

						node_id = f2.add_fragment(j, fragment_database[j].free_valence_list)
						f2.add_bond(x, node_id, atom, l)

					# if j corresponds to left fragment i is right framgent
					elif j == fragment_id:

						node_id = f2.add_fragment(i, fragment_database[i].free_valence_list)
						f2.add_bond(x, node_id, atom, k)

					else:
						sys.error('fragmend_id not in bond', bond)

					if depth is not None:
						total = len(output_mol_list)

						if total % 10000 == 0:
							print(f'DEPTH {depth} TOTAL {total}')

					output_mol_list.append(f2)

	return output_mol_list

def extend_molecule_list_depth(FragmentMolecule_list, bond_frequencies, fragment_database, depth):

	for i in range(depth):

		FragmentMolecule_list = extend_molecule_list(FragmentMolecule_list, bond_frequencies, fragment_database, i + 1)

		print(f'FINAL DEPTH {i+1} {len(FragmentMolecule_list)}')

	return FragmentMolecule_list

def save_mol_to_sdf(mol, sdffile):

	with open(sdffile, 'a') as f:
		lines = molecule_to_sdf(mol)
		for line in lines:
			f.write(line)
		f.write('$$$$\n')

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Build Molecules using the FragmentMolecule class')
	parser.add_argument('-a','--fragments_sdf', help='SDF file of fragments',required=True)
	parser.add_argument('-d','--frequencies_txt', help='Bond frequencies dictionary in txt file',required=True)
	parser.add_argument('--parent_id', type=int, help='Parent id in the fragment database',required=True)
	parser.add_argument('--atom', type=int, help='Atom to build on parent',required=True)
	parser.add_argument('--depth', type=int, help='Depth to build up to',required=True)
	parser.add_argument('-o','--output', help='Output inchi file name',required=True)

	args = parser.parse_args()

	bond_frequencies = get_bond_frequencies(args.frequencies_txt)
	bond_frequencies = bond_frequencies_to_np(bond_frequencies)

	fragment_database = get_fragment_database(args.fragments_sdf)

	parent = FragmentMolecule()

	parent.add_fragment(args.parent_id, [args.atom])

	output_mol_list = extend_molecule_list_depth([parent], bond_frequencies, fragment_database, args.depth)

	print(len(output_mol_list))

	with open(args.output, 'w') as f:
		print('writing to', args.output)

		for j in output_mol_list:
			mol = convert_fragment_molecule_to_mol(j, fragment_database)
			inchi = molecule_to_inchi(mol)
			f.write('%s\n' %inchi)

