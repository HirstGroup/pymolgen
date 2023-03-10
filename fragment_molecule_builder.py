import argparse
import copy
import os
import sys

from pymolgen.fragment_molecule import *
from pymolgen.generate import SDFDatasetLargeRAM
from pymolgen.molecule_formats import *
from pymolgen.fragment_builder import bond_frequencies_to_np, get_bond_frequencies, get_fragment_database, get_fragment_bond_frequencies_np

from functools import partial
print = partial(print, flush=True)


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
			f.add_fragment(fragment_id, mol.attach_points)

			if i == fragment_id:

				f.add_fragment(j, fragment_database[j].attach_points)
				f.add_bond(0, 1, k, l)

			elif j == fragment_id:

				f.add_fragment(j, fragment_database[i].attach_points)
				f.add_bond(0, 1, l, k)

			else:
				sys.error('fragmend_id not in bond', bond)


			output_mol_list.append(f)

	return output_mol_list


def extend_molecule_list(FragmentMolecule_list, bond_frequencies, fragment_database_graph, depth=False):

	output_mol_list = []

	for f in FragmentMolecule_list:
		free_valence_list = f.list_free_valence_points()

		for x in range(len(free_valence_list)):

			fragment_id = f.get_frag_id(x)

			for atom in free_valence_list[x]:
				atom_can = fragment_database_graph.fragments[fragment_id].get_canonical_mapping()[atom]
				print(atom, atom_can)
				fragment_bonds, fragment_bond_frequencies = get_fragment_bond_frequencies_np(fragment_id, atom_can, bond_frequencies)
				print(fragment_bonds)
				for bond in fragment_bonds:
					i = bond[0]
					j = bond[1]
					k = bond[2]
					l = bond[3]

					f2 = copy.deepcopy(f)

					# if i corresponds to left fragment j is right fragment
					if i == fragment_id and k == atom_can:
						node_id = f2.add_fragment(j, fragment_database_graph.fragments[j].attachment_points)
						f2.add_bond(x, node_id, atom, l)
						print('new bond82', x, node_id, atom, l)

					# if j corresponds to left fragment i is right framgent
					elif j == fragment_id and l == atom_can:
						node_id = f2.add_fragment(i, fragment_database_graph.fragments[i].attachment_points)
						f2.add_bond(x, node_id, atom, k)
						print('new bond88', x, node_id, atom, k)

					else:
						sys.error('fragmend_id and atom_can not in bond', bond, atom_can)

					if depth is not None:
						total = len(output_mol_list)

						if total % 10000 == 0:
							print(f'DEPTH {depth} TOTAL {total}')

					output_mol_list.append(f2)
					print(f2, f2.list_bonds())

	return output_mol_list

def extend_molecule_list_count(FragmentMolecule_list, bond_frequencies, fragment_database, depth=None):

	total = 0

	for f in FragmentMolecule_list:

		free_valence_list = f.list_free_valence_points()

		for x in range(len(free_valence_list)):

			fragment_id = f.get_frag_id(x)

			for atom in free_valence_list[x]:

				atom_can = atom 

				fragment_bonds, fragment_bond_frequencies = get_fragment_bond_frequencies_np(fragment_id, atom_can, bond_frequencies)

				total += len(fragment_bonds)

				if total % 10000 == 0:
					print(f'DEPTH {depth} TOTAL {total}')

	return total

def extend_molecule_list_depth(FragmentMolecule_list, bond_frequencies, fragment_database_graph, depth):

	for i in range(depth):

		FragmentMolecule_list = extend_molecule_list(FragmentMolecule_list, bond_frequencies, fragment_database_graph, i + 1)

		print(f'FINAL DEPTH {i+1} TOTAL {len(FragmentMolecule_list)}')

	return FragmentMolecule_list

def extend_molecule_list_depth_count(FragmentMolecule_list, bond_frequencies, fragment_database_graph, depth):

	for i in range(depth - 1):

		FragmentMolecule_list = extend_molecule_list(FragmentMolecule_list, bond_frequencies, fragment_database_graph, i + 1)

		print(f'FINAL DEPTH {i+1} TOTAL {len(FragmentMolecule_list)}')

	total = extend_molecule_list_count(FragmentMolecule_list, bond_frequencies, fragment_database)

	print(f'FINAL DEPTH {depth} TOTAL {total}')

	return total

def save_mol_to_sdf(mol, sdffile):

	with open(sdffile, 'a') as f:
		lines = molecule_to_sdf(mol)
		for line in lines:
			f.write(line)
		f.write('$$$$\n')

def save_mol_list_to_sdf(mol_list, sdffile):

	with open(sdffile) as f:
		print('saving to', sdffile)

	for mol in mol_list:
		save_mol_to_sdf(mol, sdffile)


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Build Molecules using the FragmentMolecule class')
	parser.add_argument('-a','--fragments_sdf', help='SDF file of fragments',required=True)
	parser.add_argument('-d','--frequencies_txt', help='Bond frequencies dictionary in txt file',required=True)
	parser.add_argument('--parent_id', type=int, help='Parent id in the fragment database',required=True)
	parser.add_argument('--atom', type=int, help='Atom to build on parent',required=True)
	parser.add_argument('--depth', type=int, help='Depth to build up to',required=True)
	parser.add_argument('--count', action='store_true', help='Count total number of molecules without making them', required=False)
	parser.add_argument('-o','--output', help='Output inchi file name',required=False)

	args = parser.parse_args()

	bond_frequencies = get_bond_frequencies(args.frequencies_txt)
	bond_frequencies = bond_frequencies_to_np(bond_frequencies)

	fragment_database = get_fragment_database(args.fragments_sdf)
	fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

	parent = FragmentMolecule()

	parent.add_fragment(args.parent_id, [args.atom])

	if args.count:
		extend_molecule_list_depth_count([parent], bond_frequencies, fragment_database_graph, args.depth)

	else:
		output_mol_list = extend_molecule_list_depth([parent], bond_frequencies, fragment_database_graph, args.depth)

	if args.output is not None:

		with open(args.output, 'w') as f:
			print('writing to', args.output)

			for j in output_mol_list:
				mol = convert_fragment_molecule_to_mol(j, fragment_database)
				inchi = molecule_to_inchi(mol)
				f.write('%s\n' %inchi)

