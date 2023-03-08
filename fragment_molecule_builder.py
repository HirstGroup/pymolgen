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

def extend_molecule_list(FragmentMolecule_list, bond_frequencies, fragment_database):

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

def extend_molecule_list_depth(FragmentMolecule_list, bond_frequencies, fragment_database, depth):

	for i in range(depth):

		FragmentMolecule_list = extend_molecule_list(FragmentMolecule_list, bond_frequencies, fragment_database)

	return FragmentMolecule_list




