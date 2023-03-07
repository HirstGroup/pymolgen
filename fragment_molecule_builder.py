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