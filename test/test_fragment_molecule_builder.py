import os
import sys

from pymolgen.fragment_molecule_builder import *

def test_fragment_molecule_builder():

	f = FragmentMolecule()

	f.add_fragment(10, [0,1]) # 0
	f.add_fragment(20, [2,2]) # 1
	f.add_fragment(20, [2,2]) # 2
	f.add_fragment(30, [4])   # 3
	f.add_fragment(40, [5])   # 4

	f.add_bond(0, 1, 0, 2)
	assert f.list_free_valence_points() == [[1], [2], [2, 2], [4], [5]]

	f.add_bond(1, 2, 2, 2)
	assert f.list_free_valence_points() == [[1], [], [2], [4], [5]]

	f.add_bond(2, 3, 2, 4)
	assert f.list_free_valence_points() == [[1], [], [], [], [5]]

	f.add_bond(0, 4, 1, 5)
	assert f.list_free_valence_points() == [[], [], [], [], []]

def test_extend_molecule():

	bond_frequencies = get_bond_frequencies('../datasets/database1000/frequencies1.txt')
	bond_frequencies = bond_frequencies_to_np(bond_frequencies)

	fragment_database = get_fragment_database('../datasets/database1000/fragments1.sdf')

	answers = [1, 2, 2, 8, 4, 4, 2, 1, 2]

	for fragment_id in range(9):

		output_mol_list = extend_molecule(fragment_id, bond_frequencies, fragment_database)

		assert len(output_mol_list) == answers[fragment_id]

