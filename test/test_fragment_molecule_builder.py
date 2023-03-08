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

def test_extend_molecule_list():

	bond_frequencies = get_bond_frequencies('../datasets/database1000/frequencies1.txt')
	bond_frequencies = bond_frequencies_to_np(bond_frequencies)

	fragment_database = get_fragment_database('../datasets/database1000/fragments1.sdf')

	ch3 = FragmentMolecule()

	ch3.add_fragment(0, [0])

	print(len(extend_molecule_list([ch3], bond_frequencies, fragment_database, [])))

def test_extend_molecule_list_2():

	bond_frequencies = get_bond_frequencies('../datasets/database1000/frequencies1.txt')
	bond_frequencies = bond_frequencies_to_np(bond_frequencies)

	fragment_database = get_fragment_database('../datasets/database1000/fragments1.sdf')

	amide = FragmentMolecule()

	amide.add_fragment(2, [1, 2])

	print(len(extend_molecule_list([amide], bond_frequencies, fragment_database, [])))

def test_extend_molecule_list_all():

	bond_frequencies = get_bond_frequencies('../datasets/database1000/frequencies1.txt')
	bond_frequencies = bond_frequencies_to_np(bond_frequencies)

	fragment_database = get_fragment_database('../datasets/database1000/fragments1.sdf')

	answers = [1, 2, 2, 8, 4, 4, 2, 1, 2]

	for i in range(len(fragment_database)):

		mol = fragment_database[i]

		mol2 = FragmentMolecule()

		mol2.add_fragment(i, mol.free_valence_list)

		output = len(extend_molecule_list([mol2], bond_frequencies, fragment_database, []))

		assert output == answers[i]

def test_extend_molecule_list_depth():

	bond_frequencies = get_bond_frequencies('../datasets/database1000/frequencies1.txt')
	bond_frequencies = bond_frequencies_to_np(bond_frequencies)

	fragment_database = get_fragment_database('../datasets/database1000/fragments1.sdf')

	ch3 = FragmentMolecule()

	ch3.add_fragment(0, [0])

	output_mol_list = extend_molecule_list_depth([ch3], bond_frequencies, fragment_database, depth=1)

	for j in output_mol_list:
		assert str(j) == '0-1'

	output_mol_list = extend_molecule_list_depth([ch3], bond_frequencies, fragment_database, depth=2)

	for j in output_mol_list:
		assert str(j) == '0-1-2'		

	output_mol_list = extend_molecule_list_depth([ch3], bond_frequencies, fragment_database, depth=3)

	for j in output_mol_list:
		assert str(j) == '0-1-2-3'

	output_mol_list = extend_molecule_list_depth([ch3], bond_frequencies, fragment_database, depth=4)

	answers = ['0-1-2-3-3', '0-1-2-3-4', '0-1-2-3-5', '0-1-2-3-6']

	for j in range(len(output_mol_list)):
		assert str(output_mol_list[j]) == answers[j]

	output_mol_list = extend_molecule_list_depth([ch3], bond_frequencies, fragment_database, depth=5)

	answers = ['0-1-2-3-3', '0-1-2-3-4', '0-1-2-3-5', '0-1-2-3-6']

	for j in range(len(output_mol_list)):
		print(str(output_mol_list[j]))
		print(output_mol_list[j]._graph.bonds)
		#assert str(output_mol_list[j]) == answers[j]

	with open('test.sdf', 'w') as f:
		print('writing to test.sdf')

	for j in output_mol_list:
		print(j)
		mol = convert_fragment_graph_to_mol(j._graph, fragment_database)
		save_mol_to_sdf(mol, 'test.sdf')

def save_mol_to_sdf(mol, sdffile):

	with open(sdffile, 'a') as f:
		lines = molecule_to_sdf(mol)
		for line in lines:
			f.write(line)
		f.write('$$$$\n')


test_extend_molecule_list_depth()