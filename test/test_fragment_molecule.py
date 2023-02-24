from pymolgen.fragment_molecule import *

def test_create_fragment_molecule():

	f = FragmentMolecule()

def test_add_fragment():

	f = FragmentMolecule()

	f.add_fragment(0, [0])	


def test_add_fragment2():

	f = FragmentMolecule()

	id = f.add_fragment(10, [0])
	assert id == 0

	id = f.add_fragment(20, [1, 1])
	assert id == 1

def test_add_bond():

	f = FragmentMolecule()

	id = f.add_fragment(10, [0])
	assert id == 0

	id = f.add_fragment(20, [1, 1])
	assert id == 1

	f.add_bond(0, 1, 0, 1)

def test_list_free_valence_points():

	f = FragmentMolecule()

	id = f.add_fragment(10, [0])
	assert id == 0

	id = f.add_fragment(20, [1, 1])
	assert id == 1

	f.add_bond(0, 1, 0, 1)

	assert f.list_free_valence_points() == [[], [1]]

