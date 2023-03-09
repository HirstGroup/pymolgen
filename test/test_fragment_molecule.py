from pymolgen.fragment_mol import print_fragments, get_canonical_mapping, map_mols, get_frag_mapping, update_bond_frequencies
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

def test_networkx_graph():

	f = FragmentMolecule()
	f.add_fragment(0, [0])	
	fragment_database = get_fragment_database('../datasets/database1000/fragments1.sdf')

	assert get_canonical_mapping(f._graph.fragments[0].get_molecule(fragment_database).graph) == {0: 0, 1: 1, 2: 1, 3: 1}

def test_networkx_graph2():

	f = FragmentMolecule()
	f.add_fragment(0, [0])	
	fragment_database = get_fragment_database('../datasets/database1000/fragments1.sdf')

	assert f.get_canonical_mapping(0, fragment_database) == {0: 0, 1: 1, 2: 1, 3: 1}

test_networkx_graph2()