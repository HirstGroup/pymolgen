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

def test_total_free_valence():

	f = FragmentMolecule()
	f.add_fragment(0, [0, 1])	
	f.add_fragment(0, [0, 1])

	assert f.get_total_free_valence() == 4

def test_convert_to_networkx():

	f = FragmentMolecule()

	f.add_fragment(0, [0])
	f.add_fragment(10, [1,2,3])
	f.add_fragment(20, [5, 6])
	f.add_fragment(40, [7, 8])
	f.add_fragment(30, [9])
	f.add_fragment(50, [4])

	f.add_bond(0, 1, 0, 1)
	f.add_bond(1, 2, 2, 5)
	f.add_bond(2, 3, 6, 7)
	f.add_bond(3, 4, 8, 9)
	f.add_bond(1, 5, 3, 4)

	g = f._graph.convert_to_networkx()

	frag_id_list = [0, 10, 20, 40, 30, 50]

	for i in range(len(g.nodes)):
		assert frag_id_list[i] == g.nodes[i]['frag_id']

	answers = [{0: 0, 10: 1}, {0: 0, 10: 1}, {10: 2, 20: 5}, {10: 3, 50: 4}, {10: 2, 20: 5}, {20: 6, 40: 7}, {20: 6, 40: 7}, {30: 9, 40: 8}, {30: 9, 40: 8}, {10: 3, 50: 4}]

	counter = 0
	for i in range(len(g.nodes)):
		for j in g[i]:
			assert g[i][j]['atoms'] == answers[counter]
			counter += 1

def test_make_canonical():

	f = FragmentMolecule()

	f.add_fragment(0, [0])
	f.add_fragment(10, [1,2,3])
	f.add_fragment(20, [5, 6])
	f.add_fragment(40, [7, 8])
	f.add_fragment(30, [9])
	f.add_fragment(50, [4])

	f.add_bond(0, 1, 0, 1)
	f.add_bond(1, 2, 2, 5)
	f.add_bond(2, 3, 6, 7)
	f.add_bond(3, 4, 8, 9)
	f.add_bond(1, 5, 3, 4)


def test_cap():

	f = FragmentMolecule()

	f.add_fragment(0, [0])
	f.add_fragment(10, [1,2,3,40])
	f.add_fragment(20, [5, 6, 70])
	f.add_fragment(40, [7, 8, 90])
	f.add_fragment(30, [9])
	f.add_fragment(50, [4])

	f.add_bond(0, 1, 0, 1)
	f.add_bond(1, 2, 2, 5)
	f.add_bond(2, 3, 6, 7)
	f.add_bond(3, 4, 8, 9)
	f.add_bond(1, 5, 3, 4)

	f2 = f.cap()

	assert f2.list_free_valence_points() == [[], [], [], [], [], [], [], [], []]

	assert f2._graph.bonds == [(0, 1, 0, 1), (1, 2, 2, 5), (2, 3, 6, 7), (3, 4, 8, 9), (1, 5, 3, 4), (1, 6, 40, 0), (2, 7, 70, 0), (3, 8, 90, 0)]
	assert f2.list_frag_id() == [0, 10, 20, 40, 30, 50, -1, -1, -1]


def test_convert_to_networkx2():

	f = FragmentMolecule()

	f.add_fragment(0, [0])
	f.add_fragment(10, [1,2,3,40])
	f.add_fragment(20, [5, 6, 70])
	f.add_fragment(40, [7, 8, 90])
	f.add_fragment(30, [9])
	f.add_fragment(50, [4])

	f.add_bond(0, 1, 0, 1)
	f.add_bond(1, 2, 2, 5)
	f.add_bond(2, 3, 6, 7)
	f.add_bond(3, 4, 8, 9)
	f.add_bond(1, 5, 3, 4)

	f2 = f.cap()

	g = f2._graph.convert_to_networkx()

	frag_id_list = [0, 10, 20, 40, 30, 50, -1, -1, -1]

	for i in range(len(g.nodes)):
		assert frag_id_list[i] == g.nodes[i]['frag_id']

	answers = [{0: 0, 10: 1}, {0: 0, 10: 1}, {10: 2, 20: 5}, {10: 3, 50: 4}, {-1: 0, 10: 40}, {10: 2, 20: 5}, {20: 6, 40: 7}, {-1: 0, 20: 70}, {20: 6, 40: 7}, {30: 9, 40: 8}, {-1: 0, 40: 90}, {30: 9, 40: 8}, {10: 3, 50: 4}, {-1: 0, 10: 40}, {-1: 0, 20: 70}, {-1: 0, 40: 90}]

	counter = 0
	for i in range(len(g.nodes)):
		for j in g[i]:
			assert g[i][j]['atoms'] == answers[counter]
			counter += 1

def test_get_hash():

	f = FragmentMolecule()

	f.add_fragment(0, [0])
	f.add_fragment(10, [1,2,3,40])
	f.add_fragment(20, [5, 6, 70])
	f.add_fragment(40, [7, 8, 90])
	f.add_fragment(30, [9])
	f.add_fragment(50, [4])

	f.add_bond(0, 1, 0, 1)
	f.add_bond(1, 2, 2, 5)
	f.add_bond(2, 3, 6, 7)
	f.add_bond(3, 4, 8, 9)
	f.add_bond(1, 5, 3, 4)

	graph_hash = f.get_hash()

	assert graph_hash == '6e9e838860c4d292c41d5e0b80a91917'

test_get_hash()
