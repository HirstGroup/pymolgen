from pymolgen.fragment_builder import bond_frequencies_to_np, get_bond_frequencies, get_fragment_database, get_fragment_bond_frequencies_np
from pymolgen.fragment_graph import *

def test_create_fragment_graph():
	f = FragmentGraph()


def test_single_fragment_graph():
	f = FragmentGraph()
	f.add_fragment(10, [1])
	assert len(f.fragments) == 1


def test_add_bond():
	f = FragmentGraph()
	f.add_fragment(0, [0])
	f.add_fragment(1, [0])
	f.add_bond(0,1,0,0)
	assert f.bonds == [(0,1,0,0)]

	assert f.free_valence_points == [[], []]


def test_add_bond2():
	f = FragmentGraph()
	f.add_fragment(0, [0])
	f.add_fragment(1, [0,1])
	f.add_bond(0,1,0,0)
	assert f.bonds == [(0,1,0,0)]

	assert f.free_valence_points == [[], [1]]


def test_add_invalid_bond():
	f = FragmentGraph()
	try:
		f.add_bond(0, 1, 0, 0)
	except AssertionError:
		return
	assert False


def test_add_invalid_bond_2():
	f = FragmentGraph()
	f.add_fragment(0, [1])
	f.add_fragment(1, [1])
	try:
		# Should fail because fragment 1 doesn't have 1 as an attachment point
		f.add_bond(0, 1, 0, 1)
	except AssertionError:
		return
	assert False


def test_add_invalid_bond_3():
	f = FragmentGraph()
	f.add_fragment(0, [1])
	f.add_fragment(1, [0])
	f.add_fragment(2, [1])
		
	f.add_bond(0, 1, 1, 0)

	try:
		# Should fail because fragment 1 has used up 0 as attachment point
		f.add_bond(1,2,0,1)
	except AssertionError:
		return
	assert False


def test_graph_build_probability():

	f = FragmentGraph()
	f.add_fragment(0, [1])
	f.add_fragment(1, [1, 1])
	f.add_bond(0, 1, 1, 1, 0.5)

	assert f.build_probability == 0.5

	f.add_fragment(2, [1])
	f.add_bond(1, 2, 1, 1, 0.5)

	assert f.build_probability == 0.25


def test_graph_build_probability2():

	f = FragmentGraph(build_probability2=1.0)
	f.add_fragment(0, [1])
	f.add_fragment(1, [1, 1])
	f.add_bond(0, 1, 1, 1, 0.5, 0.5)

	assert f.build_probability == 0.5
	assert f.build_probability2 == 0.5

	f.add_fragment(2, [1])
	f.add_bond(1, 2, 1, 1, 0.5, 0.5)

	assert f.build_probability == 0.25
	assert f.build_probability2 == 0.25


def test_molecular_weight_node():

	n = FragmentGraphNode([0,1], 100.0)

	assert n.molecular_weight == 100.0


def test_molecular_weight_graph():

	f = FragmentGraph()

	f.add_fragment(0, [0,1], {0:0}, 100.0)

	assert f.molecular_weight == 100.0

	f.add_fragment(1, [2,3], {1:1}, 100.0)

	assert f.molecular_weight == 200.0


def test_convert_fragment_database_to_graph():

    fragment_database = get_fragment_database('../datasets/database1000/fragments10.sdf')
    fragment_database_graph = convert_fragment_database_to_graph(fragment_database)   

    print(fragment_database_graph.molecular_weight)

    assert (fragment_database_graph.molecular_weight - 3156.9331529999995)**2 < 0.01

    for frag in fragment_database_graph.fragments.values():
    	print(frag.molecular_weight)

    mw_list = [x.molecular_weight for x in fragment_database_graph.fragments.values()]

    print(mw_list)

    assert mw_list == [15.03491, 67.04707, 43.02507, 14.02694, 64.0588, 113.11860999999999, 76.09787999999999, 67.07031, 18.998403, 45.01777, 141.12900999999997, 41.07285, 78.09357999999999, 453.5318100000001, 28.010399999999997, 69.10649, 13.01897, 58.06298, 16.02264, 44.03304, 77.10584999999999, 76.09787999999999, 17.007369999999998, 66.03909999999999, 127.12530999999997, 66.06233999999999, 85.10588999999999, 35.453, 121.14188999999998, 15.9994, 154.12474, 75.08990999999999, 109.17136999999998, 67.06730999999999, 76.09787999999999, 84.14139999999999, 83.08994999999999, 83.13342999999999, 79.07347, 116.14251999999998, 77.08560999999999]