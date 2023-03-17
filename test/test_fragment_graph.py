from pymolgen.fragment_builder import bond_frequencies_to_np, get_bond_frequencies, get_fragment_database, get_fragment_bond_frequencies_np
from pymolgen.fragment_graph import FragmentGraph

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

test_graph_build_probability()