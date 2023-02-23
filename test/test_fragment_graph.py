from pymolgen.fragment_graph import FragmentGraph

def test_create_fragment_graph():
	f = FragmentGraph()

def test_single_fragment_graph():
	f = FragmentGraph()
	f.add_fragment(0)

	assert 0 in f.fragments