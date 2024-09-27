from pymolgen.analysis.fragment_database_equivalent import *

def test_init():
	
    fragment_database = get_fragment_database('../../datasets/database1000/fragments10.sdf')
    fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    fragment_database_equivalent = FragmentDatabaseEquivalent(fragment_database_graph)

    print(fragment_database_equivalent.fragments)

test_init()	