from pymolgen.fragment_mol_combine import *

def test_update_limit():

	fragment_database_mol = get_fragment_database('outputs/fragments.sdf')

	fragment_database = []

	for i in fragment_database_mol:
		fragment_database.append(i.graph)

	frequencies = get_bond_frequencies('outputs/frequencies.txt')

	frag_frequencies = get_frag_frequencies('outputs/frag_frequencies.txt')

	frag_mapping = get_frag_mapping('outputs/fragments.txt')

	frequencies = update_bond_frequencies(frequencies, frag_mapping)

	new_fragment_database, new_bond_frequencies, new_frag_frequencies, new_frag_mapping = update_limit(2, fragment_database, frequencies, frag_frequencies, frag_mapping)

	assert new_frag_frequencies == [10, 7, 13, 19, 2, 2, 4, 6, 5, 6, 4, 2, 3, 6, 2]

	assert new_frag_mapping == [{0: 0, 1: 1, 2: 2, 3: 3}, {0: 0, 1: 1}, {0: 0, 1: 1}, {0: 0, 1: 1, 2: 2}, {0: 0, 1: 1, 2: 2}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9}, {0: 0}, {0: 0, 1: 1}, {0: 0, 1: 1}, {0: 0, 1: 1, 2: 2}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9}, {0: 0}, {0: 0}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9}]

	print(frag_frequencies)
	print(new_bond_frequencies)

	assert new_bond_frequencies == {(1, 2, 1, 0): 6, (2, 3, 0, 2): 4, (3, 4, 2, 0): 1, (3, 5, 2, 1): 1, (2, 7, 0, 1): 2, (1, 8, 1, 0): 2, (3, 8, 2, 0): 5, (2, 8, 0, 0): 2, (3, 3, 2, 2): 5, (1, 3, 1, 2): 1, (3, 9, 2, 2): 1, (2, 9, 0, 2): 3, (3, 10, 2, 8): 2, (3, 11, 2, 0): 1, (7, 11, 1, 0): 1, (1, 5, 1, 1): 1, (5, 12, 1, 0): 1, (0, 13, 0, 0): 3, (0, 8, 0, 0): 2, (8, 13, 0, 0): 2, (8, 14, 0, 0): 1, (13, 14, 0, 0): 2, (3, 13, 2, 0): 1, (7, 8, 1, 0): 1, (12, 14, 0, 0): 1, (1, 11, 1, 0): 1, (4, 11, 0, 0): 1, (1, 4, 1, 0): 1}

