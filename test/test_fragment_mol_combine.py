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

	print(frag_frequencies)

	new_fragment_database, new_bond_frequencies, new_frag_frequencies = update_limit(2, fragment_database, frequencies, frag_frequencies)

	print(new_frag_frequencies)

	print(frequencies)

	print(new_bond_frequencies)

	check = {(2, 3, 1, 0): 6, (3, 4, 0, 2): 4, (4, 5, 2, 0): 1, (4, 7, 2, 1): 1, (3, 10, 0, 1): 2, (2, 16, 1, 0): 2, (4, 16, 2, 0): 5, (3, 16, 0, 0): 2, (4, 4, 2, 2): 5, (2, 4, 1, 2): 1, (4, 18, 2, 2): 1, (3, 18, 0, 2): 3, (4, 19, 2, 8): 2, (4, 20, 2, 0): 1, (10, 20, 1, 0): 1, (2, 7, 1, 1): 1, (7, 25, 1, 0): 1, (0, 27, 0, 0): 3, (0, 16, 0, 0): 2, (16, 27, 0, 0): 2, (16, 32, 0, 0): 1, (27, 32, 0, 0): 2, (4, 27, 2, 0): 1, (10, 16, 1, 0): 1, (25, 32, 0, 0): 1, (2, 20, 1, 0): 1, (5, 20, 0, 0): 1, (2, 5, 1, 0): 1}

	assert new_bond_frequencies == check

	check = [10, 7, 13, 19, 2, 2, 4, 6, 5, 6, 4, 2, 3, 6, 2]

	assert new_frag_frequencies == check

