from pymolgen.fragment_mol_combine import *


def test_update_limit():

	fragment_database_mol = get_fragment_database('inputs/fragments10.sdf')

	fragment_database = []

	for i in fragment_database_mol:
		fragment_database.append(i.graph)

	frequencies = get_bond_frequencies('inputs/frequencies10.txt')

	frag_frequencies = get_frag_frequencies('inputs/frag_frequencies10.txt')

	frag_mapping = get_frag_mapping('inputs/fragments10.txt')

	frequencies = update_bond_frequencies(frequencies, frag_mapping)

	new_fragment_database, new_bond_frequencies, new_frag_frequencies, new_frag_mapping = update_limit(2, fragment_database, frequencies, frag_frequencies, frag_mapping)

	assert new_frag_frequencies == [10, 6, 19, 2, 4, 2, 2, 5, 3, 4, 2, 4, 3, 6, 2]

	assert new_frag_mapping == [{0: 0, 1: 1, 2: 2, 3: 3}, {0: 0, 1: 1, 2: 2, 3: 3}, {0: 0, 1: 1, 2: 2}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9}, {0: 0}, {0: 0, 1: 1, 2: 2, 3: 3}, {0: 0, 1: 1}, {0: 0, 1: 1}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9}, {0: 0, 1: 1}, {0: 0}, {0: 0}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9}]

	print('frag_frequencies', frag_frequencies)
	print('new_frag_frequencies', new_frag_frequencies)
	print(new_bond_frequencies)
	print('original frag_mapping', frag_mapping)
	print('new_frag_mapping', new_frag_mapping)
	assert new_bond_frequencies == {(1, 2, 2, 2): 1, (2, 3, 2, 1): 1, (1, 7, 1, 0): 1, (2, 7, 2, 0): 5, (1, 7, 2, 0): 1, (2, 2, 2, 2): 5, (1, 2, 1, 2): 1, (2, 8, 2, 1): 3, (2, 9, 2, 8): 2, (2, 10, 2, 0): 1, (10, 11, 0, 1): 1, (1, 3, 1, 1): 1, (3, 12, 1, 0): 1, (0, 13, 0, 0): 3, (0, 7, 0, 0): 2, (7, 13, 0, 0): 2, (7, 14, 0, 0): 1, (13, 14, 0, 0): 2, (2, 13, 2, 0): 1, (7, 11, 0, 1): 1, (12, 14, 0, 0): 1, (1, 10, 1, 0): 1, (5, 7, 1, 0): 1}


def test_sort_fragments():

	fragment_database_mol = get_fragment_database('inputs/fragments10.sdf')

	fragment_database = []

	for i in fragment_database_mol:
		fragment_database.append(i.graph)

	frequencies = get_bond_frequencies('inputs/frequencies10.txt')

	frag_frequencies = get_frag_frequencies('inputs/frag_frequencies10.txt')

	frag_mapping = get_frag_mapping('inputs/fragments10.txt')

	frequencies = update_bond_frequencies(frequencies, frag_mapping)

	new_fragment_database, new_bond_frequencies, new_frag_frequencies, new_frag_mapping = sort_fragments(fragment_database, frequencies, frag_frequencies, frag_mapping)

	print('new_bond_frequencies', new_bond_frequencies)
	assert new_bond_frequencies == {(1, 16, 0, 0): 1, (3, 16, 1, 2): 1, (0, 3, 2, 2): 1, (0, 29, 2, 0): 1, (29, 33, 0, 1): 1, (0, 33, 2, 3): 1, (0, 10, 2, 1): 1, (10, 27, 1, 3): 1, (6, 33, 0, 9): 1, (6, 33, 0, 11): 1, (11, 38, 1, 0): 1, (37, 38, 0, 2): 1, (36, 38, 4, 6): 1, (6, 38, 0, 7): 1, (1, 35, 0, 0): 1, (13, 35, 1, 15): 1, (13, 34, 1, 7): 1, (3, 34, 2, 11): 1, (3, 4, 1, 0): 1, (0, 4, 2, 0): 5, (3, 4, 2, 0): 1, (0, 0, 2, 2): 5, (0, 32, 2, 5): 1, (0, 3, 2, 1): 1, (0, 31, 2, 2): 1, (0, 35, 2, 22): 1, (0, 9, 2, 1): 3, (0, 35, 2, 26): 1, (0, 35, 2, 31): 1, (0, 7, 2, 8): 2, (0, 35, 2, 35): 1, (0, 12, 2, 0): 1, (5, 12, 1, 0): 1, (1, 28, 0, 0): 1, (1, 28, 0, 3): 1, (26, 28, 3, 4): 1, (26, 39, 10, 0): 1, (0, 39, 2, 3): 1, (25, 26, 6, 12): 1, (7, 25, 8, 10): 1, (3, 10, 1, 1): 1, (8, 10, 0, 1): 1, (3, 15, 2, 5): 1, (7, 15, 8, 11): 1, (1, 2, 0, 0): 3, (2, 17, 0, 0): 1, (8, 17, 0, 0): 1, (5, 17, 1, 5): 2, (2, 24, 0, 0): 1, (3, 24, 2, 3): 1, (3, 19, 1, 0): 1, (0, 19, 2, 7): 1, (0, 18, 2, 3): 1, (2, 24, 0, 6): 1, (1, 4, 0, 0): 2, (2, 4, 0, 0): 2, (4, 14, 0, 0): 1, (2, 14, 0, 0): 2, (0, 2, 2, 0): 1, (4, 5, 0, 1): 1, (0, 20, 2, 8): 1, (1, 21, 0, 0): 1, (13, 21, 1, 2): 1, (13, 22, 1, 2): 1, (2, 22, 0, 9): 1, (8, 14, 0, 0): 1, (3, 40, 2, 5): 1, (3, 12, 1, 0): 1, (12, 23, 0, 1): 1, (4, 23, 0, 4): 1, (4, 11, 0, 1): 1, (0, 30, 2, 4): 1, (6, 40, 0, 8): 1}

	print('new_frag_frequencies', new_frag_frequencies)
	assert new_frag_frequencies == [19, 10, 6, 6, 5, 4, 4, 4, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

	print('new_frag_mapping', new_frag_mapping)
	assert new_frag_mapping == [{0: 0, 1: 1, 2: 2}, {0: 0, 1: 1, 2: 2, 3: 3}, {0: 0}, {0: 0, 1: 1, 2: 2, 3: 3}, {0: 0, 1: 1}, {0: 0, 1: 1}, {0: 0}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10}, {0: 0}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9}, {0: 0, 1: 1, 2: 2, 3: 3}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9}, {0: 0, 1: 1}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10, 11: 11, 12: 12, 13: 13, 14: 14, 15: 15}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10, 11: 11, 12: 12, 13: 13}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10, 11: 11, 12: 12, 13: 13, 14: 14, 15: 15, 16: 16, 17: 17, 18: 18}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10, 11: 11, 12: 12, 13: 13, 14: 14, 15: 15}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10, 11: 11, 12: 12, 13: 13, 14: 14}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10, 11: 11, 12: 12}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10, 11: 11, 12: 12}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4}, {0: 0, 1: 1, 2: 2}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10, 11: 11, 12: 12, 13: 13, 14: 14}, {0: 0, 1: 1, 2: 2}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10, 11: 11}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10, 11: 11}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10, 11: 11, 12: 12, 13: 13, 14: 14, 15: 15, 16: 16, 17: 17, 18: 18, 19: 19, 20: 20, 21: 21, 22: 22, 23: 23, 24: 24, 25: 25, 26: 26, 27: 27, 28: 28, 29: 29, 30: 30, 31: 31, 32: 32, 33: 33, 34: 34, 35: 35, 36: 36, 37: 37, 38: 38, 39: 39, 40: 40, 41: 41, 42: 42, 43: 43, 44: 44, 45: 45, 46: 46, 47: 47, 48: 48, 49: 49, 50: 50, 51: 51, 52: 52}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10, 11: 11, 12: 12, 13: 13}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8}]


def test_filter_ids_1():

	fragment_database_mol = get_fragment_database('../datasets/database1000/fragments1.sdf')

	fragment_database = []

	for i in fragment_database_mol:
		fragment_database.append(i.graph)

	frequencies = get_bond_frequencies('../datasets/database1000/frequencies1.txt')

	frag_frequencies = get_frag_frequencies('../datasets/database1000/frag_frequencies1.txt')

	frag_mapping = get_frag_mapping('../datasets/database1000/fragments1.txt')

	frequencies = update_bond_frequencies(frequencies, frag_mapping)

	filter_ids = [0, 1]

	# new frequencies not containing removed fragments
	frequencies = remove_filter_ids(filter_ids, fragment_database_mol, frequencies, folder='outputs')

	print(frequencies)

	# assert removed fragments not in frequencies
	for key, val in frequencies.items():
		i = key[0]
		j = key[1]

		assert i != 0
		assert j != 1

def test_filter_ids_2():
	# test remove fragment 0 from command line

	os.system('cp ../datasets/database1000/fragments1.sdf outputs/')
	os.system('cp ../datasets/database1000/fragments1.txt outputs/')
	os.system('cp ../datasets/database1000/frequencies1.txt outputs/')
	os.system('cp ../datasets/database1000/frag_frequencies1.txt outputs/')

	os.chdir('outputs')

	subprocess.check_output('python ../../fragment_mol_combine.py -n 1 -i 1 -o 1out --filter_ids 0', shell=True)

	os.chdir('../')
