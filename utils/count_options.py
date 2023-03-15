import argparse

from pymolgen.fragment_builder import *

parser = argparse.ArgumentParser(description='Count total possible number of molecules')
parser.add_argument('-a','--fragments_sdf', help='SDF file of fragments',required=True)
parser.add_argument('--atom', type=int, help='Atom to build from in parent fragment',required=True)
parser.add_argument('-d','--frequencies_txt', help='Bond frequencies dictionary in txt file',required=True)
#parser.add_argument('-p','--parent_file', help='Parent Structure File in SDF format',required=True)
parser.add_argument('-f','--fragments_txt', help='List of fragments in TXT file',required=True)
parser.add_argument('--parent_frag_i', type=int, help='Index i in fragment database of fragment to build from',required=True)

args = parser.parse_args()

#parent_mol = molecule_from_sdf(args.parent_file)

fragment_database = get_fragment_database(args.fragments_sdf)
frag_mapping = get_frag_mapping(args.fragments_txt)
bond_frequencies = get_bond_frequencies(args.frequencies_txt)
bond_frequencies = update_bond_frequencies(bond_frequencies, frag_mapping)
bond_frequencies = bond_frequencies_to_np(bond_frequencies)

parent_frag_i = args.parent_frag_i

atom = args.atom

total = 0

bond_freq_i = get_fragment_bond_frequencies_np(parent_frag_i, atom, bond_frequencies)[0]

count_dict = {}

dead_end_list = []

for x in range(len(bond_freq_i)):

	j = bond_freq_i[x]

	frag_list = []

	i1 = j[0]
	j1 = j[1]
	k1 = j[2]
	l1 = j[3]

	if i1 == parent_frag_i:
		parent_frag_j = j1
		atom_j = l1
	elif j1 == parent_frag_i:
		parent_frag_j = i1
		atom_j = k1
	else:
		sys.error('parent_frag_i not in bond_freq_i') 

	count_dict[parent_frag_j] = 0

	frag_list.append(parent_frag_j)

	#bond_freq_j = get_fragment_bond_frequencies_np(parent_frag_j, atom_j, bond_frequencies)[0]

	j_mol = fragment_database[parent_frag_j]

	j_val = j_mol.free_valence_list
	j_val.remove(atom_j)

	dead_end = True

	for k in j_val:

		canonical_mapping = get_canonical_mapping(j_mol.graph)

		atom_i_can = canonical_mapping[k]

		fragment_bond_frequencies = get_fragment_bond_frequencies_np(parent_frag_j, atom_i_can, bond_frequencies)[0]

		if len(fragment_bond_frequencies) > 0:
			dead_end = False

		if dead_end is True:
			print('fragment %s has free valence points but no bonds' %parent_frag_j)

		total += len(fragment_bond_frequencies)

		count_dict[parent_frag_j] += len(fragment_bond_frequencies)

		if total % 100 == 0: print(total)

	if dead_end is True:
		dead_end_list.append(parent_frag_j)


print(total)

count_dict = dict(sorted(count_dict.items(), key=lambda item: item[1]))

for key, val in count_dict.items():
	print(key, val)

with open('dead_end_list.sdf', 'w') as f:
	for i in dead_end_list:
		mol = fragment_database[i]
		lines = molecule_to_sdf(mol)
		for line in lines:
			f.write(line)
		f.write('$$$$\n')
