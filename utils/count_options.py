import argparse

from pymolgen.fragment_builder import *

parser = argparse.ArgumentParser(description='Count total possible number of molecules')
parser.add_argument('-a','--fragments_sdf', help='SDF file of fragments',required=True)
parser.add_argument('-d','--frequencies_txt', help='Bond frequencies dictionary in txt file',required=True)
#parser.add_argument('-p','--parent_file', help='Parent Structure File in SDF format',required=True)
parser.add_argument('-f','--fragments_txt', help='List of fragments in TXT file',required=True)

args = parser.parse_args()

#parent_mol = molecule_from_sdf(args.parent_file)

fragment_database = get_fragment_database(args.fragments_sdf)
frag_mapping = get_frag_mapping(args.fragments_txt)
bond_frequencies = get_bond_frequencies(args.frequencies_txt)
bond_frequencies = update_bond_frequencies(bond_frequencies, frag_mapping)
bond_frequencies = bond_frequencies_to_np(bond_frequencies)

parent_frag_i = 14

atom = 2

total = 0

bond_freq_i = get_fragment_bond_frequencies_np(parent_frag_i, atom, bond_frequencies)[0]

print('parent bond_freq_i =', len(bond_freq_i))

for j in bond_freq_i:

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

	bond_freq_j = get_fragment_bond_frequencies_np(parent_frag_j, atom_j, bond_frequencies)[0]

	j_mol = fragment_database[parent_frag_j]

	j_val = j_mol.free_valence_list
	j_val.remove(atom_j)

	for k in j_val:

		canonical_mapping = get_canonical_mapping(j_mol.graph)

		atom_i_can = canonical_mapping[k]

		fragment_bond_frequencies = get_fragment_bond_frequencies_np(parent_frag_j, atom_i_can, bond_frequencies)[0]

		total += len(fragment_bond_frequencies)

		if total % 100 == 0: print(total)


print(total)


