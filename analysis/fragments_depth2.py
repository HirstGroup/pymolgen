import argparse

from pymolgen.fragment_molecule_builder import *

parser = argparse.ArgumentParser(description='Build Molecules using the FragmentMolecule class')
parser.add_argument('-a','--fragments_sdf', help='SDF file of fragments',required=True)
parser.add_argument('-d','--frequencies_txt', help='Bond frequencies dictionary in txt file',required=True)
parser.add_argument('--parent_id', type=int, help='Parent id in the fragment database',required=True)
parser.add_argument('--atom', type=int, help='Atom to build on parent',required=True)
parser.add_argument('--depth', type=int, help='Depth to build up to',required=True)
parser.add_argument('-o','--output', help='Output inchi file name',required=False)

args = parser.parse_args()

bond_frequencies = get_bond_frequencies(args.frequencies_txt)
bond_frequencies = bond_frequencies_to_np(bond_frequencies)

fragment_database = get_fragment_database(args.fragments_sdf)
#fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

parent = FragmentMolecule()

parent.add_fragment(args.parent_id, [args.atom])

fragment_bonds, fragment_bond_frequencies = get_fragment_bond_frequencies_np(args.parent_id, args.atom, bond_frequencies)

print(fragment_bonds)

frag_a = []

for bond in fragment_bonds:

	i = bond[0]
	j = bond[1]
	k = bond[2]
	l = bond[3]

	if i == args.parent_id:
		fragment_from = i
		fragment_to = j
		attach_from = k
		attach_to = l
		frag_a.append(j)
	elif j == args.parent_id:
		frag_a.append(i)
		fragment_from = j
		fragment_to = i
		attach_from = l
		attach_to = k
	else:
		print('parent id not in bond', bond)

	fragment_bonds_a, fragment_bond_frequencies_a = get_fragment_bond_frequencies_np(, args.atom, bond_frequencies)
