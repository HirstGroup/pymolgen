import argparse

from pymolgen.fragment_builder import *

parser = argparse.ArgumentParser(description='Count total possible number of molecules')
parser.add_argument('-a','--fragments_sdf', help='SDF file of fragments',required=True)
parser.add_argument('-d','--frequencies_txt', help='Bond frequencies dictionary in txt file',required=True)
parser.add_argument('-f','--fragments_txt', help='List of fragments in TXT file',required=True)
parser.add_argument('-o','--output', help='Output file name',required=True)
parser.add_argument('-p','--parent_file', help='Parent Structure File in SDF format',required=True)
parser.add_argument('--parent_mapping_1', nargs='+', type=int, help='Parent Fragment i dict list space-separated to search fragment database in SDF format',required=True)
parser.add_argument('-r','--remove_hydrogens', type=int, nargs='+', help='Space-separated hydrogen atoms that will be created as attachment points, numbered from 0',required=False)
parser.add_argument('-R','--remove_hydrogens_parent_fragment', type=int, nargs='+', help='Space-separated hydrogen atoms that will be created as attachment points for the parent fragment in database, numbered from 0',required=True)
parser.add_argument('-x','--parent_fragment_file_list', nargs='+', help='Parent Fragment Structure File list space-separated to search fragment database in SDF format',required=True)

args = parser.parse_args()

parent_mapping_1 = args.parent_mapping_1
remove_hydrogens = args.remove_hydrogens
remove_hydrogens_parent_fragment = args.remove_hydrogens_parent_fragment

fragment_database = get_fragment_database(args.fragments_sdf)
frag_mapping = get_frag_mapping(args.fragments_txt)
bond_frequencies = get_bond_frequencies(args.frequencies_txt)
bond_frequencies = update_bond_frequencies(bond_frequencies, frag_mapping)
bond_frequencies = bond_frequencies_to_np(bond_frequencies)

parent_mol = molecule_from_sdf(args.parent_file)

attachment_points = []

new_dict = {}

for i in range(0, len(parent_mapping_1), 2):
    new_dict[parent_mapping_1[i]] = parent_mapping_1[i+1]

parent_mapping_1 = new_dict

# remove hydrogens from parent and determine atoms that will have open valence
for i in remove_hydrogens:
    parent_mol = parent_mol.remove_atom(i)
    for j in parent_mol.free_valence_list:
        if j not in attachment_points:
            attachment_points.append(j)

# make list of equivalent fragments to build on parent
parent_fragment_list = [molecule_from_sdf(x) for x in args.parent_fragment_file_list]

# remove hydrogens from equivalent fragments
for i in range(len(parent_fragment_list)):
    parent_fragment_list[i] = parent_fragment_list[i].remove_atom(remove_hydrogens_parent_fragment[i])

# the original equivalent fragments will be mapped to those in the database to account for the different atom numberings
parent_fragment_original_list = [x for x in parent_fragment_list]

# make a dictionary parent_fragment_i_dict that will map each attachment point to the equivalent fragment id in the database
# make a list parent_fragment_i_list that will contain all equivalent fragments ids
parent_fragment_i_dict = {}
parent_fragment_i_list = []
for i in range(len(parent_fragment_list)):
    j = find_fragment(parent_fragment_list[i], fragment_database)
    print(attachment_points)
    print('id in database of equivalent fragment, j =', j)
    parent_fragment_i_dict[attachment_points[i]] = j
    parent_fragment_i_list.append(j)

    lines = molecule_to_sdf(fragment_database[j])

    with open('parent_fragment%s.sdf' %i, 'w') as outfile:
        for line in lines:
            outfile.write(line)

        outfile.write('$$$$\n')

    if j is False:
        sys.exit('Parent fragment not found')

# make list of all fragments as molecule objects
parent_fragment_list = [fragment_database[x] for x in parent_fragment_i_list]

# map all atoms in each equivalent fragment to the atom numbers in the database
parent_mapping_2 = []
for i in range(len(parent_fragment_list)):
    parent_mapping_2.append(map_mols(parent_fragment_original_list[i].graph, parent_fragment_list[i].graph))

# parent_mapping will map the atoms in the parent with those atom numbers in the equivalent fragments in the database
parent_mapping = {}
n = 0
for key, val in parent_mapping_1.items():
    parent_mapping[key] = parent_mapping_2[n][val]
    n += 1

parent_frag_i = parent_fragment_i_list[0]

atom = parent_fragment_list[0].free_valence_list[0]
atom_parent = 14 #parent_mol.free_valence_list[0]

total = 0

bond_freq_i = get_fragment_bond_frequencies_np(parent_frag_i, atom, bond_frequencies)[0]

if args.output is not None:
	with open(args.output, 'w') as f:
		print('Writing output to', args.output)

all_inchis = []

# loop through all bonds that fragment1 can make
for j in bond_freq_i:

	#frag_list = [parent_frag_i]
	frag_mol_list = [parent_mol]
	frag_bond_list = []
	#frag_free_valence_list = [[]]

	i1 = j[0]
	j1 = j[1]
	k1 = j[2]
	l1 = j[3]

	if i1 == parent_frag_i:
		parent_frag_j = j1
		atom_i = atom_parent
		atom_j = l1
	elif j1 == parent_frag_i:
		parent_frag_j = i1
		atom_i = atom_parent
		atom_j = k1
	else:
		sys.error('parent_frag_i not in bond_freq_i') 

	#frag_list.append(parent_frag_j)
	frag_mol_list.append(fragment_database[parent_frag_j])
	frag_bond_list.append((0,1,atom_i, atom_j))

	mol = combine_all_fragments(frag_mol_list, frag_bond_list)
	mol.hydrogenate()
	inchi = molecule_to_inchi(mol)
	print(inchi)
	#bond_freq_j = get_fragment_bond_frequencies_np(parent_frag_j, atom_j, bond_frequencies)[0]

	j_mol = fragment_database[parent_frag_j]

	j_val = j_mol.free_valence_list
	j_val.remove(atom_j)
	print('parent_frag_j =', parent_frag_j)
	print('j_val =', j_val)
	# loop through all attachment points of fragment2
	for k in j_val:

		canonical_mapping = get_canonical_mapping(j_mol.graph)

		atom_i_can = canonical_mapping[k]

		fragment_bond_frequencies = get_fragment_bond_frequencies_np(parent_frag_j, atom_i_can, bond_frequencies)[0]


		for l in fragment_bond_frequencies:

			total += 1

			i2 = l[0]
			j2 = l[1]
			k2 = l[2]
			l2 = l[3]					

			if i2 == parent_frag_j:
				parent_frag_l = j2
				atom_i2 = k2
				atom_j2 = l2
			elif j2 == parent_frag_j:
				parent_frag_l = i2
				atom_i2 = l2
				atom_j2 = k2
			else:
				sys.error('parent_frag_j not in bond_freq_l')

			#frag_list2 = [x for x in frag_list]
			frag_mol_list2 = [x for x in frag_mol_list]
			frag_bond_list2 = [x for x in frag_bond_list]

			#frag_list2.append(parent_frag_l)
			frag_mol_list2.append(fragment_database[parent_frag_l])
			frag_bond_list2.append((1,2,k, atom_j2))			 

			mol = combine_all_fragments(frag_mol_list2, frag_bond_list2)
			mol.hydrogenate()
			inchi = molecule_to_inchi(mol)
			print(inchi)
			all_inchis.append(inchi)

			if total % 100 == 0: 
				print(total)

				if args.output is not None:
					with open(args.output, 'a') as f:
						for inchi in all_inchis:
							f.write('%s\n' %inchi)
					all_inchis = []


if args.output is not None:
	with open(args.output, 'a') as f:
		for inchi in all_inchis:
			f.write('%s\n' %inchi)

print(total)


