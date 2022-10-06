import sys,os
import numpy as np

import networkx

from pymolgen.fragment_mol import *
from pymolgen.fragment_builder import *

from functools import partial
print = partial(print, flush=True)

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def combine_fragment_databases(fragment_database, frequencies, frag_frequencies, frag_mapping, fragments_sdf_2, fragments_txt_2, frequencies_txt_2, frag_frequencies_txt_2, fragments_sdf_out, fragments_txt_out, frequencies_txt_out, frag_frequencies_txt_out, limit=None):

    fragment_database_mol2 = get_fragment_database(fragments_sdf_2)

    fragment_database2 = []

    for i in fragment_database_mol2:
        fragment_database2.append(i.graph)

    frequencies2 = get_bond_frequencies(frequencies_txt_2)

    frag_mapping2 = get_frag_mapping(fragments_txt_2)

    frequencies2 = update_bond_frequencies(frequencies2, frag_mapping2)

    frag_frequencies2 = get_frag_frequencies(frag_frequencies_txt_2)

    if limit is not None:
        print('Before limit ', len(fragment_database2))
        fragment_database2, frequencies2, frag_frequencies2, frag_mapping2 = update_limit(limit, fragment_database2, frequencies2, frag_frequencies2, frag_mapping2)
        print('After limit ', len(fragment_database2))     

    #mapping of fragment atom indeces from 2 to 1 (or 2 to 2 if new fragment)
    frag_mapping2to1 = []

    #map of fragment index from 2 to final database
    frag_index_mapping = []

    fragment_database2_len = len(fragment_database2)

    for i in range(fragment_database2_len):

        if i % (fragment_database2_len // 10) == 0:
            print(i)

        fragment = fragment_database2[i]

        frag1_is_new, frag1_index, frag1_map = get_fragment_index(fragment, fragment_database)

        if frag1_is_new: 

            frag_frequencies.append(frag_frequencies2[i])

            frag_index_mapping.append(len(fragment_database))

            fragment_database.append(fragment)

        else:
            frag_frequencies[frag1_index] += frag_frequencies2[i]

            frag_index_mapping.append(frag1_index)

        frag_mapping2to1.append(frag1_map)

    for key, val in frequencies2.items():

        frag1_index = key[0]
        frag2_index = key[1]
        frag1_bond = key[2]
        frag2_bond = key[3]

        #get mapping for atom numbers
        frag1_map = frag_mapping2to1[frag1_index]
        frag2_map = frag_mapping2to1[frag2_index]

        #convert frag indeces 
        frag1_index = frag_index_mapping[key[0]]
        frag2_index = frag_index_mapping[key[1]]

        update_freq(frequencies, frag1_index, frag2_index, frag1_map, frag2_map, frag1_bond, frag2_bond, val)

    save_frequencies_txt(frequencies, frequencies_txt_out)

    save_fragments_sdf(fragment_database, fragments_sdf_out)

    save_frag_frequencies_txt(frag_frequencies, frag_frequencies_txt_out)

    save_fragments_txt(fragment_database, fragments_txt_out)

"""
cycles = networkx.cycle_basis(mol.graph)

        for cycle in cycles:
            if len(cycle) > 10:

                with open('cycles.sdf', 'a') as outfile:
                    lines = molecule_to_sdf(mol)

                    for line in lines:
                        outfile.write(line)

                    outfile.write('$$$$\n')
            continue
"""

def remove_bond_frequencies(bond_frequencies, fragment_list):
    """
    Remove bond frequencies when removing framgents from a database)

    fragment_list: list of fragments to keep
    """

    fragment_list = set(fragment_list)

    d = {}

    for key, val in bond_frequencies.items():

        i = key[0]
        j = key[1]

        k = key[2]
        l = key[3]

        if i in fragment_list and j in fragment_list:

            d[(i,j,k,l)] = val

    return d

def list_elements(mol):

    elements = set()

    for i in mol.graph.nodes:
        elements.add(mol.graph.nodes[i]["element"])

    return elements

def filter_database(fragment_database_mol, inchi_filter=None, pains=False):

    elements = set()

    filter_list = []

    with open('filter.sdf', 'w') as outfile:
        print('Writing to filter.sdf')

    with open('aliphatic3.sdf', 'w') as outfile:
        print('Writing to aliphatic3.sdf')

    with open('aliphatic4.sdf', 'w') as outfile:
        print('Writing to aliphatic4.sdf')

    with open('cycles.sdf', 'w') as outfile:
        print('Writing to cycles.sdf')

    with open('cages3.sdf', 'w') as outfile:
        print('Writing to cages3.sdf')

    with open('cages4.sdf', 'w') as outfile:
        print('Writing to cages4.sdf')

    with open('sulfone.sdf', 'w') as outfile:
        print('Writing to sulfone.sdf')

    with open('sulfur_out.sdf', 'w') as outfile:
        print('Writing to sulfone.sdf')

    with open('sulfur_cyclic.sdf', 'w') as outfile:
        print('Writing to sulfur_cyclic.sdf')

    with open('check.sdf', 'w') as outfile:
        print('Writing to check.sdf')

    if inchi_filter is not None:
        inchi_list = set()

        with open(inchi_filter) as infile:
            for line in infile:
                inchi_list.add(line.strip('\n'))

    if pains is True:

        from openeye import oechem
        from pymolgen.properties_pymolgen import gen_pains_database, pains_filter

        pains_database = gen_pains_database()

        with open('pains.sdf', 'w') as outfile:
            print('Writing to pains.sdf')

    with open('thioether.sdf', 'w') as outfile:
        print('Writing to thioether.sdf')

    for i in range(len(fragment_database_mol)):

        mol = fragment_database_mol[i]

        mol_elements = list_elements(mol)

        elements.update(mol_elements)

        allowed_elements = {'Cl', 'O', 'C', 'F', 'N', 'S', 'H'}

        allowed = True

        for j in mol_elements:
            if j not in allowed_elements:
                allowed = False

        if allowed is False:
            continue

        mol_h = mol.copy()
        mol_h.hydrogenate()

        inchi = molecule_to_inchi(mol_h)
        smi = molecule_to_smiles(mol_h)

        if inchi_filter is not None:
            if inchi in inchi_list:
                print('INCHI FILTER', inchi)
                continue

        if pains is True:

            oemol = oechem.OEGraphMol()
            oechem.OESmilesToMol(oemol, smi)

            oechem.OEAddExplicitHydrogens(oemol)

            if pains_filter(oemol, pains_database) is False:

                save_mol_to_sdf('pains.sdf', mol)

                continue

        rdmol = Chem.MolFromInchi(inchi)

        if rdmol is None:
            continue

        if rdmol is not None:
            ri = rdmol.GetRingInfo()
            largest_ring_size = max((len(r) for r in ri.AtomRings()), default=0)
            if largest_ring_size > 8:        

                save_mol_to_sdf('cycles.sdf', mol)

                continue

        n = rdMolDescriptors.CalcNumAliphaticRings(rdmol)

        if n == 3:

            if rdMolDescriptors.CalcNumHeavyAtoms(rdmol) < 12:

                save_mol_to_sdf('cages3.sdf', mol)

            else:

                save_mol_to_sdf('aliphatic3.sdf', mol)

                continue

        if n == 4:

            if rdMolDescriptors.CalcNumHeavyAtoms(rdmol) < 16:

                save_mol_to_sdf('cages4.sdf', mol)

                continue

            save_mol_to_sdf('aliphatic4.sdf', mol)

            continue

        if n > 4:

            save_mol_to_sdf('filter.sdf', mol)

            continue

        if is_sulfur(mol):

            if is_thioether(mol):
                save_mol_to_sdf('thioether.sdf', mol)
                continue

            if is_sulfone(mol):
                save_mol_to_sdf('sulfone.sdf', mol)

            else:

                if is_cyclic_sulfur(mol):
                    save_mol_to_sdf('sulfur_cyclic.sdf', mol)

                else:

                    save_mol_to_sdf('sulfur_out.sdf', mol)

                    continue

        filter_list.append(i)

        save_mol_to_sdf('check.sdf', mol) 

    print(elements)

    return filter_list

def save_mol_to_sdf(outfile_name, mol):

    with open(outfile_name, 'a') as outfile:
        lines = molecule_to_sdf(mol)

        for line in lines:
            outfile.write(line)

        outfile.write('$$$$\n')    

def is_sulfur(mol):

    for i in mol.graph.nodes:
        if mol.graph.nodes[i]["element"] == 'S':
            return True


def is_cyclic_sulfur(mol):

    for i in mol.graph.nodes:
        if mol.graph.nodes[i]["element"] == 'S':
            if mol.is_cyclic(i):
                return True
            else:
                for j in mol.graph[i]:
                    if mol.is_cyclic(j):
                        return True


    return False

def is_sulfone(mol):

    for i in mol.graph.nodes:
        if mol.graph.nodes[i]["element"] == 'S':
            o_count = 0
            for j in mol.graph[i]:
                if mol.graph.nodes[j]["element"] == 'O':
                    o_count += 1
            if o_count == 2:
                return True

    return False


def is_thioether(mol):

    for i in mol.graph.nodes:
        if mol.graph.nodes[i]["element"] == 'S':
            for j in mol.graph[i]:
                if mol.graph.nodes[j]["element"] == 'S':
                    return True

    return False

def loop(n, fragments_sdf_in, fragments_txt_in, frequencies_txt_in, frag_frequencies_txt_in, fragments_sdf_out, fragments_txt_out, frequencies_txt_out, frag_frequencies_txt_out, limit=None, filter=False, first=None, inchi_filter=None, pains=False, sort=False, test=False):

    print('Loading %s' %first)

    if first is not None:
        fragment_database_mol = get_fragment_database('%s_%s.sdf' %(fragments_sdf_in, first))

        frequencies = get_bond_frequencies('%s_%s.txt' %(frequencies_txt_in, first) )

        frag_frequencies = get_frag_frequencies('%s_%s.txt' %(frag_frequencies_txt_in, first))

        frag_mapping = get_frag_mapping('%s_%s.txt' %(fragments_txt_in, first) )

    else:
        fragment_database_mol = get_fragment_database('%s.sdf' %(fragments_sdf_in) )

        frequencies = get_bond_frequencies('%s.txt' %(frequencies_txt_in) )

        frag_frequencies = get_frag_frequencies('%s.txt' %(frag_frequencies_txt_in))

        frag_mapping = get_frag_mapping('%s.txt' %(fragments_txt_in) )


    fragment_database = []

    for i in fragment_database_mol:
        fragment_database.append(i.graph)

    frequencies = update_bond_frequencies(frequencies, frag_mapping)

    if sort:
        fragment_database, frequencies, frag_frequencies, frag_mapping = sort_fragments(fragment_database, frequencies, frag_frequencies, frag_mapping)

        save_frequencies_txt(frequencies, frequencies_txt_out)

        save_fragments_sdf(fragment_database, fragments_sdf_out)

        save_frag_frequencies_txt(frag_frequencies, frag_frequencies_txt_out)

        save_fragments_txt(fragment_database, fragments_txt_out)

        sys.exit('Fragment data sorted')


    if filter:
        filter_list = filter_database(fragment_database_mol, inchi_filter, pains)

        with open('filter_in.sdf', 'w') as outfile:

            for i in filter_list:

                mol = fragment_database_mol[i]

                lines = molecule_to_sdf(mol)

                for line in lines:
                    outfile.write(line)

                outfile.write('$$$$\n')               

        frequencies = remove_bond_frequencies(frequencies, filter_list)

        save_frequencies_txt(frequencies, 'filter_bond_frequencies.txt')       

        sys.exit('Fragment database filetered')

    if limit is not None:
        print('Before limit ', len(fragment_database))
        fragment_database, frequencies, frag_frequencies, frag_mapping = update_limit(limit, fragment_database, frequencies, frag_frequencies, frag_mapping)
        print('After limit ', len(fragment_database))

    if first is None:
        raise Exception ('Need to define first argument')

    for i in range(first+1, first+n):

        print('Loading %s' %i)

        fragments_sdf_2 = '%s_%s.sdf' %(fragments_sdf_in, i)
        fragments_txt_2 = '%s_%s.txt' %(fragments_txt_in, i) 
        frequencies_txt_2 = '%s_%s.txt' %(frequencies_txt_in, i) 
        frag_frequencies_txt_2 = '%s_%s.txt' %(frag_frequencies_txt_in, i)

        if test:
            frag_frequencies2 = get_frag_frequencies(frag_frequencies_txt_2)
            print('Before limit', len(frag_frequencies2))

            limit_count = 0

            for i in frag_frequencies2:
                if i >= limit:
                    limit_count += 1

            print('After limit', limit_count)

            continue

        print(fragments_sdf_2, fragments_txt_2, frequencies_txt_2, frag_frequencies_txt_2, fragments_sdf_out, fragments_txt_out, frequencies_txt_out, frag_frequencies_txt_out)

        combine_fragment_databases(fragment_database, frequencies, frag_frequencies, frag_mapping, fragments_sdf_2, fragments_txt_2, frequencies_txt_2, frag_frequencies_txt_2, fragments_sdf_out, fragments_txt_out, frequencies_txt_out, frag_frequencies_txt_out, limit)

def update_limit(limit, fragment_database, bond_frequencies, frag_frequencies, frag_mapping):

    # set mapping list and loop through elements of frag_frequencies, if element < limit then set mapping to -1
    # create new fragment database with fragments that are within limit
    mapping = []
    new_fragment_database = []
    new_frag_frequencies = []
    new_frag_mapping = []
    check = {}
    j = 0
    for i in range(len(frag_frequencies)):
        if frag_frequencies[i] < limit:
            mapping.append(-1)
        else:
            mapping.append(j)
            new_fragment_database.append(fragment_database[i])
            new_frag_frequencies.append(frag_frequencies[i])
            new_frag_mapping.append(frag_mapping[i])
            j += 1

    # update bond_frequencies: remove frequencies for removed fragments and update fragment numbers according to mapping

    new_bond_frequencies = {}

    for key, val in bond_frequencies.items():
        i = key[0]
        j = key[1]
        k = key[2]
        l = key[3]
        if mapping[i] != -1 and mapping[j] != -1:
            check[key] = val
            new_bond_frequencies[mapping[i], mapping[j], k, l] = val

    return new_fragment_database, new_bond_frequencies, new_frag_frequencies, new_frag_mapping

def sort_fragments(fragment_database, bond_frequencies, frag_frequencies, frag_mapping):

    frag_frequencies_np = np.array(frag_frequencies)
    sort_index = list(np.argsort(-1*frag_frequencies_np))

    mapping = {}

    for i in range(len(sort_index)):
        mapping[sort_index[i]] = i

    new_fragment_database = []
    new_bond_frequencies = {}
    new_frag_frequencies = []
    new_frag_mapping = []

    # add fragments to new_fragment_database according to sorted order, same for new_frag_mapping
    for i in range(len(sort_index)):
        new_fragment_database.append(fragment_database[sort_index[i]])
        new_frag_mapping.append(frag_mapping[sort_index[i]])

    new_frag_frequencies = sorted(frag_frequencies, reverse=True)

    # update bond frequencies for new fragment indeces, larger index should be j
    for key, val in bond_frequencies.items():

        i = mapping[key[0]]
        j = mapping[key[1]]
        k = key[2]
        l = key[3]

        if i <= j:
            new_bond_frequencies[i,j,k,l] = val            
        else:
            new_bond_frequencies[j,i,l,k] = val

    return new_fragment_database, new_bond_frequencies, new_frag_frequencies, new_frag_mapping

def renumber_frequencies(fragments_txt_in, frequencies_txt_in, frequencies_txt_out):

    frequencies = get_bond_frequencies(frequencies_txt_in)

    frag_mapping = get_frag_mapping(fragments_txt_in)

    frequencies = update_bond_frequencies(frequencies, frag_mapping)

    save_frequencies_txt(frequencies, frequencies_txt_out)

if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Combine fragmented molecules')
    parser.add_argument('-n','--n_files', help='Number of fragment files to combine',required=True, type=int)
    parser.add_argument('-i','--in_sub', help='Input subscript',required=True)
    parser.add_argument('-o','--out_sub', help='Output subscript',required=True)
    parser.add_argument('-f','--first', help='First file index to consider',required=False, type=int)
    parser.add_argument('-l','--limit', help='Limit for minimum fragment frequency to consider',required=False, type=int)
    parser.add_argument('--filter', action='store_true', help='Filter fragment database', required=False)
    parser.add_argument('--inchi_filter', help='Inchi list to filter',required=False)
    parser.add_argument('--pains', action='store_true', help='Filter fragments with pains', required=False)
    parser.add_argument('--sort', action='store_true', help='Sort fragment data and exit', required=False)
    parser.add_argument('--test', action='store_true', help='Test run', required=False)

    args = parser.parse_args()

    n = args.n_files

    in_sub = args.in_sub
    out_sub = args.out_sub

    if args.first is not None:
        first = args.first * n

    fragments_sdf_in = 'fragments%s' %in_sub
    fragments_txt_in = 'fragments%s' %in_sub
    frequencies_txt_in = 'frequencies%s' %in_sub
    frag_frequencies_txt_in = 'frag_frequencies%s' %in_sub
    fragments_sdf_out = 'fragments%s.sdf' %out_sub
    fragments_txt_out = 'fragments%s.txt' %out_sub
    frequencies_txt_out = 'frequencies%s.txt' %out_sub
    frag_frequencies_txt_out = 'frag_frequencies%s.txt' %out_sub

    loop(n, fragments_sdf_in, fragments_txt_in, frequencies_txt_in, frag_frequencies_txt_in, fragments_sdf_out, fragments_txt_out, frequencies_txt_out, frag_frequencies_txt_out, filter=args.filter, first=args.first, limit=args.limit, inchi_filter=args.inchi_filter, pains=args.pains, sort=args.sort, test=args.test)

    print('Normal termination')