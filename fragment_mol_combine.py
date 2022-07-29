import sys,os

from pymolgen.fragment_mol import *
from pymolgen.fragment_builder import *

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

def loop(n, first, fragments_sdf_in, fragments_txt_in, frequencies_txt_in, frag_frequencies_txt_in, fragments_sdf_out, fragments_txt_out, frequencies_txt_out, frag_frequencies_txt_out, limit=None, test=None):

    print('Loading %s' %first)

    fragment_database_mol = get_fragment_database('%s%s.sdf' %(fragments_sdf_in, first))

    fragment_database = []

    for i in fragment_database_mol:
        fragment_database.append(i.graph)

    frequencies = get_bond_frequencies('%s%s.txt' %(frequencies_txt_in, first) )

    frag_frequencies = get_frag_frequencies('%s%s.txt' %(frag_frequencies_txt_in, first))

    frag_mapping = get_frag_mapping('%s%s.txt' %(fragments_txt_in, first) )

    frequencies = update_bond_frequencies(frequencies, frag_mapping)

    if limit is not None:
        print('Before limit ', len(fragment_database))
        fragment_database, frequencies, frag_frequencies, frag_mapping = update_limit(limit, fragment_database, frequencies, frag_frequencies, frag_mapping)
        print('After limit ', len(fragment_database))

    for i in range(first+1, first+n):

        print('Loading %s' %i)

        fragments_sdf_2 = '%s%s.sdf' %(fragments_sdf_in, i)
        fragments_txt_2 = '%s%s.txt' %(fragments_txt_in, i) 
        frequencies_txt_2 = '%s%s.txt' %(frequencies_txt_in, i) 
        frag_frequencies_txt_2 = '%s%s.txt' %(frag_frequencies_txt_in, i)

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

    # update bond_frequencies

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
    parser.add_argument('-f','--first', help='First file index to consider',required=True, type=int)
    parser.add_argument('-l','--limit', help='Limit for minimum fragment frequency to consider',required=False, type=int)
    parser.add_argument('--test', action='store_true', help='Test run', required=False)

    args = parser.parse_args()

    n = args.n_files

    in_sub = args.in_sub
    out_sub = args.out_sub
    first = args.first * n

    fragments_sdf_in = 'fragments%s_' %in_sub
    fragments_txt_in = 'fragments%s_' %in_sub
    frequencies_txt_in = 'frequencies%s_' %in_sub
    frag_frequencies_txt_in = 'frag_frequencies%s_' %in_sub
    fragments_sdf_out = 'fragments%s.sdf' %out_sub
    fragments_txt_out = 'fragments%s.txt' %out_sub
    frequencies_txt_out = 'frequencies%s.txt' %out_sub
    frag_frequencies_txt_out = 'frag_frequencies%s.txt' %out_sub

    loop(n, first, fragments_sdf_in, fragments_txt_in, frequencies_txt_in, frag_frequencies_txt_in, fragments_sdf_out, fragments_txt_out, frequencies_txt_out, frag_frequencies_txt_out, limit=args.limit, test=args.test)