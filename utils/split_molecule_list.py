import argparse
import numpy as np

def split_molecule_list(molecule_list, n):

    """
    Split a molecule list in lists of equal size so that each list adds up to a similar value of build probabilities
    """
    size = len(molecule_list)

    build_probability_list = np.zeros(size)

    for i in range(size):
        build_probability_list[i] = molecule_list[i].split(':')[2]

    sort_index = np.argsort(build_probability_list)

    output_mol_list = []
    for i in range(n):
        output_mol_list.append([])

    remainder = size % 2

    sort_index_1 = sort_index[:size // 2 + remainder]
    sort_index_2 = sort_index[size // 2 + remainder:][::-1]

    for i in range(len(sort_index_1)):

        list_index = i % n

        output_mol_list[list_index].append(molecule_list[sort_index_1[i]])

    for i in range(len(sort_index_2)):

        list_index = i % n

        output_mol_list[list_index].append(molecule_list[sort_index_2[i]])

    return output_mol_list


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Split a molecule list file in FragmentMolecule format in files of equal total build probability')

    # required arguments
    parser.add_argument('-i','--input', help='Input file with molecules in FragmentMolecule format',required=True)
    parser.add_argument('-n','--n_files', type=int, help='Number of output files to split into',required=True)
    parser.add_argument('-o','--output', help='Base name of output files, without txt extension',required=True)

    args = parser.parse_args()

    fragment_molecule_list = []

    with open(args.input) as infile:

        for line in infile:

            fragment_molecule_list.append(line.strip())

    fragment_molecule_list_list = split_molecule_list(fragment_molecule_list, args.n_files)

    for n, fragment_molecule_list in enumerate(fragment_molecule_list_list):

        with open(f'{args.output}_{n}.txt', 'w') as outfile:

            for fragment_molecule in fragment_molecule_list:

                outfile.write(f'{fragment_molecule}\n')

