import argparse
import heapq
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


def divide_into_n_lists(strings, n):
    # Parse the strings and extract the C values
    items = [(s, float(s.split(':')[-1])) for s in strings]
    
    # Sort items based on the C value in descending order
    items.sort(key=lambda x: x[1], reverse=True)
    
    # Create a list of n empty lists and initialize their sums
    partitions = [[] for _ in range(n)]
    sums = [0] * n  # To keep track of the sums of the partitions

    # Use a min heap to keep track of which partition has the smallest sum
    heap = [(0, i) for i in range(n)]  # Each element is a tuple (sum, partition_index)
    heapq.heapify(heap)

    # Distribute the strings based on their C value
    for string, c_value in items:
        smallest_sum, partition_index = heapq.heappop(heap)  # Get the partition with the smallest sum
        partitions[partition_index].append(string)           # Add the string to this partition
        new_sum = smallest_sum + c_value                     # Update the sum of this partition
        heapq.heappush(heap, (new_sum, partition_index))     # Push the updated partition back into the heap

    return partitions


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

    fragment_molecule_list_list = divide_into_n_lists(fragment_molecule_list, args.n_files)

    for n, fragment_molecule_list in enumerate(fragment_molecule_list_list):

        with open(f'{args.output}_{n}.txt', 'w') as outfile:

            for fragment_molecule in fragment_molecule_list:

                outfile.write(f'{fragment_molecule}\n')

