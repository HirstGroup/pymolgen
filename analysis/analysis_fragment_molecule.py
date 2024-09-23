import argparse
import builtins
import inspect
import networkx

import pandas as pd

from functools import partial
from networkx.algorithms import isomorphism
from rdkit import Chem

from pymolgen.canonicalise_tautomer import canonicalise_tautomer
from pymolgen.molecule_formats import molecule_from_smiles
from pymolgen.fragment_mol import get_fragments_dataset, get_canonical_mapping, get_atom_list, compound_dict, node_compare_element
from pymolgen.fragment_molecule import *
from pymolgen.fragment_molecule_builder import *


def print_with_line(*args, **kwargs):
    """
    Define the custom print function
    """
    # Get the current frame and extract the line number from the caller
    frame = inspect.currentframe()
    caller_frame = frame.f_back  # Get the frame of the caller
    line_number = caller_frame.f_lineno  # Extract the line number
    
    # Call the original print function with the line number prepended
    builtins.print(f"[Line {line_number}]", *args, **kwargs)


# Override the built-in print function using partial
print = partial(print_with_line)


def analyse_molecule(inchi, fragment_database, fragment_database_graph, bond_frequencies=None, root=None, version=None):
    """
    Analyse a molecule in inchi format in terms of its fragment molecule structure (as a graph).

    Parameters
    ----------
    inchi : str
        Inchi string to analyse
    fragment_database : list of Molecule objects
        Fragment database in Molecule format
    fragment_database_graph : FragmentMolecule object
        Fragment database in FragmentMolecule format
    bond_frequencies : dict of (i,k) into dict of (j,l):frequencies
        Bond frequencies in dictionary format    
    root : int, optional
        Index number of root fragment in database.
        Root fragment means the fragment from which the rest of the molecule is built.
        (Fragments are obtained according to fragmentation rules, so could be different to parent fragment)
    version : int, optional
        Version for calculation of build probability factor

    Returns
    -------
    """

    rdmol = Chem.MolFromInchi(inchi)

    smi = Chem.MolToSmiles(rdmol)

    smi = canonicalise_tautomer(smi)

    mol = molecule_from_smiles(smi)

    f = FragmentMolecule()

    fragments, pairs, bonds = get_fragments_dataset(mol, carbonyl=True, fluorine=True)

    index_list = []
    mapping_list = []

    for fragment in fragments:

        index, mapping = get_fragment_index(fragment, fragment_database)

        index_list.append(index)
        mapping_list.append(mapping)

        f.add_fragment(index, fragment_database[index].free_valence_list)

    assert len(pairs) == len(bonds)

    for pair, bond in zip(pairs, bonds):

        i = pair[0]
        j = pair[1]
        k = mapping_list[i][bond[0]]
        l = mapping_list[j][bond[1]]

        f.add_bond(i, j, k, l)

    if root is not None:
        f._graph._build_probability = calculate_build_probability_version2(bond_frequencies, fragment_database_graph, f, root, version)

    return str(f)


def calculate_build_probability(bond_frequencies, fragment_database_graph, fragment_molecule, root):
    """
    Calculate build probability for molecule by traversing breath first search through nodes

    Parameters
    ----------
    bond_frequencies : dict of (i,k) into dict of (j,l):frequencies
        Bond frequencies in dictionary format
    fragment_database_graph : FragmentMolecule object
        Fragment database in FragmentMolecule format
    fragment_molecule : FragmentMolecule object
        Molecule in FragmentMolecule format
    root : int
        Index of molecule in fragment databse

    Returns
    -------
    build_probability : float
        Build probability of molecule
    """

    build_probability = 1.0

    root_index = fragment_molecule.list_frag_id().index(root)

    ordered_bonds = get_ordered_bonds(fragment_molecule, root_index)

    bonds = fragment_molecule._graph.bonds

    bond_dict = make_bond_dict(bonds)

    for bond in ordered_bonds:

        i, j = bond

        k, l = bond_dict[bond]

        i_id = fragment_molecule.list_frag_id()[i]
        j_id = fragment_molecule.list_frag_id()[j]

        k_can = fragment_database_graph.fragments[i_id].get_canonical_mapping()[k]
        l_can = fragment_database_graph.fragments[j_id].get_canonical_mapping()[l]

        freq = bond_frequencies[(i_id,k_can)][(j_id,l_can)]

        fragment_bonds = bond_frequencies[(i_id, k_can)]

        total_freq = sum(fragment_bonds.values())

        bond_freq = freq/total_freq

        build_probability *= bond_freq

    print(fragment_molecule._graph._fragments)

    list_free_valence_points = fragment_molecule.list_free_valence_points()

    # calculate multiplication factor for build_probability to take number of attachment points into account

    factor = 1

    for fragment_index, fragment in fragment_molecule._graph._fragments.items():
        print(fragment_index, fragment.attachment_points)

        n_attachment_points = len(fragment.attachment_points)
        n_attachment_points_free = len(list_free_valence_points[fragment_index])
        n_neighbours = n_attachment_points - n_attachment_points_free

        print(fragment_index, n_neighbours)

        if fragment_index == root_index:
            n_max = n_attachment_points
        else:
            n_max = n_attachment_points - 1

        n_min = n_attachment_points_free

        for i in range(n_max, n_min, -1):
            print(i)
            factor *= 1/i

    print('factor', factor)

    print(fragment_molecule.list_free_valence_points())

    build_probability *= factor

    return build_probability


def traverse_least_neighbors(fragment_molecule, root):
    """
    Get bonds by custom traversal from the root, choosing the neighbor with the fewest further neighbors

    Parameters
    ----------
    fragment_molecule : FragmentMolecule object
        Fragment molecule in FragmentMolecule format

    Returns
    -------
    bonds : list of (i,j,k,l)
        List of bonds as tuples
    """

    graph = fragment_molecule._graph.convert_to_networkx()

    visited = set()  # To keep track of visited nodes
    to_visit = [(root, None)]  # Store tuples of (current node, parent node)
    
    bonds = []

    while to_visit:
        # Pop the current node (DFS-style: process one node before continuing others)
        node, parent = to_visit.pop()
        
        if node not in visited:
            if parent is not None:
                print(f"Visited Node: {node} (bonded to {parent})")
                bonds.append((parent, node))
            else:
                print(f"Visited Node: {node} (root)")

            visited.add(node)
            
            # Get unvisited neighbors of the current node
            unvisited_neighbors = [n for n in graph.neighbors(node) if n not in visited]
            
            # Sort neighbors by the number of their further neighbors (degree)
            sorted_neighbors = sorted(unvisited_neighbors, key=lambda n: graph.degree(n))
            
            # Add sorted neighbors to the to_visit list in reverse order,
            # passing the current node as the parent
            to_visit.extend([(n, node) for n in sorted_neighbors[::-1]])

    return bonds


def calculate_build_probability_version2(bond_frequencies, fragment_database_graph, fragment_molecule, root, version=1):
    """
    2nd version of calculate build probability function. It traverses nodes so that the total number of available attachment
    points is minimised during the fragment addition.
    
    It takes version as a parameter. 
    If version=1, then the multiplication factor takes the total number of attachment points into account. This version 
    gives the minimum expected build probability for a molecule built using version=1 (i.e. the probability can be larger
    since build probabilities are added up for all different ways to achieve the same molecule).

    If version=2, then the multiplication factor will take the number of available attachment points of a node times the number
    of fragments in the molecule. In this version, the build probability will be the same regardless of the order of fragment
    addition.

    Parameters
    ----------
    bond_frequencies : dict of (i,k) into dict of (j,l):frequencies
        Bond frequencies in dictionary format
    fragment_database_graph : FragmentMolecule object
        Fragment database in FragmentMolecule format
    fragment_molecule : FragmentMolecule object
        Molecule in FragmentMolecule format
    root : int
        Index of molecule in fragment databse
    version : int, optional
        Version for calculation of build probability factor

    Returns
    -------
    build_probability : float
        Build probability of molecule
    """
    
    build_probability = 1.0

    print(fragment_molecule)
    print(fragment_molecule.list_frag_id())
    root_index = fragment_molecule.list_frag_id().index(root)

    ordered_bonds = traverse_least_neighbors(fragment_molecule, root_index)

    bonds = fragment_molecule._graph.bonds

    bond_dict = make_bond_dict(bonds)

    original_attachment_points = []

    for fragment_id, fragment in fragment_molecule._graph.fragments.items():
        print(fragment.attachment_points)
        original_attachment_points.append(len(fragment.attachment_points))

    print(original_attachment_points)

    current_attachment_points = [0 for i in original_attachment_points]

    current_attachment_points[root_index] = original_attachment_points[root_index]

    #n_attachment_points = original_attachment_points[root_index]

    n_fragments = 1

    for bond in ordered_bonds:

        i, j = bond
        k, l = bond_dict[bond]
        print(i,j,k,l)

        i_frag_id = fragment_molecule.list_frag_id()[i]
        j_frag_id = fragment_molecule.list_frag_id()[j]

        print(i_frag_id, j_frag_id)

        k_can = fragment_database_graph.fragments[i_frag_id].get_canonical_mapping()[k]
        l_can = fragment_database_graph.fragments[j_frag_id].get_canonical_mapping()[l]

        print(k_can, l_can)

        bond_freq_all = bond_frequencies[(i_frag_id, k_can)]
        print(bond_freq_all)
        bond_freq = bond_freq_all[(j_frag_id, l_can)]
        print(bond_freq)

        if version == 1:
            factor = sum(current_attachment_points)
        elif version == 2:
            factor = current_attachment_points[i] * n_fragments
        else:
            raise Exception(f'Version {version} not implemented')

        current_attachment_points[i] += -1
        current_attachment_points[j] += original_attachment_points[j] - 1

        n_fragments += 1

        fragment_bonds = bond_frequencies[(i_frag_id, k_can)]

        total_freq = sum(fragment_bonds.values())        

        attachment_probability = bond_freq / (total_freq * factor)

        build_probability *= attachment_probability

    print(build_probability)

    return build_probability


def get_ordered_bonds(fragment_molecule, root):
    """
    Get bonds from molecule by traversing in breath first search order.

    Parameters
    ----------
    fragment_molecule : FragmentMolecule object
        Fragment molecule in FragmentMolecule format

    Returns
    -------
    bonds : list of (i,j,k,l)
        List of bonds as tuples
    """

    networkx_graph = fragment_molecule._graph.convert_to_networkx()

    visited = set()

    bfs_nodes = list(networkx.bfs_tree(networkx_graph, root))

    bonds = []

    for node in bfs_nodes:

        visited.add(node)

        for neighbor in networkx_graph.neighbors(node):
            if neighbor not in visited:
                bonds.append((node, neighbor))

    return bonds


def make_bond_dict(bonds):
    """
    Return bond dictionary for bonds

    Parameters
    ----------
    bonds : list of (i,j,k,l)

    Returns
    -------
    bond_dict : dict of (i,k) -> (k,l)
    """

    bond_dict = {}

    for bond in bonds:

        i,j,k,l = bond

        assert (i,j) not in bond_dict

        bond_dict[(i,j)] = (k,l)

    return bond_dict


def convert_parent(fragment_molecule, parent_mol):

    g1 = fragment_molecule.convert_to_networkx()
    g2 = parent_mol.convert_to_networkx()

    gm = isomorphism.GraphMatcher(g1, g2, node_match=lambda n1,n2:n1['frag_id']==n2['frag_id'], edge_match= lambda e1,e2: e1['atoms'] == e2['atoms'])

    if gm.subgraph_is_isomorphic():
        matching = gm.mapping
        print(matching)

    else:
        raise Exception('Parent fragment not present in molecule')


def get_fragment_index(fragment, fragment_database, fragment_database_len=None, atom_list_all=None):
    """
    Get fragment index by searching fragment database

    Parameters
    ----------
    fragment : Networkx object
    fragment_database : list of Molecule objects

    Returns
    -------
    index[0] : int
        Index for fragment in fragment database
    map : dict 
        Dictionary mapping the original atom numbers in fragment to those in the fragment databse
    """

    index = []

    fragment_len = len(fragment)

    fragment_atom_list = get_atom_list(fragment)

    for i in range(len(fragment_database)):

        if fragment_database_len is not None:
            fragment_database_len_i = fragment_database_len[i]
        else:
            fragment_database_len_i = len(fragment_database[i])

        if atom_list_all is not None:
            atom_list_all_i = atom_list_all[i]
        else:
            atom_list_all_i = get_atom_list(fragment_database[i].graph)

        if fragment_len == fragment_database_len_i and fragment_atom_list == atom_list_all_i:

            gm = isomorphism.GraphMatcher(fragment, fragment_database[i].graph, node_match=node_compare_element)

            if gm.is_isomorphic():
                index.append(i)
                newmap = gm.mapping

                #map = get_canonical_mapping(fragment_database[i].graph)

                #map = compound_dict(newmap, map)

    if len(index) > 1:
        print(index)
        raise Exception('fragment in fragment_database more than once')

    return index[0], newmap


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Analyse molecules in terms of FragmentMolecule class')

    # required arguments
    parser.add_argument('-a','--fragments_sdf', help='SDF file of fragments', required=True)
    parser.add_argument('-i','--input', help='Input inchi file', required=True)
    parser.add_argument('-o','--output', help='Output file name in FragmentMolecule format', required=True)
    parser.add_argument('-p','--parent_file', help='Parent Structure File in SDF format', required=True)
    parser.add_argument('-r','--remove_hydrogens', type=int, nargs='+', help='Space-separated hydrogen atoms that will be created as attachment points, numbered from 0', required=True)

    # optional arguments
    parser.add_argument('-d','--frequencies_txt', help='Bond frequencies dictionary in txt file', required=False)
    parser.add_argument('-rd', '--read_bond_frequencies_dict', help='Read bond frequencies dict from file', required=False)
    parser.add_argument('-rf', '--read_fragment_database', help='Read fragment database from file containing attachment points and canonical mapping', required=False)
    parser.add_argument('--root', type=int, help='Index for root fragment, this variable also triggers build probability calculation', required=False)
    parser.add_argument('--version', default=1, type=int, help='Version for build probability factor, version 1 gives different build probabilities according to the order of fragment addition, version 2 gives same build probabilities for any order', required=False)

    args = parser.parse_args()

    if args.input == args.output:
        raise Exception('Same input and output')

    fragment_database = get_fragment_database(args.fragments_sdf)

    if args.read_fragment_database is not None:
        fragment_database_graph = read_fragment_database_graph(args.read_fragment_database)
    else:    
        fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    if args.read_bond_frequencies_dict is None:
        bond_frequencies = get_bond_frequencies(args.frequencies_txt)
        bond_frequencies = bond_frequencies_to_np(bond_frequencies)
        bond_frequencies = convert_bond_freq_np_to_dict(fragment_database_graph, bond_frequencies)
    else:
        bond_frequencies = read_bond_frequencies_dict(args.read_bond_frequencies_dict)

    if args.root is None:
        bond_frequencies = None

    parent_mol = molecule_from_sdf(args.parent_file)

    parent, bond_frequencies, fragment_database, fragment_database_graph = prepare_parent(bond_frequencies, fragment_database, fragment_database_graph, args.parent_file, args.parent_fragment_file_list, args.parent_mapping_1, args.remove_hydrogens, args.remove_hydrogens_parent_fragment)

    with open(args.input) as infile, open(args.output, 'w') as outfile:

        for line in infile:

            inchi = line.strip().split()[0]

            string_representation = analyse_molecule(inchi, fragment_database, fragment_database_graph, bond_frequencies, args.root, args.version)

            outfile.write(f'{string_representation}\n')

            f = generate_fragment_molecule_from_string(string_representation, fragment_database_graph)

            mol = convert_fragment_molecule_to_mol(f, fragment_database)

            assert inchi == molecule_to_inchi(mol)











