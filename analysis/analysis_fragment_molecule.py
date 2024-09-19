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


# Define the custom print function
def print_with_line(*args, **kwargs):
    # Get the current frame and extract the line number from the caller
    frame = inspect.currentframe()
    caller_frame = frame.f_back  # Get the frame of the caller
    line_number = caller_frame.f_lineno  # Extract the line number
    
    # Call the original print function with the line number prepended
    builtins.print(f"[Line {line_number}]", *args, **kwargs)


# Override the built-in print function using partial
print = partial(print_with_line)


def analyse_molecule(inchi, fragment_database, fragment_database_graph):
    """
    Analyse a molecule in inchi format in terms of its fragment molecule structure (as a graph).
    If a parent_mol is given, then the constituent fragments for the parent strucrue present in the final molecule will be converted into a single parent fragment.
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

    calculate_build_probability(bond_frequencies, fragment_database_graph, fragment_molecule, root)

    return str(f)


def calculate_build_probability(bond_frequencies, fragment_database_graph, fragment_molecule, root):

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


def get_ordered_bonds(fragment_molecule, root):

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

    #map = get_canonical_mapping(fragment)

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
    parser.add_argument('-rf', '--read_fragment_database', help='Read fragment database from file containing attachment points and canonical mapping', required=False)

    args = parser.parse_args()

    if args.input == args.output:
        raise Exception('Same input and output')

    fragment_database = get_fragment_database(args.fragments_sdf)

    if args.read_fragment_database is not None:
        fragment_database_graph = read_fragment_database_graph(args.read_fragment_database)
    else:    
        fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    parent_mol = molecule_from_sdf(args.parent_file)

    attachment_points = []

    # remove hydrogens from parent and determine atoms that will have open valence
    for i in args.remove_hydrogens:
        parent_mol = parent_mol.remove_atom(i)
        for j in parent_mol.free_valence_list:
            if j not in attachment_points:
                attachment_points.append(j)

    # include parent in fragment_database and fragment_database_graph
    parent_id = len(fragment_database)
    fragment_database.add_mol(parent_mol)
    fragment_database_graph.add_fragment(parent_id, attachment_points)
    fragment_database_graph.fragments[parent_id].set_attribute('frag_id', parent_id)
    fragment_database_graph.fragments[parent_id].set_canonical_mapping(fragment_database)

    with open(args.input) as infile, open(args.output, 'w') as outfile:

        for line in infile:

            inchi = line.strip().split()[0]

            string_representation = analyse_molecule(inchi, fragment_database, fragment_database_graph, parent_mol)

            outfile.write(f'{string_representation}\n')

            f = generate_fragment_molecule_from_string(string_representation, fragment_database_graph)

            mol = convert_fragment_molecule_to_mol(f, fragment_database)

            assert inchi == molecule_to_inchi(mol)











