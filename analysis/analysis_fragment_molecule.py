import argparse

import pandas as pd

from networkx.algorithms import isomorphism
from rdkit import Chem

from pymolgen.canonicalise_tautomer import canonicalise_tautomer
from pymolgen.molecule_formats import molecule_from_smiles
from pymolgen.fragment_mol import get_fragments_dataset, get_canonical_mapping, get_atom_list, compound_dict, node_compare_element
from pymolgen.fragment_molecule import *
from pymolgen.fragment_molecule_builder import *


def get_fragment_index(fragment, fragment_database, fragment_database_len=None, atom_list_all=None):
    """
    Get fragment index by searching fragment database

    Parameters
    ----------
    fragment : Molecule object
    fragment_database : list of Molecule objects

    Returns
    -------
    index[0] : int
        Index for fragment in fragment database
    map : dict 
        Dictionary mapping the original atom numbers in fragment to those in the fragment databse
    """

    index = []

    map = get_canonical_mapping(fragment)

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
            atom_list_all_i = get_atom_list(fragment_database[i])

        if fragment_len == fragment_database_len_i and fragment_atom_list == atom_list_all_i:

            gm = isomorphism.GraphMatcher(fragment, fragment_database[i], node_match=node_compare_element)

            if gm.is_isomorphic():
                index.append(i)
                newmap = gm.mapping

                map = get_canonical_mapping(fragment_database[i])

                map = compound_dict(newmap, map)

    if len(index) > 1:
        print(index)
        raise Exception('fragment in fragment_database more than once')

    return index[0], map


def analyse_molecule(inchi, fragment_database):

    rdmol = Chem.MolFromInchi(inchi)

    smi = Chem.MolToSmiles(rdmol)

    smi = canonicalise_tautomer(smi)

    mol = molecule_from_smiles(smi)

    fragments, pairs, bonds = get_fragments_dataset(mol, carbonyl=True, fluorine=True)

    for fragment in fragments:

        get_fragment_index(fragment, fragment_database)        




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

            analyse_molecule(inchi, fragment_database)









