#!/usr/bin/env python

import argparse
import copy
import numpy as np
import os
import sys

from pymolgen.fragment_molecule import *
from pymolgen.generate import SDFDatasetLargeRAM
from pymolgen.molecule_formats import *
from pymolgen.fragment_builder import bond_frequencies_to_np, get_bond_frequencies, get_fragment_database, get_fragment_bond_frequencies_np

from functools import partial
print = partial(print, flush=True)


def extend_molecule(fragment_id, bond_frequencies, fragment_database):

    output_mol_list = []

    mol = fragment_database[fragment_id]

    free_valence_list = mol.free_valence_list
    for atom in free_valence_list:

        fragment_bonds, fragment_bond_frequencies = get_fragment_bond_frequencies_np(fragment_id, atom, bond_frequencies)

        for bond in fragment_bonds:
            i = bond[0]
            j = bond[1]
            k = bond[2]
            l = bond[3]

            f = FragmentMolecule()
            f.add_fragment(fragment_id, mol.attach_points)

            if i == fragment_id:

                f.add_fragment(j, fragment_database[j].attach_points)
                f.add_bond(0, 1, k, l)

            elif j == fragment_id:

                f.add_fragment(j, fragment_database[i].attach_points)
                f.add_bond(0, 1, l, k)

            else:
                sys.error('fragmend_id not in bond', bond)


            output_mol_list.append(f)

    return output_mol_list


def extend_molecule_list(FragmentMolecule_list, bond_frequencies, fragment_database_graph, depth=False, threshold=None):

    output_mol_list = []

    for f in FragmentMolecule_list:
        free_valence_list = f.list_free_valence_points()
        total_free_valence = f.get_total_free_valence()

        for x in range(len(free_valence_list)):

            fragment_id = f.get_frag_id(x)

            for atom in free_valence_list[x]:

                atom_can = fragment_database_graph.fragments[fragment_id].get_canonical_mapping()[atom]

                fragment_bonds, fragment_bond_frequencies = get_fragment_bond_frequencies_np(fragment_id, atom_can, bond_frequencies)

                total_freq = np.sum(fragment_bond_frequencies)

                for bond, bond_freq in zip(fragment_bonds, fragment_bond_frequencies):

                    attachment_probability = bond_freq / ( total_freq * total_free_valence)

                    new_build_probability = attachment_probability * f.get_build_probability()

                    if threshold is not None and new_build_probability < threshold:
                        # do not build molecule if its build probability is below the threshold
                        continue

                    i = bond[0]
                    j = bond[1]
                    k = bond[2]
                    l = bond[3]

                    f2 = copy.deepcopy(f)

                    # if i corresponds to left fragment j is right fragment
                    if i == fragment_id and k == atom_can:
                        node_id = f2.add_fragment(j, fragment_database_graph.fragments[j].attachment_points)
                        f2.add_bond(x, node_id, atom, l, attachment_probability)

                    # if j corresponds to left fragment i is right framgent
                    elif j == fragment_id and l == atom_can:
                        node_id = f2.add_fragment(i, fragment_database_graph.fragments[i].attachment_points)
                        f2.add_bond(x, node_id, atom, k, attachment_probability)

                    else:
                        sys.error('fragmend_id and atom_can not in bond', bond, atom_can)

                    if depth is not None:
                        total = len(output_mol_list)

                        if total % 10000 == 0:
                            print(f'DEPTH {depth} TOTAL {total}')

                    output_mol_list.append(f2)

    return output_mol_list

def extend_molecule_list_count(FragmentMolecule_list, bond_frequencies, fragment_database_graph, depth=None):

    total = 0

    for f in FragmentMolecule_list:

        free_valence_list = f.list_free_valence_points()

        for x in range(len(free_valence_list)):

            fragment_id = f.get_frag_id(x)

            for atom in free_valence_list[x]:

                atom_can = fragment_database_graph.fragments[fragment_id].get_canonical_mapping()[atom]

                fragment_bonds, fragment_bond_frequencies = get_fragment_bond_frequencies_np(fragment_id, atom_can, bond_frequencies)

                total += len(fragment_bonds)

                if total % 10000 == 0:
                    print(f'DEPTH {depth} TOTAL {total}')

    return total


def extend_molecule_list_depth(FragmentMolecule_list, bond_frequencies, fragment_database_graph, depth, output=None, unique=True, threshold=None):

    for i in range(depth):

        FragmentMolecule_list = extend_molecule_list(FragmentMolecule_list, bond_frequencies, fragment_database_graph, i + 1, threshold)

        print(f'FINAL DEPTH {i+1} TOTAL {len(FragmentMolecule_list)}')

        print('unique = ', unique)

        if unique is True:

            FragmentMolecule_list = get_unique_molecule_list(FragmentMolecule_list, fragment_database=fragment_database_graph)

            print(f'FINAL DEPTH {i+1} TOTAL UNIQUE {len(FragmentMolecule_list)}')

        if len(FragmentMolecule_list) == 0:
            # stop building molecules
            return FragmentMolecule_list

        if output is not None:
            with open('%s-depth%s.inchi' %(output, i+1), 'w') as f:
                for j in FragmentMolecule_list:
                    mol = convert_fragment_molecule_to_mol(j, fragment_database)
                    inchi = molecule_to_inchi(mol)
                    f.write('%s %s\n' %(inchi, j.get_build_probability()  ) )

    return FragmentMolecule_list


def get_unique_molecule_list(FragmentMolecule_list, sort_list=True, fragment_database=None):
    
    unique_dict = {}

    check = {}
    mol_check = []

    for idx, i in enumerate(FragmentMolecule_list):
        if i in check:
            if check[i] != i._graph._build_probability:
                print('PROBABILITIES NOT THE SAME', )
                mol_check.append(i)
        else:
            check[i] = i._graph._build_probability


    for idx, i in enumerate(FragmentMolecule_list):
        if i in mol_check:
            print(idx, i._graph._build_probability, i.__hash__())
            mol = convert_fragment_molecule_to_mol(i, fragment_database)
            save_mol_to_sdf(mol, 'different.sdf')



    for i in FragmentMolecule_list:
        if i in unique_dict:
            unique_dict[i]._graph._build_probability += i._graph._build_probability
        else:
            unique_dict[i] = i

    if sort_list is False:
        return unique_dict.keys()

    sorted_list = dict(sorted(unique_dict.items(), key=lambda item: item[1]._graph._build_probability, reverse=True)).keys()

    return sorted_list


def extend_molecule_list_depth_count(FragmentMolecule_list, bond_frequencies, fragment_database_graph, depth, threshold=None):

    if threshold is not None:
        sys.error('Cannot count with threshold')

    for i in range(depth - 1):

        FragmentMolecule_list = extend_molecule_list(FragmentMolecule_list, bond_frequencies, fragment_database_graph, i + 1, threshold)

        print(f'FINAL DEPTH {i+1} TOTAL {len(FragmentMolecule_list)}')

    total = extend_molecule_list_count(FragmentMolecule_list, bond_frequencies, fragment_database_graph, depth)

    print(f'FINAL DEPTH {depth} TOTAL {total}')

    return total

def save_mol_to_sdf(mol, sdffile):

    with open(sdffile, 'a') as f:
        lines = molecule_to_sdf(mol)
        for line in lines:
            f.write(line)
        f.write('$$$$\n')

def save_mol_list_to_sdf(mol_list, sdffile):

    with open(sdffile, 'w') as f:
        print('saving to', sdffile)

    for mol in mol_list:
        save_mol_to_sdf(mol, sdffile)

def read_fragment_database_graph(filename):

    with open(filename) as f:
        lines = f.readlines()

    attach_points_sel = False
    canonical_mapping_sel = False

    attach_points_list = []
    canonical_mapping_list = []

    for line in lines:
        if line.startswith('CANONICAL MAPPING'):
            attach_points_sel = False
        if attach_points_sel is True:
            attach_points_list.append(eval(line))
        if canonical_mapping_sel is True:
            canonical_mapping_list.append(eval(line))
        if line.startswith('ATTACHMENT POINTS'):
            attach_points_sel = True
        if line.startswith('CANONICAL MAPPING'):
            canonical_mapping_sel = True

    f = FragmentGraph()

    assert len(attach_points_list) == len(canonical_mapping_list)

    for i in range(len(attach_points_list)):
        f.add_fragment(i, attach_points_list[i])
        f.fragments[i].set_attribute('frag_id', i)
        f.fragments[i].manual_canonical_mapping(canonical_mapping_list[i])

    return f

def write_fragment_database_graph(fragment_database, filename):

    with open(filename, 'w') as f:
        f.write('ATTACHMENT POINTS\n')
        for i in range(len(fragment_database.fragments)):
            f.write(f'{fragment_database.fragments[i].attachment_points}\n')
        f.write('CANONICAL MAPPING\n')
        for i in range(len(fragment_database.fragments)):
            f.write(f'{fragment_database.fragments[i].get_canonical_mapping()}\n')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Build Molecules using the FragmentMolecule class')
    parser.add_argument('-a','--fragments_sdf', help='SDF file of fragments',required=True)
    parser.add_argument('--atom', type=int, help='Atom to build on parent',required=True)
    parser.add_argument('-d','--frequencies_txt', help='Bond frequencies dictionary in txt file',required=True)
    parser.add_argument('--depth', type=int, help='Depth to build up to',required=True)
    parser.add_argument('--parent_id', type=int, help='Parent id in the fragment database',required=True)    
    
    parser.add_argument('--count', action='store_true', default=False, help='Count total number of molecules without making them', required=False)
    parser.add_argument('-o','--output', help='Output inchi file name', required=False)
    parser.add_argument('-r', '--read_fragment_database', help='Read fragment database from file containing attachment points and canonical mapping', required=False)
    parser.add_argument('-t','--threshold', help='Log10 of build probability threshold of molecules to be built', type=float, required=False)
    parser.add_argument('-w', '--write_fragment_database', help='Write fragment database to file containing attachment points and canonical mapping', required=False)

    with open('different.sdf', 'w') as f:
        print('Writing to different.sdf')

    args = parser.parse_args()

    bond_frequencies = get_bond_frequencies(args.frequencies_txt)
    bond_frequencies = bond_frequencies_to_np(bond_frequencies)

    fragment_database = get_fragment_database(args.fragments_sdf)

    if args.read_fragment_database is not None:
        fragment_database_graph = read_fragment_database_graph(args.read_fragment_database)
    else:    
        fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    if args.write_fragment_database is not None:
        write_fragment_database_graph(fragment_database_graph, args.write_fragment_database)

    parent = FragmentMolecule()

    parent.add_fragment(args.parent_id, [args.atom])

    if args.threshold is not None:
        threshold = 10 ** args.threshold
    else:
        threshold = None

    if args.count:
        extend_molecule_list_depth_count([parent], bond_frequencies, fragment_database_graph, args.depth, args.threshold)

    else:
        output_mol_list = extend_molecule_list_depth([parent], bond_frequencies, fragment_database_graph, depth=args.depth, output=args.output, threshold=threshold)

"""
    if args.output is not None:

        print('writing to', args.output)

        outfile_format = args.output.split('.')[-1].lower()

        if outfile_format == 'sdf':
            output_mol_list_mol = []
            for j in output_mol_list:
                mol = convert_fragment_molecule_to_mol(j, fragment_database)
                output_mol_list_mol.append(mol)
            save_mol_list_to_sdf(output_mol_list_mol, args.output)

        else:
            with open(args.output, 'w') as f:

                for j in output_mol_list:
                    mol = convert_fragment_molecule_to_mol(j, fragment_database)
                    inchi = molecule_to_inchi(mol)
                    f.write('%s %s\n' %(inchi, j.get_build_probability()  ) )
"""

