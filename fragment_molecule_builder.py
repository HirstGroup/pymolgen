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


def extend_molecule_list(FragmentMolecule_list, bond_frequencies, fragment_database_graph, depth=False, threshold=None):

    output_mol_list = []

    # loop through all molecules
    for f in FragmentMolecule_list:
        free_valence_list = f.list_free_valence_points()
        total_free_valence = f.get_total_free_valence()

        # loop through fragments in molecule
        for x in range(len(free_valence_list)):

            fragment_id = f.get_frag_id(x)

            # loop through attachment points in each fragment 
            for atom in free_valence_list[x]:

                atom_can = fragment_database_graph.fragments[fragment_id].get_canonical_mapping()[atom]
                fragment_bonds = bond_frequencies[(fragment_id, atom_can)]
                total_freq = sum(fragment_bonds.values())

                for bond, bond_freq in fragment_bonds.items():

                    attachment_probability = bond_freq / ( total_freq * total_free_valence)

                    new_build_probability = attachment_probability * f.get_build_probability()

                    if threshold is not None and new_build_probability < threshold:
                        # do not build molecule if its build probability is below the threshold
                        continue

                    j = bond[0]
                    l = bond[1]

                    f2 = copy.deepcopy(f)

                    node_id = f2.add_fragment(j, fragment_database_graph.fragments[j].attachment_points, fragment_database_graph.fragments[j].get_canonical_mapping())
                    f2.add_bond(x, node_id, atom, l, attachment_probability)

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

                fragment_bonds = bond_frequencies[(fragment_id, atom_can)]

                total += len(fragment_bonds)

                if total % 10000 == 0:
                    print(f'DEPTH {depth} TOTAL {total}')

    return total


def extend_molecule_list_depth(FragmentMolecule_list, bond_frequencies, fragment_database_graph, depth, fragment_database=None, output=None, unique=True, threshold=None):

    for i in range(depth):

        FragmentMolecule_list = extend_molecule_list(FragmentMolecule_list, bond_frequencies, fragment_database_graph, i + 1, threshold)

        print(f'FINAL DEPTH {i+1} TOTAL {len(FragmentMolecule_list)}')

        if unique is True:

            FragmentMolecule_list = get_unique_molecule_list(FragmentMolecule_list, fragment_database=fragment_database)

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

    with open('different.sdf', 'w') as f:
        print('Writing to different.sdf')
    
    unique_dict = {}
    if fragment_database is not None:
        check = {}
        mol_check = []

        for idx, i in enumerate(FragmentMolecule_list):
            if i in check:
                if check[i] != i._graph._build_probability:
                    print('PROBABILITIES NOT THE SAME')
                    mol_check.append(i)
            else:
                check[i] = i._graph._build_probability

        new_check = []

        for i in check.keys():
            for j in FragmentMolecule_list:
                if i == j:
                    new_check.append(j)
                    print('CHECK', i.__hash__(), i.get_build_probability(), j.get_build_probability())

        for i in new_check:
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

    print('Reading fragment database graph ...')

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

    print('Reading fragment database graph FINISHED')

    return f

def write_fragment_database_graph(fragment_database, filename):

    print('Writing fragment database graph ...')

    with open(filename, 'w') as f:
        f.write('ATTACHMENT POINTS\n')
        for i in range(len(fragment_database.fragments)):
            f.write(f'{fragment_database.fragments[i].attachment_points}\n')
        f.write('CANONICAL MAPPING\n')
        for i in range(len(fragment_database.fragments)):
            f.write(f'{fragment_database.fragments[i].get_canonical_mapping()}\n')

    print('Writing fragment database graph FINISHED')

def prepare_parent():

    attachment_points = []

    # remove hydrogens from parent and determine atoms that will have open valence
    for i in remove_hydrogens:
        parent_mol = parent_mol.remove_atom(i)
        for j in parent_mol.free_valence_list:
            if j not in attachment_points:
                attachment_points.append(j)

    parent_mw = Molecule.molecular_weight(parent_mol)

    # make list of equivalent fragments to build on parent
    parent_fragment_list = [molecule_from_sdf(x) for x in parent_fragment_file_list]

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
        print(attachment_points); print('j =', j)
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

def convert_bond_freq_np_to_dict(fragment_database_graph, bond_frequencies, sort_dict=True):

    print('Converting bond frequencies to dictionary ...')

    bond_frequencies_dict = {}

    for frag_id in range(len(fragment_database_graph.fragments)):

        fragment = fragment_database_graph.fragments[frag_id]

        attachment_points_can = [fragment.get_canonical_mapping()[x] for x in fragment.attachment_points]

        attachment_points_can_sorted = []

        for i in attachment_points_can:
            if i not in attachment_points_can_sorted:
                attachment_points_can_sorted.append(i)

        for atom in sorted(set(attachment_points_can)):

            frag_id_atom_dict = {}

            atom_can = fragment.get_canonical_mapping()[atom]

            bonds, freq = get_fragment_bond_frequencies_np(frag_id, atom_can, bond_frequencies)

            if sort_dict is True:
                # sort in descending order
                sort_index = np.argsort(-freq)

                bonds = [bonds[i] for i in sort_index]

                freq = [freq[i] for i in sort_index]

            for x in range(len(bonds)):
                i = bonds[x][0]
                j = bonds[x][1]
                k = bonds[x][2]
                l = bonds[x][3]

                if i == frag_id and k == atom_can:
                    frag_id_atom_dict[j,l] = freq[x]

                elif j == frag_id and l == atom_can:
                    frag_id_atom_dict[i,k] = freq[x]

                else:
                    print('frag_id %s and atom %s not found' %(frag_id, atom))

            bond_frequencies_dict[frag_id, atom_can] = frag_id_atom_dict

    print('Converting bond frequencies to dictionary FINISHED')

    return bond_frequencies_dict


def read_bond_frequencies_dict(infile):

    print('Reading bond frequencies dict ...')

    bond_frequencies_dict = {}

    with open(infile) as f:
        for line in f:

            key = eval(line.split(':')[0])
            val = eval('{' + line.strip('\n').split(': {')[1])

            bond_frequencies_dict[key] = val

    print('Reading bond frequencies dict FINISHED')            

    return bond_frequencies_dict

def write_bond_frequencies_dict(bond_frequencies_dict, outfile):

    print('Writing bond frequencies dict to %s ...' %outfile)

    with open(outfile, 'w') as f:
        for key, val in bond_frequencies_dict.items():
            f.write(f'{key}: {val}\n')

    print('Writing bond frequencies dict to %s FINISHED' %outfile)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Build Molecules using the FragmentMolecule class')
    parser.add_argument('-a','--fragments_sdf', help='SDF file of fragments',required=True)
    parser.add_argument('--atom', type=int, help='Atom to build on parent',required=True)
    parser.add_argument('--depth', type=int, help='Depth to build up to',required=True)
    parser.add_argument('--parent_id', type=int, help='Parent id in the fragment database',required=True)    
    
    parser.add_argument('--count', action='store_true', default=False, help='Count total number of molecules without making them', required=False)
    parser.add_argument('-d','--frequencies_txt', help='Bond frequencies dictionary in txt file',required=False)
    parser.add_argument('-o','--output', help='Output inchi file name', required=False)
    parser.add_argument('-r', '--read_fragment_database', help='Read fragment database from file containing attachment points and canonical mapping', required=False)
    parser.add_argument('-rd', '--read_bond_frequencies_dict', help='Read bond frequencies dict from file', required=False)
    parser.add_argument('-t','--threshold', help='Log10 of build probability threshold of molecules to be built', type=float, required=False)
    parser.add_argument('-w', '--write_fragment_database', help='Write fragment database to file containing attachment points and canonical mapping', required=False)
    parser.add_argument('-wd', '--write_bond_frequencies_dict', help='Write bond frequencies dict to file', required=False)

    args = parser.parse_args()

    fragment_database = get_fragment_database(args.fragments_sdf)

    if args.read_fragment_database is not None:
        fragment_database_graph = read_fragment_database_graph(args.read_fragment_database)
    else:    
        fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    if args.write_fragment_database is not None:
        write_fragment_database_graph(fragment_database_graph, args.write_fragment_database)

    if args.read_bond_frequencies_dict is None:
        bond_frequencies = get_bond_frequencies(args.frequencies_txt)
        bond_frequencies = bond_frequencies_to_np(bond_frequencies)
        bond_frequencies = convert_bond_freq_np_to_dict(fragment_database_graph, bond_frequencies)
    else:
        bond_frequencies = read_bond_frequencies_dict(args.read_bond_frequencies_dict)

    if args.write_bond_frequencies_dict is not None:
        write_bond_frequencies_dict(bond_frequencies, args.write_bond_frequencies_dict)

    parent = FragmentMolecule()

    parent.add_fragment(args.parent_id, [args.atom], {args.atom:args.atom})

    if args.threshold is not None:
        threshold = 10 ** args.threshold
    else:
        threshold = None

    if args.count:
        extend_molecule_list_depth_count([parent], bond_frequencies, fragment_database_graph, args.depth, args.threshold)

    else:
        output_mol_list = extend_molecule_list_depth([parent], bond_frequencies, fragment_database_graph, depth=args.depth, output=args.output, threshold=threshold)


