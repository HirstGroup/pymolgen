#!/usr/bin/env python

import argparse
import ast
import copy
import numpy as np
import os
import random
import sys

from multiprocessing import Pool

from pymolgen.fragment_molecule import *
from pymolgen.generate import SDFDatasetLargeRAM
from pymolgen.molecule_formats import *
from pymolgen.fragment_builder import bond_frequencies_to_np, find_fragment, get_bond_frequencies, get_fragment_database, get_fragment_bond_frequencies_np, map_mols

from functools import partial
print = partial(print, flush=True)

sys.setrecursionlimit(10000)


def extend_molecule_list(FragmentMolecule_list, bond_frequencies, fragment_database_graph, depth=False, threshold=None, version=1):
    """
    Extend a list of molecules by adding one further fragment to each molecule

    Parameters
    ----------
    FragmentMolecule_list : list of FragmentMolecule objects
    bond_frequencies : dict of (int,int) to (int,int):float dict
        Bond frequencies in dict format
    fragment_database_graph : FragmentMolecule object
        Fragment database in FragmentMolecule format
    depth : int, optional
        If set, current depth for generation
    threshold : float, optional
        Build probability threshold for generation
    version : int
        Version for build probability calculation.
        Version 1: calculates build probability factor for new fragment addition as 1/total_free_valence,
        this gives rise to different build probabilities for same molecule if built in different order,
        since adding a fragment with many available attachment points early will have greater effect
        than adding it later on.
        Version 2: calculates build probability factor for new fragment addition as 1/number of fragments * 
        1/number of attachment points for chosen fragment. This gives same build probability for same molecule
        regardless of order of fragment addition.

    Returns
    -------
    output_mol_list : list of FragmentMolecule objects
    """

    output_mol_list = []

    # loop through all molecules
    for f in FragmentMolecule_list:
        free_valence_list = f.list_free_valence_points()
        total_free_valence = f.get_total_free_valence()
        f_build_probability = f.get_build_probability()
        n_fragments = len(f.list_frag_id())

        # loop through fragments in molecule
        for x in range(len(free_valence_list)):

            fragment_id = f.get_frag_id(x)

            # loop through attachment points in each fragment 
            for atom in free_valence_list[x]:

                atom_can = fragment_database_graph.fragments[fragment_id].get_canonical_mapping()[atom]
                fragment_bonds = bond_frequencies[(fragment_id, atom_can)]
                total_freq = sum(fragment_bonds.values())

                for bond, bond_freq in fragment_bonds.items():

                    if version == 1:
                        factor = total_free_valence
                    elif version == 2:
                        factor = len(f.list_free_valence_points()[x]) * n_fragments
                    else:
                        raise Exception(f'Version {version} not implemented')
                    
                    attachment_probability = bond_freq / (total_freq * factor)

                    new_build_probability = attachment_probability * f_build_probability

                    if threshold is not None and new_build_probability < threshold:
                        # do not build molecule if its build probability is below the threshold
                        break

                    j = bond[0]
                    l = bond[1]

                    f2 = copy.deepcopy(f)
                    f2.bp_factor = f.bp_factor * factor

                    node_id = f2.add_fragment(j, fragment_database_graph.fragments[j].attachment_points, fragment_database_graph.fragments[j].get_canonical_mapping())
                    f2.add_bond(x, node_id, atom, l, attachment_probability)

                    if depth is not None:
                        total = len(output_mol_list)

                        if total % 10000 == 0:
                            print(f'DEPTH {depth} TOTAL {total}')

                    output_mol_list.append(f2)

    return output_mol_list


def extend_molecule_random(FragmentMolecule, bond_frequencies, fragment_database_graph, depth=None, version=1):

    f = copy.deepcopy(FragmentMolecule)

    if depth is None:
        depth = 999

    n_fragments = len(f.list_frag_id())

    while f.get_total_free_valence() > 0 and n_fragments <= depth:

        free_valence_list = f.list_free_valence_points()
        total_free_valence = f.get_total_free_valence()
        f_build_probability = f.get_build_probability()

        x = random.randrange(0, len(free_valence_list))

        if len(free_valence_list[x]) == 0:
            continue

        fragment_id = f.get_frag_id(x)

        atom = random.choice(free_valence_list[x])

        atom_can = fragment_database_graph.fragments[fragment_id].get_canonical_mapping()[atom]
        fragment_bonds = bond_frequencies[(fragment_id, atom_can)]
        total_freq = sum(fragment_bonds.values())
        
        bond = random.choices(population=list(fragment_bonds.keys()), weights=fragment_bonds.values(), k=1)[0]

        bond_freq = fragment_bonds[bond]

        if version == 1:
            factor = total_free_valence
        elif version == 2:
            factor = len(f.list_free_valence_points()[x]) * n_fragments
        else:
            raise Exception(f'Version {version} not implemented')
        
        attachment_probability = bond_freq / (total_freq * factor)

        new_build_probability = attachment_probability * f_build_probability

        j = bond[0]
        l = bond[1]

        f.bp_factor = f.bp_factor * factor

        node_id = f.add_fragment(j, fragment_database_graph.fragments[j].attachment_points, fragment_database_graph.fragments[j].get_canonical_mapping())
        f.add_bond(x, node_id, atom, l, attachment_probability)

        n_fragments = len(f.list_frag_id())

        yield f


def extend_molecule_recursive(fragment_molecule, bond_frequencies, fragment_database_graph, threshold):
    """
    Extend single FragmentMolecule recursively up to a specified threshold of build_probability
    """

    f = fragment_molecule

    yield f

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

                if threshold is not None and new_build_probability >= threshold:
                    # do not build molecule if its build probability is below the threshold

                    j = bond[0]
                    l = bond[1]

                    f2 = copy.deepcopy(f)

                    node_id = f2.add_fragment(j, fragment_database_graph.fragments[j].attachment_points, fragment_database_graph.fragments[j].get_canonical_mapping())
                    f2.add_bond(x, node_id, atom, l, attachment_probability)

                    for new_mol in extend_molecule_recursive(f2, bond_frequencies, fragment_database_graph, threshold=threshold):
                        yield new_mol


def extend_molecule_list_count(FragmentMolecule_list, bond_frequencies, fragment_database_graph, depth=None):
    """
    Extend list of FragmentMolecules but only count number of possible molecules without generating them
    """

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


def extend_molecule_list_depth(FragmentMolecule_list, bond_frequencies, fragment_database_graph, depth, fragment_database=None, output=None, parallel=None, restart=None, restart_file=None, return_all=False, saveinchi=False, savesdf=False, sort=False, unique=True, threshold=None, version=1):
    """
    Extend list of FragmentMolecules by adding a further fragment to each molecule up to a certain depth.
    """

    if return_all is True:
        return_all_list = copy.deepcopy(FragmentMolecule_list)

    if restart is not None:
        assert restart_file is not None
        print('Restarting from %s ...' %restart_file)
        FragmentMolecule_list  = read_fragment_molecule_file(restart_file, fragment_database_graph)
        print('Restarting from %s FINISHED' %restart_file)
    else:
        restart = 1

    if depth is None:
        depth = 100

    for i in range(restart, depth + 1):

        if parallel is None:
            # Serial
            FragmentMolecule_list = extend_molecule_list(FragmentMolecule_list, bond_frequencies, fragment_database_graph, i, threshold, version)
            if return_all is True:
                return_all_list.extend(FragmentMolecule_list)

        else:
            # Parallel
            print('Parallel run with %s processes' %parallel)

            args = [(l, bond_frequencies, fragment_database_graph, i, threshold) for l in split_molecule_list(FragmentMolecule_list, parallel)]

            with Pool(processes=parallel) as p:
                extended = p.starmap(extend_molecule_list, args)

            FragmentMolecule_list = []
            for sublist in extended:
                FragmentMolecule_list += sublist

        print(f'FINAL DEPTH {i} TOTAL {len(FragmentMolecule_list)}')

        if sort is True:

            FragmentMolecule_list.sort(reverse=True)

        if unique is True:

            FragmentMolecule_list = get_unique_molecule_list(FragmentMolecule_list, fragment_database=fragment_database)

            print(f'FINAL DEPTH {i} TOTAL UNIQUE {len(FragmentMolecule_list)}')

        if len(FragmentMolecule_list) == 0:
            # stop building molecules
            if return_all is True:
                return return_all_list
            else:
                return FragmentMolecule_list

        if output is not None:

            with open('%s-depth%s.txt' %(output, i), 'w') as f:
                for j in FragmentMolecule_list:
                    f.write(f'{str(j)}\n')

            if saveinchi is True:
                with open('%s-depth%s.inchi' %(output, i), 'w') as f:
                    for j in FragmentMolecule_list:
                        mol = convert_fragment_molecule_to_mol(j, fragment_database)
                        inchi = molecule_to_inchi(mol)
                        f.write('%s %s\n' %(inchi, j.get_build_probability()  ) )

                if savesdf is True:
                    with open('%s-depth%s.sdf' %(output, i), 'w') as f2:
                        for j in FragmentMolecule_list:
                            mol = convert_fragment_molecule_to_mol(j, fragment_database)

                            lines = molecule_to_sdf(mol)

                            for line in lines:
                                f2.write(line)
                            f2.write('$$$$\n')

    if return_all is True:
        return return_all_list
    else:
        return FragmentMolecule_list


def get_unique_molecule_list(FragmentMolecule_list, debug=False, fragment_database=None, sort_list=True):
    """
    Remove duplicates from a FragmentMolecule list
    """

    unique_dict = {}

    # code is generated same molecule with different probabilities, probably due to fragment database having repeated fragments
    if debug is True:
        with open('different.sdf', 'w') as f:
            print('Writing to different.sdf')
        
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

    sorted_list = list(dict(sorted(unique_dict.items(), key=lambda item: item[1]._graph._build_probability, reverse=True)).keys())

    return sorted_list


def split_molecule_list(molecule_list, n):

    """
    Split a molecule list in lists of equal size so that each list adds up to a similar value of build probabilities
    """
    size = len(molecule_list)

    build_probability_list = np.zeros(size)

    for i in range(size):
        build_probability_list[i] = molecule_list[i].get_build_probability()

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


def read_fragment_molecule_file(filename, fragment_database_graph):
    """
    Read a fragment molecule file, used when restarting run
    """

    mol_list = []

    with open(filename) as infile:

        for line in infile:
            # read fragment molecule information
            frag_id_list = [int(x) for x in line.split(':')[0].split('-')]
            bond_list = [ast.literal_eval(x) for x in line.split(':')[1].split(';')]
            build_probability = float(line.strip().split(':')[2])

            f = FragmentMolecule(build_probability)

            # add fragments
            for frag_id in frag_id_list:
                attachment_points = fragment_database_graph.fragments[frag_id].attachment_points
                canonical_mapping = fragment_database_graph.fragments[frag_id].get_canonical_mapping()
                f.add_fragment(frag_id, attachment_points, canonical_mapping)

            # add bonds
            for bond in bond_list:
                f.add_bond(bond[0], bond[1], bond[2], bond[3])

            mol_list.append(f)

    return mol_list


def extend_molecule_list_depth_count(FragmentMolecule_list, bond_frequencies, fragment_database_graph, depth, threshold=None, version=1):
    """
    Extend a FramgentMolecule list up to a certain depth
    """

    if threshold is not None:
        sys.exit('ERROR: Cannot count with threshold')

    for i in range(depth - 1):

        FragmentMolecule_list = extend_molecule_list(FragmentMolecule_list, bond_frequencies, fragment_database_graph, i + 1)

        print(f'FINAL DEPTH {i+1} TOTAL {len(FragmentMolecule_list)}')

    total = extend_molecule_list_count(FragmentMolecule_list, bond_frequencies, fragment_database_graph, depth)

    print(f'FINAL DEPTH {depth} TOTAL {total}')

    return total


def save_mol_to_sdf(mol, sdffile):
    """
    Save molecule to SDF format
    """

    with open(sdffile, 'a') as f:
        lines = molecule_to_sdf(mol)
        for line in lines:
            f.write(line)
        f.write('$$$$\n')


def save_mol_list_to_sdf(mol_list, sdffile):
    """
    Save list of molecules to SDF format
    """

    with open(sdffile, 'w') as f:
        print('saving to', sdffile)

    for mol in mol_list:
        save_mol_to_sdf(mol, sdffile)


def read_fragment_database_graph(filename):
    """
    Read a fragment molecule file in internal FragmentMolecule format
    """

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

    if len(f) == 0:
        raise Exception('ERROR: no fragments read from file', filename)

    print('Reading fragment database graph FINISHED')

    return f


def write_fragment_database_graph(fragment_database, filename):
    """
    Write a fragment molecule file in internal FragmentMolecule format
    """

    print('Writing fragment database graph ...')

    with open(filename, 'w') as f:
        f.write('ATTACHMENT POINTS\n')
        for i in range(len(fragment_database.fragments)):
            f.write(f'{fragment_database.fragments[i].attachment_points}\n')
        f.write('CANONICAL MAPPING\n')
        for i in range(len(fragment_database.fragments)):
            f.write(f'{fragment_database.fragments[i].get_canonical_mapping()}\n')

    print('Writing fragment database graph FINISHED')


def prepare_parent(bond_frequencies, fragment_database, fragment_database_graph, parent_file, parent_fragment_file_list, parent_mapping_1, remove_hydrogens, remove_hydrogens_parent_fragment):
    """
    Prepare parent structure for run
    """

    assert bond_frequencies is not None
    assert fragment_database is not None
    assert parent_fragment_file_list is not None
    assert parent_mapping_1 is not None
    assert remove_hydrogens is not None
    assert remove_hydrogens_parent_fragment is not None

    # convert parent_mapping_1 list to dictionary
    new_dict = {}
    for i in range(0, len(parent_mapping_1), 2):
        new_dict[parent_mapping_1[i]] = parent_mapping_1[i+1]
    parent_mapping_1 = new_dict

    parent_mol = molecule_from_sdf(parent_file)

    attachment_points = []

    # remove hydrogens from parent and determine atoms that will have open valence
    for i in remove_hydrogens:
        parent_mol = parent_mol.remove_atom(i)
        for j in parent_mol.free_valence_list:
            if j not in attachment_points:
                attachment_points.append(j)

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
        parent_fragment_i_dict[attachment_points[i]] = j
        parent_fragment_i_list.append(j)

        lines = molecule_to_sdf(fragment_database[j])

        if j is False:
            sys.exit('Parent fragment not found')

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

    # include parent in fragment_database and fragment_database_graph
    parent_id = len(fragment_database)
    fragment_database.add_mol(parent_mol)
    fragment_database_graph.add_fragment(parent_id, attachment_points)
    fragment_database_graph.fragments[parent_id].set_attribute('frag_id', parent_id)
    fragment_database_graph.fragments[parent_id].set_canonical_mapping(fragment_database)

    # create parent FragmentMolecule object
    parent = FragmentMolecule()
    parent.add_fragment(parent_id, attachment_points, fragment_database_graph.fragments[parent_id].get_canonical_mapping())

    # update bond frequencies for parent
    for attachment_point, equivalent_frag_id in parent_fragment_i_dict.items():
        equivalent_frag = fragment_database[equivalent_frag_id]
        if len(equivalent_frag.attach_points) > 1:
            sys.error('Equivalent fragment cannot have more than 1 attachment point')
        equivalent_frag_atom = equivalent_frag.attach_points[0]
        bond_frequencies[(parent_id, attachment_point)] = bond_frequencies[(equivalent_frag_id, equivalent_frag_atom)]

    return parent, bond_frequencies, fragment_database, fragment_database_graph


def convert_bond_freq_np_to_dict(fragment_database_graph, bond_frequencies, sort_dict=True):
    """
    Convert bond frequencies to dictionary for fast access,
    i.e. fragment id and atom number will be keys
    and available bonds will be values in dictionary
    """

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
    """
    Read bond frequencies from file in dictionary format, 
    i.e. fragment id and atom number will be keys
    and available bonds will be values in dictionary
    """

    print('Reading bond frequencies dict ...')

    bond_frequencies_dict = {}

    with open(infile) as f:
        for line in f:

            key = ast.literal_eval(line.split(':')[0])
            val = ast.literal_eval('{' + line.strip('\n').split(': {')[1])

            bond_frequencies_dict[key] = val

    if len(bond_frequencies_dict) == 0:
        raise Exception('ERROR: no bonds read from file', infile)

    print('Reading bond frequencies dict FINISHED')            

    return bond_frequencies_dict


def write_bond_frequencies_dict(bond_frequencies_dict, outfile):
    """
    Write bond frequencies from file in dictionary format, 
    i.e. fragment id and atom number will be keys
    and available bonds will be values in dictionary
    """

    print('Writing bond frequencies dict to %s ...' %outfile)

    with open(outfile, 'w') as f:
        for key, val in bond_frequencies_dict.items():
            f.write(f'{key}: {val}\n')

    print('Writing bond frequencies dict to %s FINISHED' %outfile)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Build Molecules using the FragmentMolecule class')

    # required arguments
    parser.add_argument('-a','--fragments_sdf', help='SDF file of fragments', required=True)

    # either 1. build from database fragment or 2. build from parent file must be used
    
    # 1. build from database fragment
    parser.add_argument('--atom', type=int, help='Atom to build on parent', required=False)
    parser.add_argument('--parent_id', type=int, help='Parent id in the fragment database',required=False)   

    # 2. build from parent file
    parser.add_argument('-p','--parent_file', help='Parent Structure File in SDF format', required=False)
    parser.add_argument('-x','--parent_fragment_file_list', nargs='+', help='Parent Fragment Structure File list space-separated to search fragment database in SDF format', required=False)
    parser.add_argument('--parent_mapping_1', nargs='+', type=int, help='Parent Fragment i dict list space-separated to search fragment database in SDF format', required=False)
    parser.add_argument('-r','--remove_hydrogens', type=int, nargs='+', help='Space-separated hydrogen atoms that will be created as attachment points, numbered from 0', required=False)
    parser.add_argument('-R','--remove_hydrogens_parent_fragment', type=int, nargs='+', help='Space-separated hydrogen atoms that will be created as attachment points for the parent fragment in database, numbered from 0', required=False)

    # optional arguments
    parser.add_argument('--count', action='store_true', default=False, help='Count total number of molecules without making them', required=False)
    parser.add_argument('-d','--frequencies_txt', help='Bond frequencies dictionary in txt file', required=False)
    parser.add_argument('--depth', type=int, help='Depth to build up to', required=False)
    parser.add_argument('--not_unique', action='store_true', help='Do not obtain unique molecule set in each depth building stage', required=False)
    parser.add_argument('-o','--output', help='Output file name without extension, txt, sdf or inchi extensions will be added', required=False)
    parser.add_argument('--parallel', type=int, help='Number of processes for parallel run', required=False)
    parser.add_argument('-rd', '--read_bond_frequencies_dict', help='Read bond frequencies dict from file', required=False)
    parser.add_argument('-rf', '--read_fragment_database', help='Read fragment database from file containing attachment points and canonical mapping', required=False)
    parser.add_argument('--recursive', action='store_true', help='Build molecules using recursive function, if False extend_molecule_list_depth function will be used instead', required=False)
    parser.add_argument('--random', type=int, help='Build N molecules using random generation, use N = 0 for no maximum number of molecules', required=False)
    parser.add_argument('--restart', type=int, help='Restart from depth', required=False)
    parser.add_argument('--restart_file', help='Restart filename containing molecules built up to restart depth', required=False)
    parser.add_argument('--saveinchi', action='store_true', default=False, help='Save generated molecules as InChi file', required=False)
    parser.add_argument('--savesdf', action='store_true', default=False, help='Save generated molecules as SDF file', required=False)
    parser.add_argument('--seed', type=int, default=False, help='Seed for random run, for test purposes', required=False)
    parser.add_argument('-t','--threshold', help='Log10 of build probability threshold of molecules to be built', type=float, required=False)
    parser.add_argument('--version', default=1, type=int, help='Version for build probability factor, version 1 gives different build probabilities according to the order of fragment addition, version 2 gives same build probabilities for any order', required=False)
    parser.add_argument('-wf', '--write_fragment_database', help='Write fragment database to file containing attachment points and canonical mapping', required=False)
    parser.add_argument('-wd', '--write_bond_frequencies_dict', help='Write bond frequencies dict to file', required=False)

    args = parser.parse_args()

    fragment_database = get_fragment_database(args.fragments_sdf)

    if args.saveinchi is False and args.savesdf is True:
        raise Exception('savesdf option requires saveinchi')

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

    if args.parent_file is not None:
        
        parent, bond_frequencies, fragment_database, fragment_database_graph = prepare_parent(bond_frequencies, fragment_database, fragment_database_graph, args.parent_file, args.parent_fragment_file_list, args.parent_mapping_1, args.remove_hydrogens, args.remove_hydrogens_parent_fragment)

    else:
        parent = FragmentMolecule()

        parent.add_fragment(args.parent_id, [args.atom], {args.atom:args.atom})

    if args.threshold is not None:
        threshold = 10 ** args.threshold
    else:
        threshold = None

    if args.count:
        extend_molecule_list_depth_count([parent], bond_frequencies, fragment_database_graph, args.depth, args.threshold, args.version)

    elif args.recursive is True:

        if args.version != 1:
            raise Exception(f'Cannot run recursive with version {version}')

        if args.output is not None:
            outfile = open(f'{args.output}.txt', 'w')

            if args.saveinchi is True:
                outfile_inchi = open(f'{args.output}.inchi', 'w')

                if args.savesdf is True:
                    outfile_sdf = open(f'{args.output}.sdf', 'w')

        for mol in extend_molecule_recursive(parent, bond_frequencies, fragment_database_graph, threshold=threshold):
            
            if args.output is not None:
                outfile.write(f'{str(mol)}\n')

                if args.saveinchi is True:

                    mol_molecule = convert_fragment_molecule_to_mol(mol, fragment_database)
                    inchi = molecule_to_inchi(mol_molecule)
                    outfile_inchi.write('%s %s\n' %(inchi, mol.get_build_probability() ) )

                    if args.savesdf is True:

                        lines = molecule_to_sdf(mol_molecule)

                        for line in lines:
                            outfile_sdf.write(line)
                        outfile_sdf.write('$$$$\n')            

            else:
                print(mol)

    elif args.random is not None:

        if args.random == 0:
            max_n = np.inf
        else:
            max_n = args.random

        if args.seed is not None:
            random.seed(args.seed)

        if args.output is not None:
            outfile = open(f'{args.output}.txt', 'w')

            if args.saveinchi is True:
                outfile_inchi = open(f'{args.output}.inchi', 'w')

                if args.savesdf is True:
                    outfile_sdf = open(f'{args.output}.sdf', 'w')

        n = 0

        while n < max_n:

            for mol in extend_molecule_random(FragmentMolecule=parent, bond_frequencies=bond_frequencies, fragment_database_graph=fragment_database_graph, depth=args.depth, version=args.version):

                if args.output is not None:

                    outfile.write(f'{str(mol)}\n')

                if args.saveinchi is True:

                    mol_molecule = convert_fragment_molecule_to_mol(mol, fragment_database)
                    inchi = molecule_to_inchi(mol_molecule)
                    outfile_inchi.write('%s %s\n' %(inchi, mol.get_build_probability() ) )

                    if args.savesdf is True:

                        lines = molecule_to_sdf(mol_molecule)

                        for line in lines:
                            outfile_sdf.write(line)
                        outfile_sdf.write('$$$$\n')                    

                n += 1

            if n == args.random:
                break

    else:
        extend_molecule_list_depth([parent], bond_frequencies, fragment_database_graph, depth=args.depth, fragment_database=fragment_database, output=args.output, parallel=args.parallel, restart=args.restart, restart_file=args.restart_file, saveinchi=args.saveinchi, savesdf=args.savesdf, sort=False, threshold=threshold, unique=not args.not_unique, version=args.version)

    print('Normal termination')

