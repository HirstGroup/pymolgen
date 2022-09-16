import sys,os
import random
import numpy as np
import argparse
import subprocess
import time

from multiprocessing import Pool
from functools import partial

import networkx
from networkx.algorithms import isomorphism

from pymolgen.generate import SDFDatasetLargeRAM
from pymolgen.molecule_formats import *
from pymolgen.fragment_mol import print_fragments, get_canonical_mapping, map_mols, get_frag_mapping, update_bond_frequencies
from pymolgen.newmol import WEIGHT_THRESHOLD

print = partial(print, flush=True)

def count_generated_molecules(outfile_name):
    """
    Count the number of generated molecules in an SDF file
    """
    generated_molecules = 0
    with open(outfile_name) as f:
        for line in f:
            if 'V2000' in line:
                generated_molecules += 1
    return generated_molecules

def node_compare_element(node_1, node_2):
    return node_1["element"] == node_2["element"] and node_1["hybridization"] == node_2["hybridization"]

def get_frag_frequencies(frag_frequencies_txt):
    frag_frequencies = []

    with open(frag_frequencies_txt) as infile:
        for line in infile:
            frag_frequencies.extend(int(i) for i in line.split())

    return frag_frequencies

def get_bond_frequencies(bond_frequencies_txt):
    bond_frequencies = {}

    with open(bond_frequencies_txt) as infile:
        for line in infile:
            i = int(line.split()[0].strip('(').strip(','))
            j = int(line.split()[1].strip(','))
            k = int(line.split()[2].strip(','))
            l = int(line.split()[3].strip(':').strip(')').strip(','))
            f = int(line.split()[4])
            bond_frequencies[(i,j,k,l)] = f

    return bond_frequencies

def get_fragment_bond_frequencies(fragment_i, atom_i, bond_frequencies):

    fragment_bond_frequencies = {}

    for key, val in bond_frequencies.items():
        if fragment_i == key[0] and atom_i == key[2]:
            fragment_bond_frequencies[key] = val
        if fragment_i == key[1] and atom_i == key[3]:
            fragment_bond_frequencies[key] = val

    return fragment_bond_frequencies

def get_fragment_bond_frequencies_np(fragment_i, atom_i_can, bond_frequencies_np):

    frag = np.array([fragment_i,atom_i_can])

    key = bond_frequencies_np[0]
    val = bond_frequencies_np[1]

    freq_left = key[:,[0,2]]

    equal_left = freq_left == frag

    equal_left = np.logical_and(equal_left[:,0], equal_left[:,1] )

    freq_right = key[:,[1,3]]

    equal_right = freq_right == frag

    equal_right = np.logical_and(equal_right[:,0], equal_right[:,1] )

    equal = np.logical_or(equal_left, equal_right )
    key_filtered = key[equal]
    val_filered = val[equal]

    return key_filtered, val_filered

def save_neighbours(fragment_i, fragment_bond_frequencies, fragment_database, outfile_name):

    neighbours = []

    with open(outfile_name, 'w') as outfile:
        pass

    for key, val in fragment_bond_frequencies.items():
        if fragment_i == key[0]: 
            neighbours.append(key[1])
        if fragment_i == key[1]:
            neighbours.append(key[0])

    for i in neighbours:

        mol = fragments_database.load_molecule(i)

        lines = molecule_to_sdf(mol)

        with open(outfile_name, 'a') as outfile:
            for line in lines:
                outfile.write(line)
            outfile.write('$$$$\n')

def get_fragment_database(fragments_sdf):

    fragment_database = SDFDatasetLargeRAM(fragments_sdf)

    return fragment_database

def find_fragment(fragment, fragment_database):

    for i in range(len(fragment_database)):

        gm = isomorphism.GraphMatcher(fragment.graph, fragment_database[i].graph, node_match=node_compare_element)

        if gm.is_isomorphic():

            return i

    return False

def get_random_neighbour(fragment_i, fragment_bond_frequencies):

    keys = []
    vals = []

    for key, val in fragment_bond_frequencies.items():
        keys.append(key)
        vals.append(val)

    if len(fragment_bond_frequencies) == 0:
        print('fragment bond frequencies =', fragment_bond_frequencies)
        print('fragment_i =', fragment_i)

    try:
        draw = random.choices(population=keys, weights=vals, k=1)[0]
    except:
        print('keys =', keys)
        print('vals =', vals)
        return None

    if fragment_i == draw[0]:

        new_frag_i = draw[1]
        fragment_i_atom = draw[2]
        new_frag_i_atom = draw[3]

    if fragment_i == draw[1]:

        new_frag_i = draw[0]
        fragment_i_atom = draw[3]
        new_frag_i_atom = draw[2]

    return new_frag_i, new_frag_i_atom

def get_random_neighbour_np(fragment_i, fragment_bond_frequencies):

    keys = fragment_bond_frequencies[0]
    vals = fragment_bond_frequencies[1]

    draw = random.choices(population=keys, weights=vals, k=1)[0]

    if fragment_i == draw[0]:

        new_frag_i = draw[1]
        fragment_i_atom = draw[2]
        new_frag_i_atom = draw[3]

    if fragment_i == draw[1]:

        new_frag_i = draw[0]
        fragment_i_atom = draw[3]
        new_frag_i_atom = draw[2]

    return new_frag_i, new_frag_i_atom

def get_length(list):

    length = 0

    for i in list:
        length += len(i)

    return length

def reverse_canonical_mapping(fragment):
    gm = isomorphism.GraphMatcher(fragment, fragment, node_match=node_compare_element)

    all_mappings = []

    for mapping in gm.isomorphisms_iter():
        all_mappings.append(mapping)

    canonical_mapping = all_mappings[0]

    for i in all_mappings:
        for key, val in i.items():
            if canonical_mapping[key] > val:
                canonical_mapping[key] = val

    return canonical_mapping

def bond_frequencies_to_np(bond_frequencies):

    n = len(bond_frequencies)

    a = np.zeros((n,4), dtype=int)

    b = np.zeros(n, dtype=int)

    n = 0
    for key, val in bond_frequencies.items():
        a[n] = np.array(key)
        b[n] = val
        n += 1

    return a, b

def build_molecule(fragments_sdf, fragments_txt, frequencies_txt, parent_file, parent_fragment_file_list, parent_mapping_1,  remove_hydrogens, remove_hydrogens_parent_fragment, n_mol=None, filters=False, unique=False, figure=None, rules=False, rules_file=None, restart=False, verbose=False, mw_check=False, use_numpy=True, batch_size=None, cpu=1, candidate_file=None, cap=False, intermediates=False):

    if filters:
        from pymolgen.newmol import filters_final_mol, filters_final_mol_return_mol

    if batch_size is None:
        batch_size = 1

    new_dict = {}

    for i in range(0, len(parent_mapping_1), 2):
        new_dict[parent_mapping_1[i]] = parent_mapping_1[i+1]

    parent_mapping_1 = new_dict

    pains_database = None
    # build pains_database if using filters
    if filters:
        from pymolgen.newmol import gen_pains_database
        try:
            pains_database = gen_pains_database()
        except:
            raise Exception("Could not generate pains database")

    #make databases and update atom numberings
    fragment_database = get_fragment_database(fragments_sdf)
    frag_mapping = get_frag_mapping(fragments_txt)
    bond_frequencies = get_bond_frequencies(frequencies_txt)   
    bond_frequencies = update_bond_frequencies(bond_frequencies, frag_mapping)

    if use_numpy:
        bond_frequencies = bond_frequencies_to_np(bond_frequencies)

    if unique:
        candidate_list = set()
    else:
        candidate_list = None

    parent_mol = molecule_from_sdf(parent_file)

    attachment_points = []

    for i in remove_hydrogens:
        parent_mol = parent_mol.remove_atom(i)
        for j in parent_mol.free_valence_list:
            if j not in attachment_points:
                attachment_points.append(j)

    parent_mw = Molecule.molecular_weight(parent_mol)
    parent_fragment_list = []

    for i in parent_fragment_file_list:
        parent_fragment_list.append(molecule_from_sdf(i))

    for i in range(len(parent_fragment_list)):
        parent_fragment_list[i] = parent_fragment_list[i].remove_atom(remove_hydrogens_parent_fragment[i])

    parent_fragment_original_list = []

    for i in parent_fragment_list:
        parent_fragment_original_list.append(i)

    parent_fragment_i_list = []

    new_dict = {}
    for i in range(len(parent_fragment_list)):
        j = find_fragment(parent_fragment_list[i], fragment_database)
        new_dict[attachment_points[i]] = j
        parent_fragment_i_list.append(j)

        lines = molecule_to_sdf(fragment_database[j])

        with open('parent_fragment%s.sdf' %i, 'w') as outfile:
            for line in lines:
                outfile.write(line)

            outfile.write('$$$$\n')

        if j is False:
            sys.exit('Parent fragment not found')

    parent_fragment_i_dict = new_dict

    parent_fragment_list = []

    for i in parent_fragment_i_list:
        parent_fragment_list.append(fragment_database[i])

    parent_mapping_2 = []

    for i in range(len(parent_fragment_list)):
        parent_mapping_2.append(map_mols(parent_fragment_original_list[i].graph, parent_fragment_list[i].graph))

    parent_mapping = {}
    n = 0
    for key, val in parent_mapping_1.items():
        parent_mapping[key] = parent_mapping_2[n][val]
        n += 1

    if restart is False:
        n = 0
    else:
        n = count_generated_molecules(outfile_name)

    if figure is not None:
        with open(figure, 'w') as outfile:
            print('Writing to figure', figure)

    output_mol_list = []

    if n_mol is None:
        n_mol = np.inf

    build_mol_single_partial = partial(build_mol_single,parent_mol, parent_fragment_list, parent_fragment_i_list, parent_fragment_i_dict, fragment_database, bond_frequencies, parent_mapping, filters, pains_database, candidate_list, figure, verbose, mw_check, use_numpy, cap, intermediates)

    start_time = time.time()
    current_time = start_time

    if candidate_file is not None:
        with open(candidate_file, 'w') as outfile:
            print('Writing candidates to', candidate_file)

    while n < n_mol:

        p = Pool(processes=cpu)

        size = min(batch_size, n_mol - n)

        output_mol_list_parallel = p.map(build_mol_single_partial, range(size) )

        for i in output_mol_list_parallel:
            if i is not None:
                output_mol_list.extend(i)

        p.close()

        output_mol_list = [i for i in output_mol_list if i is not None]

        if candidate_list is not None:

            output_mol_list, new_inchi_set = unique_mol_list(output_mol_list)

            candidate_list.update(new_inchi_set)

            if candidate_file is not None:
                with open(candidate_file, 'a') as outfile:
                    for inchi in new_inchi_set:
                        outfile.write('%s\n' %inchi)


        if len(output_mol_list) > 0:
            previous_time = current_time
            current_time = time.time() - start_time
            interval_time = (current_time - previous_time) / len(output_mol_list)
            print('TIME %.2f' %interval_time)

        if len(output_mol_list) > min(batch_size, n_mol - n):
                    
            if filters:

                filters_final_mol_return_mol_partial = partial(filters_final_mol_return_mol, pains_database)

                p = Pool(processes=cpu)

                new_output_mol_list = p.map(filters_final_mol_return_mol_partial, output_mol_list )

                p.close()

                output_mol_list = [i for i in new_output_mol_list if i is not None]

            if rules:
                output_mol_list = rules_batch(output_mol_list, rules_file)                

            for mol in output_mol_list:

                n += 1

                if verbose:
                    smi = molecule_to_smiles(mol)
                    mw = mol.molecular_weight()
                    print('NEW_CANDIDATE %s %s %.1f' % (n, smi, mw))            
                else:
                    print('NEW_CANDIDATE %s' %n )

            yield output_mol_list

            output_mol_list = []

def unique_mol_list(mol_list):

    inchi_set = set()

    output_mol_list = []

    for mol in mol_list:
        inchi = molecule_to_inchi(mol)
        if inchi not in inchi_set:
            inchi_set.add(inchi)
            output_mol_list.append(mol)
        else:
            print('Not unique', inchi)

    return output_mol_list, inchi_set



def rules_batch(output_mol_list, rules_file):

    n = 0

    with open(rules_file, 'w') as outfile:
        pass

    for mol in output_mol_list:

        smi = molecule_to_smiles(mol)
        with open(rules_file, 'a') as outfile:
            outfile.write('%s %s\n' %(smi, n) )
        n += 1

    home = os.path.expanduser('~/')

    result = subprocess.run([home + 'Lilly-Medchem-Rules/Lilly_Medchem_Rules.rb %s' %rules_file], shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8')

    new_output_mol_list = []

    for line in result.split('\n'):
        if not line.strip():
            continue
        i_mol = int(line.split()[1])

        new_output_mol_list.append(output_mol_list[i_mol])

    return new_output_mol_list

def build_mol_single_batch(parent_mol, parent_fragment_list, parent_fragment_i_list, parent_fragment_i_dict, fragment_database, bond_frequencies, parent_mapping, filters=False, pains_database=None, candidate_list=None, figure=None, verbose=False, mw_check=False, use_numpy=True, batch_size=None):

    if batch_size is None:
        batch_size = 1

    output_mol_list = []

    while len(output_mol_list) < batch_size:

        mol = build_mol_single(parent_mol, parent_fragment_list, parent_fragment_i_list, parent_fragment_i_dict, fragment_database, bond_frequencies, parent_mapping, filters, pains_database, candidate_list, figure, verbose, mw_check, use_numpy)

        if mol is not None:

            if candidate_list is not None:

                inchi = molecule_to_inchi(mol)

                if inchi not in candidate_list:
                    output_mol_list.append(mol)

            else:
                output_mol_list.append(mol)

    return output_mol_list        


def build_mol_single(parent_mol, parent_fragment_list, parent_fragment_i_list, parent_fragment_i_dict, fragment_database, bond_frequencies, parent_mapping, filters=False, pains_database=None, candidate_list=None, figure=None, verbose=False, mw_check=False, use_numpy=True, cap=False, intermediates=False, dummy=None):

    return_mol_list = []

    #prepare parent fragment
    frag_list = []
    frag_mol_list = [parent_mol]
    frag_bond_list = []
    frag_free_valence_list = []

    frag_free_valence_list.append([])

    for i in parent_mol.free_valence_list:
        frag_free_valence_list[0].append(i)

    #for i in range(len(parent_fragment_i_list)):
    #    frag_list.append(-(i+1))
    frag_list.append(-1)

    if mw_check:
        mw = parent_mol.molecular_weight()

    counter = 0
    while get_length(frag_free_valence_list) != 0:

        counter += 1
        if counter == 1000:
            if verbose: print('MAX counter')
            return None

        # choose random position of constituent fragments in molecule
        i = random.randrange(len(frag_list))

        # if chosen fragment has free valence points
        if len(frag_free_valence_list[i]) > 0:

            j = len(frag_list)

            # get atom from fragment_i to build on
            atom_i = random.choice(frag_free_valence_list[i])
            # get fragment_i (index in fragment_database)
            fragment_i = frag_list[i]

            if fragment_i == -1:

                fragment_i = parent_fragment_i_dict[atom_i]
                
                # get mol for fragment_i
                fragment_i_mol = fragment_database[fragment_i]

                # get mapped atom_i since fragment_bond_frequencies are stored for canonical atoms
                atom_i_can = parent_mapping[atom_i]

            else:
                # get mol for fragment_i
                fragment_i_mol = fragment_database[fragment_i]

                # get canonical mapping
                canonical_mapping = get_canonical_mapping(fragment_i_mol.graph)

                # get mapped atom_i since fragment_bond_frequencies are stored for canonical atoms
                atom_i_can = canonical_mapping[atom_i]

            # get bond frequencies for fragment_i
            if use_numpy:
                fragment_bond_frequencies = get_fragment_bond_frequencies_np(fragment_i, atom_i_can, bond_frequencies)

                # return none molecule if fragment_bond_frequencies has length 0 (cannot build on fragment)
                # this shouldn't happen since all fragments come from molecules so they shuold all have bonds
                # but there could be errors in the database
                if len(fragment_bond_frequencies[0]) == 0:
                    if verbose: print('len(fragment_bond_frequencies[0]) == 0')
                    return None

                # choose random neighbour
                get_random_neighbour_out = get_random_neighbour_np(fragment_i, fragment_bond_frequencies)

            else:
                fragment_bond_frequencies = get_fragment_bond_frequencies(fragment_i, atom_i_can, bond_frequencies)

                # return none molecule if fragment_bond_frequencies has length 0 (cannot build on fragment)
                # this shouldn't happen since all fragments come from molecules so they shuold all have bonds
                # but there could be errors in the database
                if len(fragment_bond_frequencies) == 0:
                    if verbose:
                        print('fragment_bond_frequencies = 0')
                        print(fragment_bond_frequencies)
                    return None

                # choose random neighbour
                get_random_neighbour_out = get_random_neighbour(fragment_i, fragment_bond_frequencies)

            if get_random_neighbour_out is not None:
                new_frag_i = get_random_neighbour_out[0]
                new_frag_i_atom = get_random_neighbour_out[1]
            else:
                if verbose:
                    print('get_random_neighbour_out is None')
                return None

            # generate molecule object from new_frag_i
            new_frag = fragment_database[new_frag_i]

            # get free valence points of new fragment
            new_free_valence_list = new_frag.free_valence_list

            if not new_frag.is_fluorine() and new_frag_i_atom in new_free_valence_list:

                if mw_check:
                    mw += new_frag.molecular_weight()
                    if mw > WEIGHT_THRESHOLD:
                        if cap:
                            for i in frag_list[1:]:
                                frag_mol_list.append(fragment_database[i])

                            mol = combine_all_fragments(frag_mol_list, frag_list, frag_bond_list)

                            mol.hydrogenate()

                            if candidate_list is not None:
                                inchi = molecule_to_inchi(mol)
                                if inchi in candidate_list:
                                    if verbose: print('Not unique')
                                    return None

                            return_mol_list.append(mol)

                            return return_mol_list
                        else:    
                            if verbose: print('Failed mw_check')
                            return None

                # add neighbour index in fragment_database to frag_list
                frag_list.append(new_frag_i)

                # add bond betweent current fragment and new fragment to list of bonds between fragments (frag_bond_list)
                frag_bond_list.append((i, j, atom_i, new_frag_i_atom))

                # remove atom from current fragment making bond to new fragment from frag_free_valence_list[i]
                frag_free_valence_list[i].remove(atom_i)

                # remove atom from new fragment making bond to current fragment from new fragment's list of free valence points
                try: new_free_valence_list.remove(new_frag_i_atom)
                except: 
                    smi = molecule_to_smiles(new_frag)
                    print(smi)
                    print_molecule(new_frag)
                    lines = molecule_to_sdf(new_frag)
                    with open('new_frag.sdf', 'w') as outfile:
                        for line in lines:
                            outfile.write(line)
                    print(lines)
                    raise Exception('Could not remove atom', new_frag_i_atom, 'from', new_free_valence_list, 'fragment_i =', fragment_i, 'new_frag_i=', new_frag_i)
                    return None

                # add new_free_valence_list to the list of available valence points in molecule being built
                frag_free_valence_list.append(new_free_valence_list)

                if intermediates:

                    frag_mol_list_int = [i for i in frag_mol_list]

                    for i in frag_list[1:]:
                        frag_mol_list_int.append(fragment_database[i])

                    mol = combine_all_fragments(frag_mol_list_int, frag_list, frag_bond_list)

                    mol.hydrogenate()

                    if candidate_list is not None:
                        inchi = molecule_to_inchi(mol)
                        if inchi not in candidate_list:
                            return_mol_list.append(mol)
                            print('INTERMEDIATE %s' %inchi)

                    else:
                        return_mol_list.append(mol)

    for i in frag_list[1:]:
        frag_mol_list.append(fragment_database[i])

    mol = combine_all_fragments(frag_mol_list, frag_list, frag_bond_list)

    if candidate_list is not None:
        inchi = molecule_to_inchi(mol)
        if inchi in candidate_list:
            if verbose: print('Not unique')
            return None

    return_mol_list.append(mol)

    return return_mol_list

def combine_all_fragments(frag_mol_list, frag_list, frag_bond_list):

    mol = Molecule()

    new_frag_bond_list = []

    frag_len_list = [len(i.graph.nodes) for i in frag_mol_list]

    added_frag_len_list = [0]

    for i in range(1,len(frag_len_list)):
        added_frag_len_list.append(sum(frag_len_list[:i]))

    for bond in frag_bond_list:
        i = bond[0]
        j = bond[1]
        k = bond[2]
        l = bond[3]

        k += added_frag_len_list[i]
        l += added_frag_len_list[j]

        new_frag_bond_list.append((i,j,k,l))

    graphs = [x.graph for x in frag_mol_list]

    mol.graph = networkx.disjoint_union_all(graphs)

    for bond in new_frag_bond_list:
        k = bond[2]
        l = bond[3]
        mol.graph.add_edge(k, l, order=1)        

    return mol

def fragment_builder(fragments_sdf, fragments_txt, frequencies_txt, parent_file, parent_fragment_file_list, parent_mapping_1,  remove_hydrogens, remove_hydrogens_parent_fragment, outfile_name, n_mol=None, filters=False, unique=False, figure=None, rules=False, rules_file=None, restart=False, verbose=False, mw_check=False, use_numpy=True, batch_size=None, cpu=1, candidate_file=None, cap=False, intermediates=False):

    if outfile_name is not None:
        with open(outfile_name, 'w') as outfile:
            print('Writing to', outfile_name)

    for mol_list in build_molecule(fragments_sdf=fragments_sdf, fragments_txt=fragments_txt, frequencies_txt=frequencies_txt, parent_file=parent_file, parent_fragment_file_list=parent_fragment_file_list, parent_mapping_1=parent_mapping_1, remove_hydrogens=remove_hydrogens, remove_hydrogens_parent_fragment=remove_hydrogens_parent_fragment, n_mol=n_mol, unique=unique, rules=rules, rules_file=rules_file, filters=filters, restart=restart, verbose=verbose, mw_check=mw_check, use_numpy=use_numpy, batch_size=batch_size, cpu=cpu, candidate_file=candidate_file, cap=cap, intermediates=intermediates):

        for mol in mol_list:

            lines = molecule_to_sdf(mol)

            with open(outfile_name, 'a') as outfile:
                for line in lines:
                    outfile.write(line)

                outfile.write('$$$$\n')

            if figure is not None:

                newatoms = []
                for i in mol.graph.nodes:
                    if i >= 44:
                        newatoms.append(i)

                fig = mol.get_fragment(newatoms)
                #fig.hydrogenate()
                smi = molecule_to_smiles(fig)
                #print('ATTACHED ', smi)
                print_molecule(fig)

                lines = molecule_to_sdf(fig)

                with open(figure, 'a') as outfile:
                    for line in lines:
                        outfile.write(line)

                    outfile.write('$$$$\n')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Pymolgen molecular generator from fragments')
    parser.add_argument('-a','--fragments_sdf', help='SDF file of fragments',required=True)
    parser.add_argument('-f','--fragments_txt', help='List of fragments in TXT file',required=True)
    parser.add_argument('-d','--frequencies_txt', help='Bond frequencies dictionary in txt file',required=True)
    parser.add_argument('-p','--parent_file', help='Parent Structure File in SDF format',required=True)
    parser.add_argument('-x','--parent_fragment_file_list', nargs='+', help='Parent Fragment Structure File list space-separated to search fragment database in SDF format',required=True)
    parser.add_argument('--parent_mapping_1', nargs='+', type=int, help='Parent Fragment i dict list space-separated to search fragment database in SDF format',required=True)
    parser.add_argument('-r','--remove_hydrogens', type=int, nargs='+', help='Space-separated hydrogen atoms that will be created as attachment points, numbered from 0',required=False)
    parser.add_argument('-R','--remove_hydrogens_parent_fragment', type=int, nargs='+', help='Space-separated hydrogen atoms that will be created as attachment points for the parent fragment in database, numbered from 0',required=True)
    parser.add_argument('-s','--seed', type=int, help='Seed for random number generator',required=False)
    parser.add_argument('-o','--outfile_name', help='Output File Name',required=True)
    parser.add_argument('-n','--n_mol', type=int, help='Number of molecules to generate',required=False)
    parser.add_argument('--unique', action='store_true', help='Generate unique set of molecules', required=False)
    parser.add_argument('--rules', action='store_true', help='Use rules to filter', required=False)
    parser.add_argument('--rules_file', help='Rules file name for rules to filter', required=False)
    parser.add_argument('--filters', action='store_true', help='Use filters', required=False)
    parser.add_argument('--restart', action='store_true', help='Restart generation from previous run')
    parser.add_argument('--verbose', action='store_true', help='Verbose output', required=False)
    parser.add_argument('--mw_check', action='store_true', help='MW filter in every fragment addition')
    parser.add_argument('--no_numpy', action='store_true', help='Do not use numpy for fragment bond frequencies')
    parser.add_argument('--batch_size', type=int, help='Batch size for rules', required=False)
    parser.add_argument('--cpu', type=int, help='Number of processes in parallel', default=1, required=False)
    parser.add_argument('--candidate_file', help='Candidate file to save all molecules generated as inchi', required=False)
    parser.add_argument('--cap', action='store_true', help='Cap intermediate molecules if new fragment goes over mass budget', required=False)
    parser.add_argument('--intermediates', action='store_true', help='Save intermediate molecules while constructing new ones', required=False)


    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)
        if args.n_mol is None:
            sys.exit('Cannot run with seed and n_mol infinite')
        elif args.n_mol > 10000:
            sys.exit('Cannot run with seed and n_mol > 1000')

    use_numpy = not args.no_numpy

    fragment_builder(fragments_sdf=args.fragments_sdf, fragments_txt=args.fragments_txt, frequencies_txt=args.frequencies_txt, parent_file=args.parent_file, parent_fragment_file_list=args.parent_fragment_file_list, parent_mapping_1=args.parent_mapping_1, remove_hydrogens=args.remove_hydrogens, remove_hydrogens_parent_fragment=args.remove_hydrogens_parent_fragment,outfile_name=args.outfile_name, n_mol=args.n_mol, unique=args.unique, rules=args.rules, rules_file=args.rules_file, filters=args.filters, restart=args.restart, verbose=args.verbose, mw_check=args.mw_check, use_numpy=use_numpy, batch_size=args.batch_size, cpu=args.cpu, candidate_file=args.candidate_file, cap=args.cap, intermediates=args.intermediates)

    print('Normal termination')








