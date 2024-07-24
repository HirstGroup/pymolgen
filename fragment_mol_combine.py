import sys,os
import numpy as np

import networkx

from pymolgen.fragment_mol import *
from pymolgen.fragment_builder import *

from functools import partial
print = partial(print, flush=True)

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


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


def remove_bond_frequencies(bond_frequencies, fragment_list):
    """
    Remove bond frequencies when removing framgents from a database)

    fragment_list: list of fragments to keep
    """

    fragment_list = set(fragment_list)

    d = {}

    for key, val in bond_frequencies.items():

        i = key[0]
        j = key[1]

        k = key[2]
        l = key[3]

        if i in fragment_list and j in fragment_list:

            d[(i,j,k,l)] = val

    return d


def remove_bond_frequencies_halogen(bond_frequencies, aromatic_list, halogen_list, verbose=False):
    """
    Remove aliphatic halogen bonds

    Parameters
    ----------
    bond_frequencies : dict
        Dictionary mapping (i,j,k,l) bond to bond frequency
    aromatic_list: list of aromatic fragments
    halogen_list: list of halogen fragments
    verbose : bool, optional
        If true, print verbose output
    """

    aromatic_list = set(aromatic_list)
    halogen_list = set(halogen_list)

    d = {}

    for key, val in bond_frequencies.items():

        i = key[0]
        j = key[1]

        k = key[2]
        l = key[3]

        if i in halogen_list:
            if j not in aromatic_list:
                if verbose and val > 1000: print('HALOGEN OUT', i, j, k, l, val)
                continue
            else:
                if verbose and val > 1000: print('HALOGEN IN', i, j, k, l, val)

        if j in halogen_list:
            if i not in aromatic_list:
                if verbose and val > 1000: print('HALOGEN OUT', i, j, k, l, val)
                continue
            else:
                if verbose and val > 1000: print('HALOGEN IN', i, j, k, l, val)

            d[(i,j,k,l)] = val

    return d


def list_elements(mol):

    elements = set()

    for i in mol.graph.nodes:
        elements.add(mol.graph.nodes[i]["element"])

    return elements


def filter_database(fragment_database_mol, inchi_filter=None, pains=False):
    """
    Filter database by removing non-wanted fragments
    Writes several sdf files with the classification of fragments

    Parameters
    ----------
    fragment_database_mol : list of molecule objects
        List of fragments represented as molecule objects
    inchi_filter : str, optional
        Inchi file with structures to remove
    pains : bool, optional
        If true, apply pains filters

    Returns
    -------
    filter_list : list of molecule objects
        List of fragment to keep, represented as molecule objects
    """

    elements = set()

    filter_list = []

    with open('filter.sdf', 'w') as outfile:
        print('Writing fragments filtered out to filter.sdf')

    with open('aliphatic3.sdf', 'w') as outfile:
        print('Writing fragments containing 3 aliphatic rings to aliphatic3.sdf')

    with open('aliphatic4.sdf', 'w') as outfile:
        print('Writing fragments containing 4 aliphatic rings to aliphatic4.sdf')

    with open('cycles.sdf', 'w') as outfile:
        print('Writing large cycles to cycles.sdf')

    with open('cages3.sdf', 'w') as outfile:
        print('Writing fragments containing cages and 3 aliphatic rings to cages3.sdf')

    with open('cages4.sdf', 'w') as outfile:
        print('Writing fragments containing cages and 4 aliphatic rings to cages4.sdf')

    with open('sulfone.sdf', 'w') as outfile:
        print('Writing fragments containing sulfones to sulfone.sdf')

    with open('sulfur_out.sdf', 'w') as outfile:
        print('Writing fragments containing sulfur that were removed to sulfone.sdf')

    with open('sulfur_cyclic.sdf', 'w') as outfile:
        print('Writing fragments containing cyclic sulfur to sulfur_cyclic.sdf')

    with open('check.sdf', 'w') as outfile:
        print('Writing fragments that will be kept after filter to check.sdf')

    if inchi_filter is not None:
        inchi_list = set()

        with open(inchi_filter) as infile:
            for line in infile:
                inchi_list.add(line.strip('\n'))

    if pains is True:

        from openeye import oechem
        from pymolgen.properties_pymolgen import gen_pains_database, pains_filter

        pains_database = gen_pains_database()

        with open('pains.sdf', 'w') as outfile:
            print('Writing fragments that do not pass pains test to pains.sdf')

    with open('thioether.sdf', 'w') as outfile:
        print('Writing fragments containing thioether to thioether.sdf')

    # loop through all fragments, if fragment is not filtered it will be added to filter_list at the end, otherwise it will be skiped by using a continue statement
    for i in range(len(fragment_database_mol)):

        mol = fragment_database_mol[i]

        mol_elements = list_elements(mol)

        elements.update(mol_elements)

        # only keep fragments with allowed elements
        allowed_elements = {'Cl', 'Br', 'I', 'O', 'C', 'F', 'N', 'S', 'H'}

        allowed = True
        
        for j in mol_elements:
            if j not in allowed_elements:
                allowed = False

        if allowed is False:
            continue

        # filter by structures in inchi_list, need to first hydrogenate to do comparison
        mol_h = mol.copy()
        mol_h.hydrogenate()

        inchi = molecule_to_inchi(mol_h)
        smi = molecule_to_smiles(mol_h)

        if inchi_filter is not None:
            if inchi in inchi_list:
                print('INCHI FILTER', inchi)
                continue

        # filter by pains
        if pains is True:

            oemol = oechem.OEGraphMol()
            oechem.OESmilesToMol(oemol, smi)

            oechem.OEAddExplicitHydrogens(oemol)

            if pains_filter(oemol, pains_database) is False:

                save_mol_to_sdf('pains.sdf', mol)

                continue

        rdmol = Chem.MolFromInchi(inchi)

        if rdmol is None:
            continue

        # filter out rings containing more than 8 atoms
        if rdmol is not None:
            ri = rdmol.GetRingInfo()
            largest_ring_size = max((len(r) for r in ri.AtomRings()), default=0)
            if largest_ring_size > 8:        

                save_mol_to_sdf('cycles.sdf', mol)

                continue

        # remove cages and rings by analysing number of aliphatic rings and heavy atoms in fragment
        n = rdMolDescriptors.CalcNumAliphaticRings(rdmol)

        if n == 3:

            if rdMolDescriptors.CalcNumHeavyAtoms(rdmol) < 12:

                save_mol_to_sdf('cages3.sdf', mol)

            else:

                save_mol_to_sdf('aliphatic3.sdf', mol)

                continue

        if n == 4:

            if rdMolDescriptors.CalcNumHeavyAtoms(rdmol) < 16:

                save_mol_to_sdf('cages4.sdf', mol)

                continue

            save_mol_to_sdf('aliphatic4.sdf', mol)

            continue

        if n > 4:

            save_mol_to_sdf('filter.sdf', mol)

            continue

        # remove thioether, keep sulfone, keep cyclic sulfur, remove non-cyclic sulfur
        if is_sulfur(mol):

            if is_thioether(mol):
                save_mol_to_sdf('thioether.sdf', mol)
                continue

            if is_sulfone(mol):
                save_mol_to_sdf('sulfone.sdf', mol)

            else:

                if is_cyclic_sulfur(mol):
                    save_mol_to_sdf('sulfur_cyclic.sdf', mol)

                else:

                    save_mol_to_sdf('sulfur_out.sdf', mol)

                    continue

        filter_list.append(i)

        save_mol_to_sdf('check.sdf', mol) 

    return filter_list


def save_mol_to_sdf(outfile_name, mol):

    with open(outfile_name, 'a') as outfile:
        lines = molecule_to_sdf(mol)

        for line in lines:
            outfile.write(line)

        outfile.write('$$$$\n')    


def is_sulfur(mol):

    for i in mol.graph.nodes:
        if mol.graph.nodes[i]["element"] == 'S':
            return True

    return False


def has_halogen(mol):

    for i in mol.graph.nodes:
        if mol.graph.nodes[i]["element"] in ['F', 'Cl', 'Br', 'I']:
            return True

    return False


def is_cyclic_sulfur(mol):

    for i in mol.graph.nodes:
        if mol.graph.nodes[i]["element"] == 'S':
            if mol.is_cyclic(i):
                return True
            else:
                for j in mol.graph[i]:
                    if mol.is_cyclic(j):
                        return True


    return False


def is_sulfone(mol):

    for i in mol.graph.nodes:
        if mol.graph.nodes[i]["element"] == 'S':
            o_count = 0
            for j in mol.graph[i]:
                if mol.graph.nodes[j]["element"] == 'O':
                    o_count += 1
            if o_count == 2:
                return True

    return False


def is_thioether(mol):

    for i in mol.graph.nodes:
        if mol.graph.nodes[i]["element"] == 'S':
            for j in mol.graph[i]:
                if mol.graph.nodes[j]["element"] == 'S':
                    return True

    return False


def copy_frequencies(fragment_database, bond_frequencies, frag_frequencies, fragment_a_i, fragment_b_sdf):

    #fragment_a_mol = molecule_from_sdf(fragment_a_sdf)
    fragment_b_mol = molecule_from_sdf(fragment_b_sdf).graph

    #fragment_a_i = find_fragment(fragment_a_mol, fragment_database)
    fragment_b_i = len(fragment_database)

    fragment_a_freq = frag_frequencies[fragment_a_i]

    frag_frequencies.append(fragment_a_freq)

    for key, val in bond_frequencies.copy().items():
        i = key[0]
        j = key[1]

        k = key[2]
        l = key[3]

        if fragment_a_i == i and fragment_a_i == j:
            bond_frequencies[i,fragment_b_i,k,l] = val
            bond_frequencies[fragment_b_i,fragment_b_i,k,l] = val

        elif fragment_a_i == i:
            bond_frequencies[j,fragment_b_i,l,k] = val

        elif fragment_a_i == j:
            bond_frequencies[i,fragment_b_i,k,l] = val            
 
    fragment_database.append(fragment_b_mol)

    return fragment_database, bond_frequencies, frag_frequencies


def exclude_aliphatic_halogen_bonds(fragment_database_mol, fragment_bond_frequencies):
    """
    Exclude aliphatic halogen bonds from database

    Parameters
    ----------
    fragment_database_mol : list of molecule objects
        List of fragments represented as molecule objects
    fragment_bond_frequencies : dict
        Dictionary mapping (i,j,k,l) bond to bond frequency


    Returns
    -------
    fragment_bond_frequencies : dict
        Dictionary mapping (i,j,k,l) bond to bond frequency with removed aliphatic halogen bonds    
    """

    with open('aromatic.sdf', 'w') as outfile:
        print('Writing fragments containing aromatic rings to aromatic.sdf')

    with open('halogen.sdf', 'w') as outfile:
        print('Writing fragments containing halogens to halogen.sdf')

    halogen_list = []

    aromatic_list = []

    for i in range(len(fragment_database_mol)):

        mol = fragment_database_mol[i]
        inchi = molecule_to_inchi(mol)
        rdmol = Chem.MolFromInchi(inchi) 

        if has_halogen(mol):
            halogen_list.append(i)
            save_mol_to_sdf('halogen.sdf', mol)
            continue

        if rdmol is None:
            continue

        aromatic_carbons = get_aromatic_carbons(rdmol)

        aromatic = False

        for j in mol.attach_points:
            if j in aromatic_carbons:
                aromatic = True

        if aromatic:
            aromatic_list.append(i)
            save_mol_to_sdf('aromatic.sdf', mol)

    remove_bond_frequencies_halogen(fragment_bond_frequencies, aromatic_list, halogen_list)

    return fragment_bond_frequencies


def get_aromatic_carbons(rdmol):

    aromatic_carbon = Chem.MolFromSmarts("c")

    aromatic_carbons = rdmol.GetSubstructMatches(aromatic_carbon)
    aromatic_carbons = [i[0] for i in aromatic_carbons]

    return aromatic_carbons


def loop(n, fragments_sdf_in, fragments_txt_in, frequencies_txt_in, frag_frequencies_txt_in, fragments_sdf_out, fragments_txt_out, frequencies_txt_out, frag_frequencies_txt_out, core=None, copy=None, limit=None, filter=False, filter_ids=None, first=None, inchi_filter=None, pains=False, sort=False, test=False):
    """
    Run main loop to combine fragment databases

    Parameters
    ----------
    n : int
        Number of fragment files to combine
    fragments_sdf_in : str
        Root name of input fragment database in SDF format
    fragments_txt_in : str
        Root name of input fragment database in TXT format
    frequencies_txt_in : str
        Root name of input bond frequencies file in TXT format
    frag_frequencies_txt_in : str
        Root name of input fragment frequencies file in TXT format
    fragments_sdf_out : str
        Root name of output fragment database in SDF format
    fragments_txt_out : str
        Root name of output fragment database in TXT format
    frequencies_txt_out : str
        Root name of output bond frequencies file in TXT format
    frag_frequencies_txt_out : str
        Root name of output fragment frequencies file in TXT format
    core : bool, optional
        Filter by frequency of core fragments (equivalent fragments of any protonation state)
    copy : pair of int, optional
        Pair of fragment ids to copy bond frequencies from one to the other, then exits
    limit : int, optional
        Limit for minimum fragment frequency to consider     
    filter : bool, optional
        Filter fragment database, then exits
    filter_ids : list of int, optional
        List of fragment ids to remove, then exits
    first : int, optional
        First file index to consider
    inchi_filter : str, optional
        Inchi file with inchis to filter out
    pains : bool, optional
        Filter fragments with pains
    sort : bool, optional
        Sort fragment database according to fragment frequencies
    test : bool, optional
        Test run, prints how many fragments would be removed based on limit

    Returns
    -------
    None (writes output files)
    """

    if first is not None:

        print('Loading %s' %first)

        fragment_database_mol = get_fragment_database('%s_%s.sdf' %(fragments_sdf_in, first))

        frequencies = get_bond_frequencies('%s_%s.txt' %(frequencies_txt_in, first) )

        frag_frequencies = get_frag_frequencies('%s_%s.txt' %(frag_frequencies_txt_in, first))

        frag_mapping = get_frag_mapping('%s_%s.txt' %(fragments_txt_in, first) )

    else:
        fragment_database_mol = get_fragment_database('%s.sdf' %(fragments_sdf_in) )

        frequencies = get_bond_frequencies('%s.txt' %(frequencies_txt_in) )

        frag_frequencies = get_frag_frequencies('%s.txt' %(frag_frequencies_txt_in))

        frag_mapping = get_frag_mapping('%s.txt' %(fragments_txt_in) )


    fragment_database = []

    for i in fragment_database_mol:
        fragment_database.append(i.graph)

    frequencies = update_bond_frequencies(frequencies, frag_mapping)

    if sort:
        fragment_database, frequencies, frag_frequencies, frag_mapping = sort_fragments(fragment_database, frequencies, frag_frequencies, frag_mapping)

        save_frequencies_txt(frequencies, frequencies_txt_out)

        save_fragments_sdf(fragment_database, fragments_sdf_out)

        save_frag_frequencies_txt(frag_frequencies, frag_frequencies_txt_out)

        save_fragments_txt(fragment_database, fragments_txt_out)

        sys.exit('Fragment data sorted')

    if copy is not None:

        fragment_a_i = int(copy[0])
        fragment_b_sdf = copy[1]

        fragment_database, frequencies, frag_frequencies = copy_frequencies(fragment_database, frequencies, frag_frequencies, fragment_a_i, fragment_b_sdf)

        save_frequencies_txt(frequencies, frequencies_txt_out)

        save_fragments_sdf(fragment_database, fragments_sdf_out)

        save_frag_frequencies_txt(frag_frequencies, frag_frequencies_txt_out)

        save_fragments_txt(fragment_database, fragments_txt_out)

        sys.exit('Fragment frequencies copied')

    if core is not None:

        core_frequencies = {}

        inchi_list = []

        for i in range(len(fragment_database_mol)):

            mol = fragment_database_mol[i]
            inchi = molecule_to_inchi(mol)
            inchi_list.append(inchi)

            if inchi in core_frequencies:
                core_frequencies[inchi] += frag_frequencies[i]
            else:
                core_frequencies[inchi] = frag_frequencies[i]

        core_frequencies = dict(sorted(core_frequencies.items(), key=lambda item: item[1], reverse=True))

        total_freq = []

        for i in range(len(fragment_database_mol)):
            total_freq.append(core_frequencies[inchi_list[i]])

        sorting = sorted(range(len(total_freq)), key=lambda k: total_freq[k], reverse=True)

        with open('cores.sdf', 'w') as f:
            for i in sorting:
                mol = fragment_database_mol[i]

                lines = molecule_to_sdf(mol)

                for line in lines:
                    f.write(line)
                f.write('$$$$\n')

        with open('cores.inchi', 'w') as f:
            for key, val in core_frequencies.items():
                f.write(f'{key} {val}\n')

        with open('cores_freq.txt', 'w') as f:
            for i in sorting:
                f.write(f'{inchi_list[i]} {total_freq[i]}\n')

        sys.exit('Core frequencies calculated')

    if filter:

        filter_list = filter_database(fragment_database_mol, inchi_filter, pains)

        with open('filter_in.sdf', 'w') as outfile:
            print('Writing fragments to keep to filter_in.sdf')

            for i in filter_list:

                mol = fragment_database_mol[i]

                lines = molecule_to_sdf(mol)

                for line in lines:
                    outfile.write(line)

                outfile.write('$$$$\n')               

        frequencies = remove_bond_frequencies(frequencies, filter_list)

        frequencies = exclude_aliphatic_halogen_bonds(fragment_database_mol, frequencies)

        print('Writing new bond frequencies to filter_bond_frequencies.txt')
        save_frequencies_txt(frequencies, 'filter_bond_frequencies.txt')       

        sys.exit('Fragment database filetered')

    if filter_ids is not None:

        remove_filter_ids(filter_ids, fragment_database_mol, frequencies)

        print('Fragment IDS filtered')
        sys.exit(0)

    if limit is not None:
        print('Before limit ', len(fragment_database))
        fragment_database, frequencies, frag_frequencies, frag_mapping = update_limit(limit, fragment_database, frequencies, frag_frequencies, frag_mapping)
        print('After limit ', len(fragment_database))

    if first is None:
        raise Exception ('Need to define first argument')

    for i in range(first+1, first+n):

        print('Loading %s' %i)

        fragments_sdf_2 = '%s_%s.sdf' %(fragments_sdf_in, i)
        fragments_txt_2 = '%s_%s.txt' %(fragments_txt_in, i) 
        frequencies_txt_2 = '%s_%s.txt' %(frequencies_txt_in, i) 
        frag_frequencies_txt_2 = '%s_%s.txt' %(frag_frequencies_txt_in, i)

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


def remove_filter_ids(filter_ids, fragment_database_mol, frequencies, folder='.'):
    """
    Remove fragments in filter_ids from fragment database frequencies

    Parameters
    ----------
    filter_ids : list of int
        list of fragment ids to remove
    fragment_database_mol : list of mol
        list of fragments as molecule objects
    frequencies : list of tuples
        list of bond frequencies as tuples
    folder : str, optional
        name of folder to save filter_ids.sdf with removed fragments

    Returns
    -------
    None (saves new frequency files with fragments removed)
    """

    # create filter_list of fragments to keep (those not in filter_ids)
    filter_list = [i for i in range(len(fragment_database_mol)) if i not in filter_ids]

    # save fragments being filtered
    with open(f'{folder}/filter_ids.sdf', 'w') as outfile:

        for i in filter_ids:

            mol = fragment_database_mol[i]

            lines = molecule_to_sdf(mol)

            for line in lines:
                outfile.write(line)

            outfile.write('$$$$\n')        

    frequencies = remove_bond_frequencies(frequencies, filter_list)

    save_frequencies_txt(frequencies, 'filter_id_bond_frequencies.txt')

    return frequencies


def update_limit(limit, fragment_database, bond_frequencies, frag_frequencies, frag_mapping):
    """
    Remove fragments that are below frag frequency limit
    Fragments are removed through the bond frequencies and frag_frequencies are updated accordinly, 
    but framents stay in the database

    Parameters
    ----------
    limit : int
        Minimum frag frequency to keep fragments
    fragment_database : list of molecules
        List of molecule objects representing fragments
    bond_frequencies : dict
        Dictionary mapping (i,j,k,l) bond to bond frequency
    frag_frequencies : list of int
        List of fragment frequencies
    frag_mapping : list of dict
        List of dictionaries mapping original atom numbers to new atom numbers in each fragment

    Returns
    -------
    new_fragment_database : list of molecules
        New list of molecule objects representing fragments
    new_bond_frequencies : dict
        New dictionary mapping (i,j,k,l) bond to bond frequency
    new_frag_frequencies : list of int
        New list of fragment frequencies
    new_frag_mapping : list of dict
        New list of dictionaries mapping original atom numbers to new atom numbers in each fragment
    """


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

    # update bond_frequencies: remove frequencies for removed fragments and update fragment numbers according to mapping

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


def sort_fragments(fragment_database, bond_frequencies, frag_frequencies, frag_mapping):
    """
    Sort fragment database according to fragment frequencies
    """

    frag_frequencies_np = np.array(frag_frequencies)
    sort_index = list(np.argsort(-1*frag_frequencies_np))

    mapping = {}

    for i in range(len(sort_index)):
        mapping[sort_index[i]] = i

    new_fragment_database = []
    new_bond_frequencies = {}
    new_frag_frequencies = []
    new_frag_mapping = []

    # add fragments to new_fragment_database according to sorted order, same for new_frag_mapping
    for i in range(len(sort_index)):
        new_fragment_database.append(fragment_database[sort_index[i]])
        new_frag_mapping.append(frag_mapping[sort_index[i]])

    new_frag_frequencies = sorted(frag_frequencies, reverse=True)

    # update bond frequencies for new fragment indeces, larger index should be j
    for key, val in bond_frequencies.items():

        i = mapping[key[0]]
        j = mapping[key[1]]
        k = key[2]
        l = key[3]

        if i <= j:
            new_bond_frequencies[i,j,k,l] = val            
        else:
            new_bond_frequencies[j,i,l,k] = val

    return new_fragment_database, new_bond_frequencies, new_frag_frequencies, new_frag_mapping


def renumber_frequencies(fragments_txt_in, frequencies_txt_in, frequencies_txt_out):

    frequencies = get_bond_frequencies(frequencies_txt_in)

    frag_mapping = get_frag_mapping(fragments_txt_in)

    frequencies = update_bond_frequencies(frequencies, frag_mapping)

    save_frequencies_txt(frequencies, frequencies_txt_out)


if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Combine fragment databases')
    
    # required arguments
    parser.add_argument('-n','--n_files', help='Number of fragment files to combine', type=int, required=True)
    parser.add_argument('-i','--in_sub', help='Input subscript', required=True)
    parser.add_argument('-o','--out_sub', help='Output subscript', required=True)

    # optional arguments
    parser.add_argument('-f','--first', help='First file index to consider', type=int, required=False)
    parser.add_argument('-l','--limit', help='Limit for minimum fragment frequency to consider', type=int, required=False)
    parser.add_argument('--core', type=int, help='Filter by frequency of core fragments (equivalent fragments of any protonation state)', required=False)
    parser.add_argument('--copy', nargs='+', help='Fragment_a and fragment_b to copy bond frequencies from a to b, then exits', required=False)
    parser.add_argument('--filter', action='store_true', help='Filter fragment database, then exits', required=False)
    parser.add_argument('--filter_ids', nargs='+', help='Space-separated list of fragment ids to remove, then exits', type=int, required=False)
    parser.add_argument('--inchi_filter', help='Inchi list to filter',required=False)
    parser.add_argument('--pains', action='store_true', help='Filter fragments with pains', required=False)
    parser.add_argument('--sort', action='store_true', help='Sort fragment data and exit', required=False)
    parser.add_argument('--test', action='store_true', help='Test run, prints how many fragments would be removed based on limit', required=False)

    args = parser.parse_args()

    n = args.n_files

    in_sub = args.in_sub
    out_sub = args.out_sub

    if args.first is not None:
        first = args.first * n

    fragments_sdf_in = 'fragments%s' %in_sub
    fragments_txt_in = 'fragments%s' %in_sub
    frequencies_txt_in = 'frequencies%s' %in_sub
    frag_frequencies_txt_in = 'frag_frequencies%s' %in_sub
    fragments_sdf_out = 'fragments%s.sdf' %out_sub
    fragments_txt_out = 'fragments%s.txt' %out_sub
    frequencies_txt_out = 'frequencies%s.txt' %out_sub
    frag_frequencies_txt_out = 'frag_frequencies%s.txt' %out_sub

    loop(n, fragments_sdf_in, fragments_txt_in, frequencies_txt_in, frag_frequencies_txt_in, fragments_sdf_out, fragments_txt_out, frequencies_txt_out, frag_frequencies_txt_out, core=args.core, copy=args.copy, filter=args.filter, filter_ids=args.filter_ids, first=args.first, limit=args.limit, inchi_filter=args.inchi_filter, pains=args.pains, sort=args.sort, test=args.test)

    print('Normal termination')