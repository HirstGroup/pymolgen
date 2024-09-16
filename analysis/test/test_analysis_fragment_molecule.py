import subprocess

from pymolgen.analysis.analysis_fragment_molecule import *


def test1():

    subprocess.run('python ../analysis_fragment_molecule.py -i input/inchi1.inchi -o output/inchi1_analysis.txt -a ../../datasets/fragments/fragments_30_50k_co_10_l5_5_sorted_filter_copy.sdf -p input/phenylisoxazole.sdf -r 20 21 -rf ../../datasets/fragments/fragments_30_50k_co_10_l5_5_sorted_filter_copy.txt', check=True, shell=True)


def test2():

    subprocess.run('python ../analysis_fragment_molecule.py -i input/inchi10.inchi -o output/inchi10_analysis.txt -a ../../datasets/fragments/fragments_30_50k_co_10_l5_5_sorted_filter_copy.sdf -p input/phenylisoxazole.sdf -r 20 21 -rf ../../datasets/fragments/fragments_30_50k_co_10_l5_5_sorted_filter_copy.txt', check=True, shell=True)
test2()


def test_get_fragment_index():

    fragment_database = get_fragment_database('../../datasets/database1000/fragments10.sdf')

    mol = molecule_from_sdf('../../datasets/database1000/ch4.sdf')

    mol = mol.remove_atom(4)

    index, newmap = get_fragment_index(mol.graph, fragment_database)

    print(index, newmap)

    mol = molecule_from_sdf('input/amide.sdf')

    mol = mol.remove_atom(5)
    mol = mol.remove_atom(3)

    index, newmap = get_fragment_index(mol.graph, fragment_database)

    assert newmap == {4: 0, 2: 1, 0: 2, 1: 3}
    assert index == 2

    print(index, newmap)


def test_analyse_molecule_1():

    fragment_database = get_fragment_database('../../datasets/database1000/fragments10.sdf')
    fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    mol = molecule_from_smiles('CC')
    inchi = molecule_to_inchi(mol)
    
    string_representation = analyse_molecule(inchi, fragment_database, fragment_database_graph)

    assert string_representation == '0-0:(0, 1, 0, 0):1.0'


def test_analyse_molecule_2():

    fragment_database = get_fragment_database('../../datasets/database1000/fragments10.sdf')
    fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    mol = molecule_from_smiles('CCC')
    inchi = molecule_to_inchi(mol)
    
    string_representation = analyse_molecule(inchi, fragment_database, fragment_database_graph)

    assert string_representation == '0-3-0:(0, 1, 0, 2);(1, 2, 2, 0):1.0'


def test_analyse_molecule_3():

    fragment_database = get_fragment_database('../../datasets/database1000/fragments10.sdf')
    fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    mol = molecule_from_smiles('CCC')
    inchi = molecule_to_inchi(mol)
    
    string_representation = analyse_molecule(inchi, fragment_database, fragment_database_graph)

    assert string_representation == '0-3-0:(0, 1, 0, 2);(1, 2, 2, 0):1.0'

    f = generate_fragment_molecule_from_string(string_representation, fragment_database_graph)

    mol = convert_fragment_molecule_to_mol(f, fragment_database)

    inchi2 = molecule_to_inchi(mol)

    print(inchi, inchi2)

    assert inchi == inchi2


def test_convert_parent():

    fragment_database = get_fragment_database('../../datasets/database1000/fragments10.sdf')
    fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    parent = FragmentMolecule()
    parent.add_fragment(0, [0], {0:0})
    parent.add_fragment(0, [0], {0:0})
    parent.add_bond(0, 1, 0, 0)

    print(parent)

    f = FragmentMolecule()
    f.add_fragment(0, [0], {0:0})
    f.add_fragment(0, [0,0], {0:0})
    f.add_fragment(0, [0], {0:0})
    f.add_bond(0, 1, 0, 0)
    f.add_bond(1, 2, 0, 0)
    print(f)    

    convert_parent(f, parent)
test_convert_parent()