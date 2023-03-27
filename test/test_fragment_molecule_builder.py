import os
import sys

from pymolgen.fragment_molecule_builder import *
from pymolgen.fragment_graph import convert_fragment_database_to_graph

def test_fragment_molecule_builder():

    f = FragmentMolecule()

    f.add_fragment(10, [0,1]) # 0
    f.add_fragment(20, [2,2]) # 1
    f.add_fragment(20, [2,2]) # 2
    f.add_fragment(30, [4])   # 3
    f.add_fragment(40, [5])   # 4

    f.add_bond(0, 1, 0, 2)
    assert f.list_free_valence_points() == [[1], [2], [2, 2], [4], [5]]

    f.add_bond(1, 2, 2, 2)
    assert f.list_free_valence_points() == [[1], [], [2], [4], [5]]

    f.add_bond(2, 3, 2, 4)
    assert f.list_free_valence_points() == [[1], [], [], [], [5]]

    f.add_bond(0, 4, 1, 5)
    assert f.list_free_valence_points() == [[], [], [], [], []]

    assert f.list_frag_id() == [10, 20, 20, 30, 40]

    assert f.list_bonds() == [(0, 1, 0, 2), (1, 2, 2, 2), (2, 3, 2, 4), (0, 4, 1, 5)]

def test_convert_fragment_molecule_to_mol():

    f = FragmentMolecule()

    f.add_fragment(0, [0])
    f.add_fragment(0, [0])

    f.add_bond(0, 1, 0, 0)

    fragment_database = get_fragment_database('../datasets/database1000/fragments1.sdf')

    mol = convert_fragment_molecule_to_mol(f, fragment_database)

    assert molecule_to_inchi(mol) == 'InChI=1S/C2H6/c1-2/h1-2H3'


def test_extend_molecule_list():

    bond_frequencies = get_bond_frequencies('../datasets/database1000/frequencies1.txt')
    bond_frequencies = bond_frequencies_to_np(bond_frequencies)

    fragment_database = get_fragment_database('../datasets/database1000/fragments1.sdf')
    fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    bond_frequencies = convert_bond_freq_np_to_dict(fragment_database_graph, bond_frequencies)

    ch3 = FragmentMolecule()

    ch3.add_fragment(0, [0])

    output_mol_list = extend_molecule_list([ch3], bond_frequencies, fragment_database_graph)

    answers = ['0-1']

    for idx, x in enumerate(output_mol_list):
        assert str(x) == answers[idx]


def test_extend_molecule_list_2():

    bond_frequencies = get_bond_frequencies('../datasets/database1000/frequencies1.txt')
    bond_frequencies = bond_frequencies_to_np(bond_frequencies)

    fragment_database = get_fragment_database('../datasets/database1000/fragments1.sdf')
    fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    bond_frequencies = convert_bond_freq_np_to_dict(fragment_database_graph, bond_frequencies)

    amide = FragmentMolecule()

    amide.add_fragment(2, [1, 2])

    output_mol_list = extend_molecule_list([amide], bond_frequencies, fragment_database_graph)
    
    answers = ['2-1', '2-3']
    for idx, x in enumerate(output_mol_list):
        assert str(x) == answers[idx]

def test_extend_molecule_list_all():

    bond_frequencies = get_bond_frequencies('../datasets/database1000/frequencies1.txt')
    bond_frequencies = bond_frequencies_to_np(bond_frequencies)

    fragment_database = get_fragment_database('../datasets/database1000/fragments1.sdf')
    fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    bond_frequencies = convert_bond_freq_np_to_dict(fragment_database_graph, bond_frequencies)
    
    all_output_mol_list = []

    for i in range(len(fragment_database)):

        mol = fragment_database[i]

        mol2 = FragmentMolecule()

        mol2.add_fragment(i, mol.free_valence_list)

        output_mol_list = extend_molecule_list([mol2], bond_frequencies, fragment_database_graph)

        all_output_mol_list.extend(output_mol_list)

    answers = ['InChI=1S/C4H5NO/c1-4-2-3-5-6-4/h2-3H,1H3', 'InChI=1S/C4H5NO/c1-4-2-3-5-6-4/h2-3H,1H3', 'InChI=1S/C4H4N2O2/c7-3-5-4-1-2-8-6-4/h1-3H,(H,5,6,7)', 'InChI=1S/C4H4N2O2/c7-3-5-4-1-2-8-6-4/h1-3H,(H,5,6,7)', 'InChI=1S/C2H5NO/c1-2(3)4/h1H3,(H2,3,4)', 'InChI=1S/C2H5NO/c1-2(3)4/h1H3,(H2,3,4)', 'InChI=1S/CH4O2S/c1-4(2)3/h4H,1H3', 'InChI=1S/C9H9N/c1-10-7-6-8-4-2-3-5-9(8)10/h2-7H,1H3', 'InChI=1S/C7H8/c1-7-5-3-2-4-6-7/h2-6H,1H3', 'InChI=1S/C2H5NO/c1-2(3)4/h1H3,(H2,3,4)', 'InChI=1S/CH4O2S/c1-4(2)3/h4H,1H3', 'InChI=1S/C9H9N/c1-10-7-6-8-4-2-3-5-9(8)10/h2-7H,1H3', 'InChI=1S/C7H8/c1-7-5-3-2-4-6-7/h2-6H,1H3', 'InChI=1S/CH4O2S/c1-4(2)3/h4H,1H3', 'InChI=1S/C8H7NO2S/c10-12(11)8-5-9-7-4-2-1-3-6(7)8/h1-5,9,12H', 'InChI=1S/CH4O2S/c1-4(2)3/h4H,1H3', 'InChI=1S/C8H7NO2S/c10-12(11)8-5-9-7-4-2-1-3-6(7)8/h1-5,9,12H', 'InChI=1S/C8H7NO2S/c10-12(11)8-5-9-7-4-2-1-3-6(7)8/h1-5,9,12H', 'InChI=1S/C9H9N/c1-10-7-6-8-4-2-3-5-9(8)10/h2-7H,1H3', 'InChI=1S/C8H6FN/c9-7-2-1-6-3-4-10-8(6)5-7/h1-5,10H', 'InChI=1S/C8H6FN/c9-7-2-1-3-8-6(7)4-5-10-8/h1-5,10H', 'InChI=1S/C7H8/c1-7-5-3-2-4-6-7/h2-6H,1H3', 'InChI=1S/C9H8N2/c1-2-4-8(5-3-1)9-6-10-11-7-9/h1-7H,(H,10,11)', 'InChI=1S/C7H8/c1-7-5-3-2-4-6-7/h2-6H,1H3', 'InChI=1S/C9H8N2/c1-2-4-8(5-3-1)9-6-10-11-7-9/h1-7H,(H,10,11)', 'InChI=1S/C9H8N2/c1-2-4-8(5-3-1)9-6-10-11-7-9/h1-7H,(H,10,11)', 'InChI=1S/C8H6FN/c9-7-2-1-6-3-4-10-8(6)5-7/h1-5,10H', 'InChI=1S/C8H6FN/c9-7-2-1-3-8-6(7)4-5-10-8/h1-5,10H']

    for idx, x in enumerate(all_output_mol_list):
        mol3 = convert_fragment_molecule_to_mol(x, fragment_database)
        inchi = molecule_to_inchi(mol3)
        assert inchi == answers[idx]

def test_attach_points():

    fragment_database = get_fragment_database('../datasets/database1000/fragments1.sdf')

    out = []

    for i in range(len(fragment_database)):
        out.append(fragment_database[i].attach_points)

    assert out == [[0], [0, 2], [1, 2], [2], [0], [1, 3, 9, 11], [1, 5], [3], [0]]

def test_free_valence_original():

    fragment_database = get_fragment_database('../datasets/database1000/fragments1.sdf')

    out = []

    for i in range(len(fragment_database)):
        out.append(fragment_database[i].free_valence_original)

    assert out == [[0], [0, 2], [1, 2], [2, 2], [0, 0], [1, 3, 9, 11], [1, 5], [3], [0]]

def test_extend_molecule_list_depth():

    bond_frequencies = get_bond_frequencies('../datasets/database1000/frequencies1.txt')
    bond_frequencies = bond_frequencies_to_np(bond_frequencies)

    fragment_database = get_fragment_database('../datasets/database1000/fragments1.sdf')
    fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    bond_frequencies = convert_bond_freq_np_to_dict(fragment_database_graph, bond_frequencies)

    ch3 = FragmentMolecule()

    ch3.add_fragment(0, [0])

    output_mol_list = extend_molecule_list_depth([ch3], bond_frequencies, fragment_database_graph, depth=1, unique=False)

    for x in output_mol_list:
        mol = convert_fragment_molecule_to_mol(x, fragment_database)
        inchi = molecule_to_inchi(mol)
        assert str(x) == '0-1'
        assert inchi == 'InChI=1S/C4H5NO/c1-4-2-3-5-6-4/h2-3H,1H3'

    output_mol_list = extend_molecule_list_depth([ch3], bond_frequencies, fragment_database_graph, depth=2, unique=False)

    for x in output_mol_list:
        mol = convert_fragment_molecule_to_mol(x, fragment_database)
        inchi = molecule_to_inchi(mol)
        assert str(x) == '0-1-2'
        assert inchi == 'InChI=1S/C5H6N2O2/c1-4-2-5(6-3-8)7-9-4/h2-3H,1H3,(H,6,7,8)'

    output_mol_list = extend_molecule_list_depth([ch3], bond_frequencies, fragment_database_graph, depth=3, unique=False)

    for x in output_mol_list:
        mol = convert_fragment_molecule_to_mol(x, fragment_database)
        inchi = molecule_to_inchi(mol)
        assert str(x) == '0-1-2-3'
        assert inchi == 'InChI=1S/C6H8N2O2/c1-4-3-6(8-10-4)7-5(2)9/h3H,1-2H3,(H,7,8,9)'

    output_mol_list = extend_molecule_list_depth([ch3], bond_frequencies, fragment_database_graph, depth=4, unique=False)

    answers = ['InChI=1S/C7H9N3O3/c1-4-2-6(10-13-4)9-7(12)3-5(8)11/h2H,3H2,1H3,(H2,8,11)(H,9,10,12)', 'InChI=1S/C6H8N2O4S/c1-4-2-5(8-12-4)7-6(9)3-13(10)11/h2,13H,3H2,1H3,(H,7,8,9)', 'InChI=1S/C14H13N3O2/c1-10-8-13(16-19-10)15-14(18)9-17-7-6-11-4-2-3-5-12(11)17/h2-8H,9H2,1H3,(H,15,16,18)', 'InChI=1S/C12H12N2O2/c1-9-7-11(14-16-9)13-12(15)8-10-5-3-2-4-6-10/h2-7H,8H2,1H3,(H,13,14,15)']

    for idx, x in enumerate(output_mol_list):
        mol = convert_fragment_molecule_to_mol(x, fragment_database)
        inchi = molecule_to_inchi(mol)
        print(inchi)
        assert inchi == answers[idx]

    output_mol_list = extend_molecule_list_depth([ch3], bond_frequencies, fragment_database_graph, depth=5, unique=False)

    answers = ['InChI=1S/C10H10N4O4/c1-6-4-8(14-18-6)12-10(16)5-9(15)11-7-2-3-17-13-7/h2-4H,5H2,1H3,(H,11,13,15)(H,12,14,16)', 'InChI=1S/C7H10N2O4S/c1-5-3-6(9-13-5)8-7(10)4-14(2,11)12/h3H,4H2,1-2H3,(H,8,9,10)', 'InChI=1S/C14H13N3O4S/c1-9-6-13(17-21-9)16-14(18)8-22(19,20)12-7-15-11-5-3-2-4-10(11)12/h2-7,15H,8H2,1H3,(H,16,17,18)', 'InChI=1S/C14H13N3O4S/c1-9-6-13(16-21-9)15-14(18)8-17-7-12(22(19)20)10-4-2-3-5-11(10)17/h2-7,22H,8H2,1H3,(H,15,16,18)', 'InChI=1S/C14H12FN3O2/c1-9-6-13(17-20-9)16-14(19)8-18-5-4-10-2-3-11(15)7-12(10)18/h2-7H,8H2,1H3,(H,16,17,19)', 'InChI=1S/C14H12FN3O2/c1-9-7-13(17-20-9)16-14(19)8-18-6-5-10-11(15)3-2-4-12(10)18/h2-7H,8H2,1H3,(H,16,17,19)', 'InChI=1S/C13H14N2O2/c1-9-4-3-5-11(6-9)8-13(16)14-12-7-10(2)17-15-12/h3-7H,8H2,1-2H3,(H,14,15,16)', 'InChI=1S/C15H14N4O2/c1-10-5-14(19-21-10)18-15(20)7-11-3-2-4-12(6-11)13-8-16-17-9-13/h2-6,8-9H,7H2,1H3,(H,16,17)(H,18,19,20)']

    for idx, x in enumerate(output_mol_list):
        mol = convert_fragment_molecule_to_mol(x, fragment_database)
        inchi = molecule_to_inchi(mol)
        assert inchi == answers[idx]


def test_extend_molecule_list_depth_simple():

    bond_frequencies = get_bond_frequencies('../datasets/simple/frequencies_simple.txt')
    bond_frequencies = bond_frequencies_to_np(bond_frequencies)

    fragment_database = get_fragment_database('../datasets/simple/fragments_simple.sdf')
    fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    bond_frequencies = convert_bond_freq_np_to_dict(fragment_database_graph, bond_frequencies)

    ch3 = FragmentMolecule()

    ch3.add_fragment(0, [0])

    output_mol_list = extend_molecule_list_depth([ch3], bond_frequencies, fragment_database_graph, depth=1, unique=False)

    answers = ['0-0', '0-1']

    for j in range(len(output_mol_list)):
        assert str(output_mol_list[j]) == answers[j]

def test_extend_molecule_list_depth_rzt():

    bond_frequencies = get_bond_frequencies('../datasets/simple/frequencies_rzt.txt')
    print(bond_frequencies)
    bond_frequencies = bond_frequencies_to_np(bond_frequencies)

    fragment_database = get_fragment_database('../datasets/simple/fragments_rzt.sdf')
    fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    bond_frequencies = convert_bond_freq_np_to_dict(fragment_database_graph, bond_frequencies, sort_dict=False)

    ch3 = FragmentMolecule()

    ch3.add_fragment(0, [0])

    output_mol_list = extend_molecule_list_depth([ch3], bond_frequencies, fragment_database_graph, depth=1, unique=False)

    for j in output_mol_list:
        mol = convert_fragment_molecule_to_mol(j, fragment_database)
        inchi = molecule_to_inchi(mol)
        assert inchi == 'InChI=1S/C7H8/c1-7-5-3-2-4-6-7/h2-6H,1H3'

    output_mol_list = extend_molecule_list_depth([ch3], bond_frequencies, fragment_database_graph, depth=2, unique=False)

    answers = ['InChI=1S/C8H10/c1-7-3-5-8(2)6-4-7/h3-6H,1-2H3', 'InChI=1S/C13H12/c1-11-7-9-13(10-8-11)12-5-3-2-4-6-12/h2-10H,1H3']

    for idx, x in enumerate(output_mol_list):
        mol = convert_fragment_molecule_to_mol(x, fragment_database)
        inchi = molecule_to_inchi(mol)
        print(inchi)
        #assert inchi == answers[idx]

    output_mol_list = extend_molecule_list_depth([ch3], bond_frequencies, fragment_database_graph, depth=3, unique=False)

    answers = ['InChI=1S/C14H14/c1-11-3-7-13(8-4-11)14-9-5-12(2)6-10-14/h3-10H,1-2H3','InChI=1S/C19H16/c1-15-7-9-17(10-8-15)19-13-11-18(12-14-19)16-5-3-2-4-6-16/h2-14H,1H3']

    for idx, x in enumerate(output_mol_list):
        mol = convert_fragment_molecule_to_mol(x, fragment_database)
        inchi = molecule_to_inchi(mol)
        assert inchi == answers[idx]

def test_convert_fragment_database_to_graph():

    fragment_database = get_fragment_database('../datasets/database1000/fragments1.sdf')

    fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    attachment_point_list = []
    canonical_mapping_list = []
    canonical_mapping_list_get = []

    for i in range(len(fragment_database_graph.fragments)):
        attachment_point_list.append(fragment_database_graph.fragments[i].attachment_points)
        canonical_mapping_list.append(fragment_database_graph.fragments[i].set_canonical_mapping(fragment_database))
        canonical_mapping_list_get.append(fragment_database_graph.fragments[i].get_canonical_mapping())

    assert attachment_point_list == [[0], [0, 2], [1, 2], [2, 2], [0, 0], [1, 3, 9, 11], [1, 5], [3], [0]]
    assert canonical_mapping_list == [{0: 0, 1: 1, 2: 1, 3: 1}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5}, {0: 0, 1: 1, 2: 2, 3: 3}, {0: 0, 2: 2, 1: 0}, {0: 0, 1: 1, 2: 1}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 7: 7, 8: 8, 5: 5, 9: 9, 10: 10, 6: 6, 11: 11}, {0: 0, 2: 2, 1: 1, 3: 3, 4: 2, 5: 1, 7: 7, 8: 0, 9: 9, 6: 6}, {0: 0, 4: 4, 3: 3, 5: 5, 6: 6, 1: 1, 7: 7, 2: 2}, {0: 0}]
    assert canonical_mapping_list_get == [{0: 0, 1: 1, 2: 1, 3: 1}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5}, {0: 0, 1: 1, 2: 2, 3: 3}, {0: 0, 2: 2, 1: 0}, {0: 0, 1: 1, 2: 1}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 7: 7, 8: 8, 5: 5, 9: 9, 10: 10, 6: 6, 11: 11}, {0: 0, 2: 2, 1: 1, 3: 3, 4: 2, 5: 1, 7: 7, 8: 0, 9: 9, 6: 6}, {0: 0, 4: 4, 3: 3, 5: 5, 6: 6, 1: 1, 7: 7, 2: 2}, {0: 0}]

    for i in range(len(fragment_database_graph)):
        print(i, fragment_database_graph.fragments[i].attachment_points, fragment_database_graph.fragments[i].get_canonical_mapping())

def test_write_fragment_database_graph():

    fragment_database = get_fragment_database('../datasets/database1000/fragments1.sdf')

    fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    write_fragment_database_graph(fragment_database_graph, '../datasets/database1000/fragment_database1.txt')

def test_read_fragment_database_graph():

    fragment_database_graph = read_fragment_database_graph('../datasets/database1000/fragment_database1.txt')

    fragment_database = get_fragment_database('../datasets/database1000/fragments1.sdf')

    fragment_database_graph2 = convert_fragment_database_to_graph(fragment_database)

    for i in range(len(fragment_database_graph)):

        assert fragment_database_graph.fragments[i].get_canonical_mapping() == fragment_database_graph2.fragments[i].get_canonical_mapping()
        assert fragment_database_graph.fragments[i].attachment_points == fragment_database_graph2.fragments[i].attachment_points
        assert fragment_database_graph.fragments[i].get_attribute('frag_id') == fragment_database_graph2.fragments[i].get_attribute('frag_id')

def test_extend_molecule_list_model1():

    f = FragmentMolecule()

    f.add_fragment(0, [1,2,3])

    bond_frequencies = {(0,0,1,1):1}

    bond_frequencies = bond_frequencies_to_np(bond_frequencies)

    fragment_database_graph = FragmentGraph()

    fragment_database_graph.add_fragment(0, [1,2,3])

    fragment_database_graph.fragments[0].manual_canonical_mapping({1:1, 2:2, 3:3})

    bond_frequencies = convert_bond_freq_np_to_dict(fragment_database_graph, bond_frequencies)

    output_mol_list = extend_molecule_list_depth([f], bond_frequencies, fragment_database_graph, depth=1, unique=False)

    assert len(output_mol_list) == 1
    assert str(output_mol_list[0]) == '0-0'
    assert output_mol_list[0].list_bonds() == [(0, 1, 1, 1)]

    output_mol_list = extend_molecule_list_depth([f], bond_frequencies, fragment_database_graph, depth=2, unique=False)

    for i in output_mol_list:
        print(i, i.list_bonds())

    assert len(output_mol_list) == 0

def test_extend_molecule_list_model2():

    f = FragmentMolecule()

    f.add_fragment(0, [1,2,3])

    bond_frequencies = {(0,0,2,2):1}

    bond_frequencies = bond_frequencies_to_np(bond_frequencies)

    fragment_database_graph = FragmentGraph()

    fragment_database_graph.add_fragment(0, [1,2,3])

    fragment_database_graph.fragments[0].manual_canonical_mapping({1:1, 2:2, 3:3})

    bond_frequencies = convert_bond_freq_np_to_dict(fragment_database_graph, bond_frequencies)

    output_mol_list = extend_molecule_list_depth([f], bond_frequencies, fragment_database_graph, depth=1, unique=False)

    for i in output_mol_list:
        print(i, i.list_bonds())    

    assert len(output_mol_list) == 1
    assert str(output_mol_list[0]) == '0-0'
    assert output_mol_list[0].list_bonds() == [(0, 1, 2, 2)]

    output_mol_list = extend_molecule_list_depth([f], bond_frequencies, fragment_database_graph, depth=2, unique=False)

    for i in output_mol_list:
        print(i, i.list_bonds())    

    assert len(output_mol_list) == 0

def test_extend_molecule_list_model3():

    f = FragmentMolecule()

    f.add_fragment(0, [1,2,3])

    bond_frequencies = {(0,0,1,2):1}

    bond_frequencies = bond_frequencies_to_np(bond_frequencies)

    fragment_database_graph = FragmentGraph()

    fragment_database_graph.add_fragment(0, [1,2,3])

    fragment_database_graph.fragments[0].manual_canonical_mapping({1:1, 2:2, 3:3})

    bond_frequencies = convert_bond_freq_np_to_dict(fragment_database_graph, bond_frequencies)

    output_mol_list = extend_molecule_list_depth([f], bond_frequencies, fragment_database_graph, depth=1, unique=False)

    assert len(output_mol_list) == 2

    mols = ['0-0', '0-0']
    bonds = [[(0, 1, 1, 2)], [(0, 1, 2, 1)]]

    for idx, x in enumerate(output_mol_list):
        print(x, x.list_bonds())
        assert str(x) == mols[idx]
        assert x.list_bonds() == bonds[idx]

    output_mol_list = extend_molecule_list_depth([f], bond_frequencies, fragment_database_graph, depth=2, unique=False)

    mols = ['0-0-0','0-0-0','0-0-0','0-0-0']
    bonds = [[(0, 1, 1, 2), (0, 2, 2, 1)],[(0, 1, 1, 2), (1, 2, 1, 2)],[(0, 1, 2, 1), (0, 2, 1, 2)],[(0, 1, 2, 1), (1, 2, 2, 1)]]

    for idx, x in enumerate(output_mol_list):
        print(x, x.list_bonds())
        assert str(x) == mols[idx]
        assert x.list_bonds() == bonds[idx]

def test_extend_molecule_list_model4():

    f = FragmentMolecule()

    f.add_fragment(0, [1,2,3])

    bond_frequencies = {(0,0,1,1):1,(0,0,1,2):1}

    bond_frequencies = bond_frequencies_to_np(bond_frequencies)

    fragment_database_graph = FragmentGraph()

    fragment_database_graph.add_fragment(0, [1,2,3])

    fragment_database_graph.fragments[0].manual_canonical_mapping({1:1, 2:2, 3:3})

    bond_frequencies = convert_bond_freq_np_to_dict(fragment_database_graph, bond_frequencies)

    output_mol_list = extend_molecule_list_depth([f], bond_frequencies, fragment_database_graph, depth=1, unique=False)

    #assert len(output_mol_list) == 2

    mols = ['0-0', '0-0','0-0']
    bonds = [[(0, 1, 1, 1)], [(0, 1, 1, 2)], [(0, 1, 2, 1)]]

    for idx, x in enumerate(output_mol_list):
        print(x, x.list_bonds())
        assert str(x) == mols[idx]
        assert x.list_bonds() == bonds[idx]

    output_mol_list = extend_molecule_list_depth([f], bond_frequencies, fragment_database_graph, depth=2, unique=False)

    bonds = [[(0, 1, 1, 1), (0, 2, 2, 1)], [(0, 1, 1, 1), (1, 2, 2, 1)], [(0, 1, 1, 2), (0, 2, 2, 1)], [(0, 1, 1, 2), (1, 2, 1, 1)], [(0, 1, 1, 2), (1, 2, 1, 2)], [(0, 1, 2, 1), (0, 2, 1, 1)], [(0, 1, 2, 1), (0, 2, 1, 2)], [(0, 1, 2, 1), (1, 2, 2, 1)]]

    for idx, x in enumerate(output_mol_list):
        print(x, x.list_bonds())
        assert x.list_bonds() == bonds[idx]


def test_extend_molecule_list_model5():

    f = FragmentMolecule()

    f.add_fragment(0, [1,2,3])

    bond_frequencies = {(0,0,1,2):1}

    bond_frequencies = bond_frequencies_to_np(bond_frequencies)

    fragment_database_graph = FragmentGraph()

    fragment_database_graph.add_fragment(0, [1,2,3])

    fragment_database_graph.fragments[0].manual_canonical_mapping({1:1, 2:1, 3:3})

    bond_frequencies = convert_bond_freq_np_to_dict(fragment_database_graph, bond_frequencies)

    output_mol_list = extend_molecule_list_depth([f], bond_frequencies, fragment_database_graph, depth=1, unique=False)

    mols = ['0-0', '0-0']
    bonds = [[(0, 1, 1, 2)], [(0, 1, 2, 2)]]

    for idx, x in enumerate(output_mol_list):
        print(x, x.list_bonds())
        assert str(x) == mols[idx]
        assert x.list_bonds() == bonds[idx]

    output_mol_list = extend_molecule_list_depth([f], bond_frequencies, fragment_database_graph, depth=2, unique=False)

    mols = ['0-0-0', '0-0-0', '0-0-0', '0-0-0']
    bonds = [[(0, 1, 1, 2), (0, 2, 2, 2)], [(0, 1, 1, 2), (1, 2, 1, 2)], [(0, 1, 2, 2), (0, 2, 1, 2)], [(0, 1, 2, 2), (1, 2, 1, 2)]]

    for idx, x in enumerate(output_mol_list):
        print(x, x.list_bonds())
        assert str(x) == mols[idx]
        assert x.list_bonds() == bonds[idx]

def test_extend_molecule_list_model6():

    f = FragmentMolecule()

    f.add_fragment(0, [1,2,3])

    bond_frequencies = {(0,0,1,2):1}

    bond_frequencies = bond_frequencies_to_np(bond_frequencies)

    fragment_database_graph = FragmentGraph()

    fragment_database_graph.add_fragment(0, [1,2,3])

    fragment_database_graph.fragments[0].manual_canonical_mapping({1:1, 2:1, 3:1})

    bond_frequencies = convert_bond_freq_np_to_dict(fragment_database_graph, bond_frequencies)

    output_mol_list = extend_molecule_list_depth([f], bond_frequencies, fragment_database_graph, depth=1, unique=False)

    mols = ['0-0', '0-0', '0-0']
    bonds = [[(0, 1, 1, 2)], [(0, 1, 2, 2)], [(0, 1, 3, 2)]]

    for idx, x in enumerate(output_mol_list):
        print(x, x.list_bonds())
        assert str(x) == mols[idx]
        assert x.list_bonds() == bonds[idx]

    output_mol_list = extend_molecule_list_depth([f], bond_frequencies, fragment_database_graph, depth=2, unique=False)

    mols = ['0-0-0'] * 12
    bonds = [[(0, 1, 1, 2), (0, 2, 2, 2)], [(0, 1, 1, 2), (0, 2, 3, 2)], [(0, 1, 1, 2), (1, 2, 1, 2)], [(0, 1, 1, 2), (1, 2, 3, 2)], [(0, 1, 2, 2), (0, 2, 1, 2)], [(0, 1, 2, 2), (0, 2, 3, 2)], [(0, 1, 2, 2), (1, 2, 1, 2)], [(0, 1, 2, 2), (1, 2, 3, 2)], [(0, 1, 3, 2), (0, 2, 1, 2)], [(0, 1, 3, 2), (0, 2, 2, 2)], [(0, 1, 3, 2), (1, 2, 1, 2)], [(0, 1, 3, 2), (1, 2, 3, 2)]]


    for idx, x in enumerate(output_mol_list):
        print(idx, x, x.list_bonds())
        assert str(x) == mols[idx]
        assert x.list_bonds() == bonds[idx]


def test_extend_molecule_list_database11_20():

    bond_frequencies = get_bond_frequencies('../datasets/database1000/frequencies_11-20.txt')
    bond_frequencies = bond_frequencies_to_np(bond_frequencies)

    fragment_database = get_fragment_database('../datasets/database1000/fragments_11-20.sdf')

    fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    bond_frequencies = convert_bond_freq_np_to_dict(fragment_database_graph, bond_frequencies, sort_dict=False)

    parent = FragmentMolecule()

    parent.add_fragment(0, [0])

    output_mol_list = extend_molecule_list_depth([parent], bond_frequencies, fragment_database_graph, 2, unique=False)

    inchi_list = ['InChI=1S/C2H6O/c1-3-2/h1-2H3', 'InChI=1S/C10H9NO/c1-12-9-4-5-10-8(7-9)3-2-6-11-10/h2-7H,1H3', 'InChI=1S/C14H11NO/c1-16-12-6-7-14-11(9-12)8-10-4-2-3-5-13(10)15-14/h2-9H,1H3', 'InChI=1S/C7H8O/c1-8-7-5-3-2-4-6-7/h2-6H,1H3', 'InChI=1S/C2H6O/c1-3-2/h1-2H3', 'InChI=1S/C2H7N/c1-3-2/h3H,1-2H3', 'InChI=1S/C2H7N/c1-3-2/h3H,1-2H3', 'InChI=1S/C7H9N/c1-8-7-5-3-2-4-6-7/h2-6,8H,1H3', 'InChI=1S/C2H7N/c1-3-2/h3H,1-2H3', 'InChI=1S/C2H7N/c1-3-2/h3H,1-2H3', 'InChI=1S/C7H9N/c1-8-7-5-3-2-4-6-7/h2-6,8H,1H3', 'InChI=1S/C4H7N3/c1-3-5-4(2)7-6-3/h1-2H3,(H,5,6,7)', 'InChI=1S/C4H7N3/c1-3-5-4(2)7-6-3/h1-2H3,(H,5,6,7)', 'InChI=1S/C9H9N3/c1-8-11-10-7-12(8)9-5-3-2-4-6-9/h2-7H,1H3', 'InChI=1S/C18H17N5O/c1-23-15(24)18(22-16(23)19)14-5-3-2-4-11(14)6-17(18)7-12-9-20-21-10-13(12)8-17/h2-5,9-10H,6-8H2,1H3,(H2,19,22)', 'InChI=1S/C24H20N4O/c1-28-15-25-24(22(28)29)21-9-17(16-5-3-2-4-6-16)7-8-18(21)10-23(24)11-19-13-26-27-14-20(19)12-23/h2-9,13-15H,10-12H2,1H3', 'InChI=1S/C5H8N2/c1-5-3-4-7(2)6-5/h3-4H,1-2H3', 'InChI=1S/C5H8N2/c1-4-3-5(2)7-6-4/h3H,1-2H3,(H,6,7)', 'InChI=1S/C8H13N3/c1-6-7(5-10-11-6)8-3-2-4-9-8/h5,8-9H,2-4H2,1H3,(H,10,11)', 'InChI=1S/C5H8N2/c1-5-3-4-7(2)6-5/h3-4H,1-2H3', 'InChI=1S/C5H8N2/c1-5-3-4-6-7(5)2/h3-4H,1-2H3', 'InChI=1S/C8H13N3/c1-11-6-7(5-10-11)8-3-2-4-9-8/h5-6,8-9H,2-4H2,1H3', 'InChI=1S/C5H8N2/c1-4-3-5(2)7-6-4/h3H,1-2H3,(H,6,7)', 'InChI=1S/C5H8N2/c1-5-3-4-6-7(5)2/h3-4H,1-2H3', 'InChI=1S/C8H13N3/c1-6-7(5-10-11-6)8-3-2-4-9-8/h5,8-9H,2-4H2,1H3,(H,10,11)', 'InChI=1S/C3H7NO/c1-4(2)3-5/h3H,1-2H3', 'InChI=1S/C5H6N2OS/c1-7(4-8)5-6-2-3-9-5/h2-4H,1H3', 'InChI=1S/C3H7NO/c1-3(5)4-2/h1-2H3,(H,4,5)', 'InChI=1S/C15H15NO2/c1-3-6-11-12-9-14(16(2)15(11)17)18-13-8-5-4-7-10(12)13/h3-8,12,14H,1,9H2,2H3', 'InChI=1S/C13H13NO3/c1-14-12-6-9(10(7-15)13(14)16)8-4-2-3-5-11(8)17-12/h2-5,7,9,12,15H,6H2,1H3', 'InChI=1S/C15H15NO2/c1-3-6-11-12-9-14(16(2)15(11)17)18-13-8-5-4-7-10(12)13/h3-8,12,14H,1,9H2,2H3', 'InChI=1S/C13H13NO3/c1-14-12-6-9(10(7-15)13(14)16)8-4-2-3-5-11(8)17-12/h2-5,7,9,12,15H,6H2,1H3', 'InChI=1S/C14H15NO2/c1-9-11-8-14(2,15(3)13(9)16)17-12-7-5-4-6-10(11)12/h4-7,11H,1,8H2,2-3H3', 'InChI=1S/C15H15NO2/c1-3-6-11-12-9-15(2,16-14(11)17)18-13-8-5-4-7-10(12)13/h3-8,12H,1,9H2,2H3,(H,16,17)', 'InChI=1S/C13H13NO3/c1-13-6-9(10(7-15)12(16)14-13)8-4-2-3-5-11(8)17-13/h2-5,7,9,15H,6H2,1H3,(H,14,16)', 'InChI=1S/C15H15NO2/c1-3-6-11-12-9-15(2,16-14(11)17)18-13-8-5-4-7-10(12)13/h3-8,12H,1,9H2,2H3,(H,16,17)', 'InChI=1S/C13H13NO3/c1-13-6-9(10(7-15)12(16)14-13)8-4-2-3-5-11(8)17-13/h2-5,7,9,15H,6H2,1H3,(H,14,16)', 'InChI=1S/C14H15NO2/c1-9-11-8-14(2,15(3)13(9)16)17-12-7-5-4-6-10(11)12/h4-7,11H,1,8H2,2-3H3', 'InChI=1S/C6H8N2/c1-5-2-3-8-6(7)4-5/h2-4H,1H3,(H2,7,8)', 'InChI=1S/C8H9NO/c1-6-4-2-3-5-7(6)8(9)10/h2-5H,1H3,(H2,9,10)', 'InChI=1S/C8H10/c1-7-4-3-5-8(2)6-7/h3-6H,1-2H3', 'InChI=1S/C8H10/c1-7-4-3-5-8(2)6-7/h3-6H,1-2H3', 'InChI=1S/C8H9NO/c1-6-2-4-7(5-3-6)8(9)10/h2-5H,1H3,(H2,9,10)', 'InChI=1S/C8H10/c1-7-4-3-5-8(2)6-7/h3-6H,1-2H3', 'InChI=1S/C8H10/c1-7-4-3-5-8(2)6-7/h3-6H,1-2H3', 'InChI=1S/C5H7NS/c1-4-3-7-5(2)6-4/h3H,1-2H3']

    print(len(output_mol_list))

    for idx, x in enumerate(output_mol_list):

        mol = convert_fragment_molecule_to_mol(x, fragment_database)
        inchi = molecule_to_inchi(mol)

        assert inchi == inchi_list[idx]

    count = extend_molecule_list_depth_count([parent], bond_frequencies, fragment_database_graph, 2)

    assert count == 46


def test_extend_molecule_list_database11_20_threshold():

    bond_frequencies = get_bond_frequencies('../datasets/database1000/frequencies_11-20.txt')
    bond_frequencies = bond_frequencies_to_np(bond_frequencies)

    fragment_database = get_fragment_database('../datasets/database1000/fragments_11-20.sdf')

    fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    bond_frequencies = convert_bond_freq_np_to_dict(fragment_database_graph, bond_frequencies, sort_dict=False)

    parent = FragmentMolecule()

    parent.add_fragment(0, [0])

    output_mol_list = extend_molecule_list_depth([parent], bond_frequencies, fragment_database_graph, depth=2, unique=False, threshold=0.1)

    inchi_list = ['InChI=1S/C2H6O/c1-3-2/h1-2H3']

    print(len(output_mol_list))

    for idx, x in enumerate(output_mol_list):

        mol = convert_fragment_molecule_to_mol(x, fragment_database)
        inchi = molecule_to_inchi(mol)

        print(inchi, x.get_build_probability())

        assert inchi == inchi_list[idx]

    output_mol_list = extend_molecule_list_depth([parent], bond_frequencies, fragment_database_graph, depth=2, unique=False, threshold=0.025)

    inchi_list = ['InChI=1S/C2H6O/c1-3-2/h1-2H3', 'InChI=1S/C10H9NO/c1-12-9-4-5-10-8(7-9)3-2-6-11-10/h2-7H,1H3', 'InChI=1S/C14H11NO/c1-16-12-6-7-14-11(9-12)8-10-4-2-3-5-13(10)15-14/h2-9H,1H3', 'InChI=1S/C7H8O/c1-8-7-5-3-2-4-6-7/h2-6H,1H3', 'InChI=1S/C2H6O/c1-3-2/h1-2H3', 'InChI=1S/C2H7N/c1-3-2/h3H,1-2H3', 'InChI=1S/C2H7N/c1-3-2/h3H,1-2H3', 'InChI=1S/C9H9N3/c1-8-11-10-7-12(8)9-5-3-2-4-6-9/h2-7H,1H3', 'InChI=1S/C18H17N5O/c1-23-15(24)18(22-16(23)19)14-5-3-2-4-11(14)6-17(18)7-12-9-20-21-10-13(12)8-17/h2-5,9-10H,6-8H2,1H3,(H2,19,22)', 'InChI=1S/C24H20N4O/c1-28-15-25-24(22(28)29)21-9-17(16-5-3-2-4-6-16)7-8-18(21)10-23(24)11-19-13-26-27-14-20(19)12-23/h2-9,13-15H,10-12H2,1H3', 'InChI=1S/C3H7NO/c1-3(5)4-2/h1-2H3,(H,4,5)', 'InChI=1S/C6H8N2/c1-5-2-3-8-6(7)4-5/h2-4H,1H3,(H2,7,8)', 'InChI=1S/C8H9NO/c1-6-4-2-3-5-7(6)8(9)10/h2-5H,1H3,(H2,9,10)', 'InChI=1S/C8H10/c1-7-4-3-5-8(2)6-7/h3-6H,1-2H3', 'InChI=1S/C8H10/c1-7-4-3-5-8(2)6-7/h3-6H,1-2H3', 'InChI=1S/C5H7NS/c1-4-3-7-5(2)6-4/h3H,1-2H3']

    print(len(output_mol_list))

    all_inchi = []

    for idx, x in enumerate(output_mol_list):

        mol = convert_fragment_molecule_to_mol(x, fragment_database)
        inchi = molecule_to_inchi(mol)
        all_inchi.append(inchi)
        print(inchi, x.get_build_probability())

        assert inchi == inchi_list[idx]


def test_unique():

    bond_frequencies = get_bond_frequencies('../datasets/database1000/frequencies_11-20.txt')
    bond_frequencies = bond_frequencies_to_np(bond_frequencies)

    fragment_database = get_fragment_database('../datasets/database1000/fragments_11-20.sdf')

    fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    bond_frequencies = convert_bond_freq_np_to_dict(fragment_database_graph, bond_frequencies)

    for i in range(len(fragment_database_graph)):
        print(fragment_database_graph.fragments[i].get_canonical_mapping())

    parent = FragmentMolecule()

    parent.add_fragment(0, [0])

    output_mol_list = extend_molecule_list_depth([parent], bond_frequencies, fragment_database_graph, depth=2, threshold=0.1)

    inchi_list = ['InChI=1S/C2H6O/c1-3-2/h1-2H3']

    print(len(output_mol_list))

    for idx, x in enumerate(output_mol_list):

        mol = convert_fragment_molecule_to_mol(x, fragment_database)
        inchi = molecule_to_inchi(mol)

        print(inchi, x.get_build_probability())

        assert inchi == inchi_list[idx]

    output_mol_list = extend_molecule_list_depth([parent], bond_frequencies, fragment_database_graph, depth=2, threshold=0.025, unique=True)

    inchi_list = ['InChI=1S/C2H6O/c1-3-2/h1-2H3', 'InChI=1S/C2H7N/c1-3-2/h3H,1-2H3', 'InChI=1S/C6H8N2/c1-5-2-3-8-6(7)4-5/h2-4H,1H3,(H2,7,8)', 'InChI=1S/C5H7NS/c1-4-3-7-5(2)6-4/h3H,1-2H3', 'InChI=1S/C8H9NO/c1-6-4-2-3-5-7(6)8(9)10/h2-5H,1H3,(H2,9,10)', 'InChI=1S/C8H10/c1-7-4-3-5-8(2)6-7/h3-6H,1-2H3', 'InChI=1S/C8H10/c1-7-4-3-5-8(2)6-7/h3-6H,1-2H3', 'InChI=1S/C10H9NO/c1-12-9-4-5-10-8(7-9)3-2-6-11-10/h2-7H,1H3', 'InChI=1S/C14H11NO/c1-16-12-6-7-14-11(9-12)8-10-4-2-3-5-13(10)15-14/h2-9H,1H3', 'InChI=1S/C7H8O/c1-8-7-5-3-2-4-6-7/h2-6H,1H3', 'InChI=1S/C2H6O/c1-3-2/h1-2H3', 'InChI=1S/C9H9N3/c1-8-11-10-7-12(8)9-5-3-2-4-6-9/h2-7H,1H3', 'InChI=1S/C18H17N5O/c1-23-15(24)18(22-16(23)19)14-5-3-2-4-11(14)6-17(18)7-12-9-20-21-10-13(12)8-17/h2-5,9-10H,6-8H2,1H3,(H2,19,22)', 'InChI=1S/C24H20N4O/c1-28-15-25-24(22(28)29)21-9-17(16-5-3-2-4-6-16)7-8-18(21)10-23(24)11-19-13-26-27-14-20(19)12-23/h2-9,13-15H,10-12H2,1H3', 'InChI=1S/C3H7NO/c1-3(5)4-2/h1-2H3,(H,4,5)']

    hash_list = [239998280776683951529592014534165366710, 91943601169760982729833042586913080476, 121312088074307383177989702408239765377, 176585861146790402643445344288443295727, 61148750547223184390181491815826269078, 23483152691085918865390409116152194287, 98704469994647752752207753422223087024, 56983286536328984378243478896262177430, 229238699357066649893448545566232091737, 202501438282453391564172265968299788733, 290703596703230795744410943336175365962, 207099755104167868250596398622667517365, 265109281295404771413057048832864289724, 311771311390935554596904190555042004370, 144495658023143108760484860213117052764]

    for idx, x in enumerate(output_mol_list):

        mol = convert_fragment_molecule_to_mol(x, fragment_database)
        inchi = molecule_to_inchi(mol)
        assert hash_list[idx] == x.__hash__()
        assert inchi == inchi_list[idx]
        lines = molecule_to_sdf(mol)
        inchi_list.append(inchi)
        hash_list.append(x.__hash__())

def test_unique2():

    bond_frequencies = get_bond_frequencies('../datasets/database1000/frequencies_11-20.txt')
    bond_frequencies = bond_frequencies_to_np(bond_frequencies)

    fragment_database = get_fragment_database('../datasets/database1000/fragments_11-20.sdf')

    fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    bond_frequencies = convert_bond_freq_np_to_dict(fragment_database_graph, bond_frequencies)

    parent = FragmentMolecule()

    parent.add_fragment(0, [0], {0:0})

    output_mol_list = extend_molecule_list_depth([parent], bond_frequencies, fragment_database_graph, depth=3, fragment_database=fragment_database, threshold=0.0025)

    all_inchi = []

    outfile = open('test.sdf', 'w')

    hash_list = []

    for idx, x in enumerate(output_mol_list):

        mol = convert_fragment_molecule_to_mol(x, fragment_database)
        inchi = molecule_to_inchi(mol)
        all_inchi.append(inchi)
        #print(inchi, x.get_build_probability())
        #assert hash_list[idx] == x.get_hash()
        hash_list.append(x.__hash__())
        #assert inchi == inchi_list[idx]
        save_mol_to_sdf(mol, 'test.sdf')

    hash_count = {}

    for i in hash_list:
        if i in hash_count:
            hash_count[i] += 1
        else:
            hash_count[i] = 1

    for i in hash_list:
        print(i, hash_count[i])

    inchi_count = {}

    for i in all_inchi:
        if i in inchi_count:
            inchi_count[i] += 1
        else:
            inchi_count[i] = 1

    for i in all_inchi:
        print(i, inchi_count[i])

    print(len(hash_list), len(all_inchi))

    for i in range(len(hash_list)):

        hash_i = hash_list[i]
        inchi_i = all_inchi[i]

        if hash_count[hash_i] != inchi_count[inchi_i]:
            print(i+1, hash_list[i], hash_count[hash_i], all_inchi[i], inchi_count[inchi_i])

    print(hash_list)
    print(all_inchi)

    hash_list_check = ['bac56cc4e599146549e99566852af1f5', 'b323c7f5d2b6446c1d86217f00b01c51', 'a6b075be9ea1a50893630b55642ff14d', 'ff58fbc35063c4be1e233496f1237ef7', '6a938454aa92bee71a64d0df2bd845d1', 'ff074a9b5a837c8261fd48abd0d33263', '79203c4e37beda86430c6136a7c2c6e2', 'f5ae7816b28c35e2a8a2256391c69c52', '5cf30d27a12b082a2442fb5109e43eeb', '88ad3b879ed3103e2e505b3108047d13', '5c4af55b4760df37a29fe15ae1df611c', '75d3ed9edaebc461707235aa8e7be0ad', '5c4af55b4760df37a29fe15ae1df611c', '75d3ed9edaebc461707235aa8e7be0ad', '7ce3102ceb175f0af6b8560ab289a514', '0cf9a99236e56eedaf35b21d8112f45e', '88ad3b879ed3103e2e505b3108047d13', '5c4af55b4760df37a29fe15ae1df611c', '75d3ed9edaebc461707235aa8e7be0ad', '5c4af55b4760df37a29fe15ae1df611c', '75d3ed9edaebc461707235aa8e7be0ad', '7ce3102ceb175f0af6b8560ab289a514', '0cf9a99236e56eedaf35b21d8112f45e', '51d50d1721f9fe22e99a694c609e33d1', '21ffe6731a9b0afd69a753c5ac3b7d8b', '51d50d1721f9fe22e99a694c609e33d1', '21ffe6731a9b0afd69a753c5ac3b7d8b', 'd2c791ad2ec13be3198d2dd298dce26c', '9698c9f64ff08baca2de34caf66d4e05', '3fdd7877d781cb5173a44efb4822295c', '3fdd7877d781cb5173a44efb4822295c', '899a1aa5453d87204c08aa3d8c57460a', '7d4ffa6b3449194229f6fbeb5e8c1772', '76e8bc5220b32f77cc16dcbf427d45a9', '907cb3217c1dfaf8847f060b0c2c3d9d', '53234a85cf852dc4162b4cd354100e0c', '1e6eb8ce351aa84125c4fae8ae5274ba', '53234a85cf852dc4162b4cd354100e0c', '957348cb63a8d38a368eb389040fe49f', '1e6eb8ce351aa84125c4fae8ae5274ba', '957348cb63a8d38a368eb389040fe49f', '0b28304aec29820bd083bdde32dfec60', '53234a85cf852dc4162b4cd354100e0c', '1e6eb8ce351aa84125c4fae8ae5274ba', '53234a85cf852dc4162b4cd354100e0c', '7bdbfcc2ab484abfaf340e9b3437671c', '1e6eb8ce351aa84125c4fae8ae5274ba', '7bdbfcc2ab484abfaf340e9b3437671c', '9991b6ebb127fef3eb80590195321f68', '53234a85cf852dc4162b4cd354100e0c', '957348cb63a8d38a368eb389040fe49f', '53234a85cf852dc4162b4cd354100e0c', '7bdbfcc2ab484abfaf340e9b3437671c', '957348cb63a8d38a368eb389040fe49f', '7bdbfcc2ab484abfaf340e9b3437671c', '997d9422cf2955775a9d37e21745a93f', '3aad144d7df8e1cfe40658d0fda42e37', '50698cd0a1122787f894e4ff10cf4cb4', '3aad144d7df8e1cfe40658d0fda42e37', '50698cd0a1122787f894e4ff10cf4cb4', 'a35494f981d9da0e9ba024a95a0d6344', '5148a12a3cf72899424f3f664a4ca764', 'a35494f981d9da0e9ba024a95a0d6344', '5148a12a3cf72899424f3f664a4ca764', 'a35494f981d9da0e9ba024a95a0d6344', '5148a12a3cf72899424f3f664a4ca764', 'a35494f981d9da0e9ba024a95a0d6344', '5148a12a3cf72899424f3f664a4ca764', 'a35494f981d9da0e9ba024a95a0d6344', '5148a12a3cf72899424f3f664a4ca764', 'a35494f981d9da0e9ba024a95a0d6344', '5148a12a3cf72899424f3f664a4ca764', 'a35494f981d9da0e9ba024a95a0d6344', '5148a12a3cf72899424f3f664a4ca764', 'a35494f981d9da0e9ba024a95a0d6344', '5148a12a3cf72899424f3f664a4ca764', '39194147f588ac5c086ddbebb985cd0a', 'f72da20937eb71daaaf4f0401c58e8be', '9a7d32289f45e5f2e9f2b9fd4c6a3bce', '7491a1a3d4d9af7f14ccf8510af389b7', 'e50ad9b59618f957f04ad0bd057c2e31', 'b03df8ce9e808009a2ef89c02f277a68', 'dd79f1310bce8420b31daa832080a379', 'fc2722480af09d31242747da83475936', 'e50ad9b59618f957f04ad0bd057c2e31', '05c00c9aceb305a3bf3ddca24723edc0', 'b03df8ce9e808009a2ef89c02f277a68', '05c00c9aceb305a3bf3ddca24723edc0', 'e50ad9b59618f957f04ad0bd057c2e31', 'e50ad9b59618f957f04ad0bd057c2e31', 'f82d3e812dbe33b96421810e512d4112', 'e50ad9b59618f957f04ad0bd057c2e31', '05c00c9aceb305a3bf3ddca24723edc0', 'e50ad9b59618f957f04ad0bd057c2e31', '05c00c9aceb305a3bf3ddca24723edc0', '58fe53799f0280fe116ce65fd3a9d8a9', 'e2677b45e450785eb31e56be9e73f56c', 'db5bb7ea47560b41f6351ad441c9149b', '03d576d3a2712aadf723cfa8b3f4387f', '89270169b3a01a3e5a72322be9458724']

    all_inchi_check = ['InChI=1S/C10H10N2O/c1-13-8-5-7-3-2-4-12-10(7)9(11)6-8/h2-6H,11H2,1H3', 'InChI=1S/C14H12N2O/c1-17-9-6-7-13-11(8-9)14(15)10-4-2-3-5-12(10)16-13/h2-8H,1H3,(H2,15,16)', 'InChI=1S/C14H10ClNO/c1-17-12-4-5-13-10(7-12)6-9-2-3-11(15)8-14(9)16-13/h2-8H,1H3', 'InChI=1S/C7H9NO/c1-9-7-4-2-6(8)3-5-7/h2-5H,8H2,1H3', 'InChI=1S/C8H10O/c1-7-3-5-8(9-2)6-4-7/h3-6H,1-2H3', 'InChI=1S/C2H7NO/c1-4-2-3/h2-3H2,1H3', 'InChI=1S/C3H8O/c1-3-4-2/h3H2,1-2H3', 'InChI=1S/C3H8O/c1-3-4-2/h3H2,1-2H3', 'InChI=1S/C8H10O/c1-9-7-8-5-3-2-4-6-8/h2-6H,7H2,1H3', 'InChI=1S/C3H9N/c1-4(2)3/h1-3H3', 'InChI=1S/C3H9N/c1-4(2)3/h1-3H3', 'InChI=1S/C8H11N/c1-9(2)8-6-4-3-5-7-8/h3-7H,1-2H3', 'InChI=1S/C3H9N/c1-4(2)3/h1-3H3', 'InChI=1S/C8H11N/c1-9(2)8-6-4-3-5-7-8/h3-7H,1-2H3', 'InChI=1S/C13H13N/c1-14(12-8-4-2-5-9-12)13-10-6-3-7-11-13/h2-11H,1H3', 'InChI=1S/C8H11N/c1-7-3-5-8(9-2)6-4-7/h3-6,9H,1-2H3', 'InChI=1S/C3H9N/c1-4(2)3/h1-3H3', 'InChI=1S/C3H9N/c1-4(2)3/h1-3H3', 'InChI=1S/C8H11N/c1-9(2)8-6-4-3-5-7-8/h3-7H,1-2H3', 'InChI=1S/C3H9N/c1-4(2)3/h1-3H3', 'InChI=1S/C8H11N/c1-9(2)8-6-4-3-5-7-8/h3-7H,1-2H3', 'InChI=1S/C13H13N/c1-14(12-8-4-2-5-9-12)13-10-6-3-7-11-13/h2-11H,1H3', 'InChI=1S/C8H11N/c1-7-3-5-8(9-2)6-4-7/h3-6,9H,1-2H3', 'InChI=1S/C10H11N3/c1-8-11-12-9(2)13(8)10-6-4-3-5-7-10/h3-7H,1-2H3', 'InChI=1S/C10H11N3/c1-8-11-12-9(2)13(8)10-6-4-3-5-7-10/h3-7H,1-2H3', 'InChI=1S/C10H11N3/c1-8-11-12-9(2)13(8)10-6-4-3-5-7-10/h3-7H,1-2H3', 'InChI=1S/C10H11N3/c1-8-11-12-9(2)13(8)10-6-4-3-5-7-10/h3-7H,1-2H3', 'InChI=1S/C9H8ClN3/c1-7-12-11-6-13(7)9-4-2-8(10)3-5-9/h2-6H,1H3', 'InChI=1S/C10H11N3/c1-8-5-3-4-6-10(8)13-7-11-12-9(13)2/h3-7H,1-2H3', 'InChI=1S/C24H21N5O/c1-29-21(30)24(28-22(29)25)20-9-16(15-5-3-2-4-6-15)7-8-17(20)10-23(24)11-18-13-26-27-14-19(18)12-23/h2-9,13-14H,10-12H2,1H3,(H2,25,28)', 'InChI=1S/C24H21N5O/c1-29-21(30)24(28-22(29)25)20-9-16(15-5-3-2-4-6-15)7-8-17(20)10-23(24)11-18-13-26-27-14-19(18)12-23/h2-9,13-14H,10-12H2,1H3,(H2,25,28)', 'InChI=1S/C41H32N8O2/c1-49-23-44-41(37(49)51)35-11-27(6-8-29(35)13-39(41)16-32-20-47-48-21-33(32)17-39)25-4-2-3-24(9-25)26-5-7-28-12-38(14-30-18-45-46-19-31(30)15-38)40(34(28)10-26)36(50)42-22-43-40/h2-11,18-23H,12-17H2,1H3,(H,42,43,50)', 'InChI=1S/C25H19N5O/c1-30-15-27-25(23(30)31)22-8-18(17-4-2-3-16(7-17)12-26)5-6-19(22)9-24(25)10-20-13-28-29-14-21(20)11-24/h2-8,13-15H,9-11H2,1H3', 'InChI=1S/C34H28N4O4/c1-38-19-35-34(32(38)39)27-10-21(5-6-24(27)13-33(34)14-25-16-36-37-17-26(25)15-33)20-3-2-4-22(9-20)31-18-41-29-11-23-7-8-40-28(23)12-30(29)42-31/h2-6,9-12,16-17,19,31H,7-8,13-15,18H2,1H3', 'InChI=1S/C30H24N4O/c1-34-19-31-30(28(34)35)27-13-23(22-9-5-8-21(12-22)20-6-3-2-4-7-20)10-11-24(27)14-29(30)15-25-17-32-33-18-26(25)16-29/h2-13,17-19H,14-16H2,1H3', 'InChI=1S/C6H10N2/c1-5-4-6(2)8(3)7-5/h4H,1-3H3', 'InChI=1S/C9H15N3/c1-7-8(6-12(2)11-7)9-4-3-5-10-9/h6,9-10H,3-5H2,1-2H3', 'InChI=1S/C6H10N2/c1-5-4-6(2)8(3)7-5/h4H,1-3H3', 'InChI=1S/C9H15N3/c1-6-9(7(2)12-11-6)8-4-3-5-10-8/h8,10H,3-5H2,1-2H3,(H,11,12)', 'InChI=1S/C9H15N3/c1-7-8(6-12(2)11-7)9-4-3-5-10-9/h6,9-10H,3-5H2,1-2H3', 'InChI=1S/C9H15N3/c1-6-9(7(2)12-11-6)8-4-3-5-10-8/h8,10H,3-5H2,1-2H3,(H,11,12)', 'InChI=1S/C9H15N3/c1-7-8(6-10-11-7)9-4-3-5-12(9)2/h6,9H,3-5H2,1-2H3,(H,10,11)', 'InChI=1S/C6H10N2/c1-5-4-6(2)8(3)7-5/h4H,1-3H3', 'InChI=1S/C9H15N3/c1-7-8(6-12(2)11-7)9-4-3-5-10-9/h6,9-10H,3-5H2,1-2H3', 'InChI=1S/C6H10N2/c1-5-4-6(2)8(3)7-5/h4H,1-3H3', 'InChI=1S/C9H15N3/c1-7-8(6-11-12(7)2)9-4-3-5-10-9/h6,9-10H,3-5H2,1-2H3', 'InChI=1S/C9H15N3/c1-7-8(6-12(2)11-7)9-4-3-5-10-9/h6,9-10H,3-5H2,1-2H3', 'InChI=1S/C9H15N3/c1-7-8(6-11-12(7)2)9-4-3-5-10-9/h6,9-10H,3-5H2,1-2H3', 'InChI=1S/C9H15N3/c1-11-5-3-4-9(11)8-6-10-12(2)7-8/h6-7,9H,3-5H2,1-2H3', 'InChI=1S/C6H10N2/c1-5-4-6(2)8(3)7-5/h4H,1-3H3', 'InChI=1S/C9H15N3/c1-6-9(7(2)12-11-6)8-4-3-5-10-8/h8,10H,3-5H2,1-2H3,(H,11,12)', 'InChI=1S/C6H10N2/c1-5-4-6(2)8(3)7-5/h4H,1-3H3', 'InChI=1S/C9H15N3/c1-7-8(6-11-12(7)2)9-4-3-5-10-9/h6,9-10H,3-5H2,1-2H3', 'InChI=1S/C9H15N3/c1-6-9(7(2)12-11-6)8-4-3-5-10-8/h8,10H,3-5H2,1-2H3,(H,11,12)', 'InChI=1S/C9H15N3/c1-7-8(6-11-12(7)2)9-4-3-5-10-9/h6,9-10H,3-5H2,1-2H3', 'InChI=1S/C9H15N3/c1-7-8(6-10-11-7)9-4-3-5-12(9)2/h6,9H,3-5H2,1-2H3,(H,10,11)', 'InChI=1S/C4H9NO/c1-4(6)5(2)3/h1-3H3', 'InChI=1S/C6H8N2OS/c1-5(9)8(2)6-7-3-4-10-6/h3-4H,1-2H3', 'InChI=1S/C4H9NO/c1-4(6)5(2)3/h1-3H3', 'InChI=1S/C6H8N2OS/c1-5(9)8(2)6-7-3-4-10-6/h3-4H,1-2H3', 'InChI=1S/C16H17NO2/c1-4-7-12-13-10-16(2,17(3)15(12)18)19-14-9-6-5-8-11(13)14/h4-9,13H,1,10H2,2-3H3', 'InChI=1S/C14H15NO3/c1-14-7-10(11(8-16)13(17)15(14)2)9-5-3-4-6-12(9)18-14/h3-6,8,10,16H,7H2,1-2H3', 'InChI=1S/C16H17NO2/c1-4-7-12-13-10-16(2,17(3)15(12)18)19-14-9-6-5-8-11(13)14/h4-9,13H,1,10H2,2-3H3', 'InChI=1S/C14H15NO3/c1-14-7-10(11(8-16)13(17)15(14)2)9-5-3-4-6-12(9)18-14/h3-6,8,10,16H,7H2,1-2H3', 'InChI=1S/C16H17NO2/c1-4-7-12-13-10-16(2,17(3)15(12)18)19-14-9-6-5-8-11(13)14/h4-9,13H,1,10H2,2-3H3', 'InChI=1S/C14H15NO3/c1-14-7-10(11(8-16)13(17)15(14)2)9-5-3-4-6-12(9)18-14/h3-6,8,10,16H,7H2,1-2H3', 'InChI=1S/C16H17NO2/c1-4-7-12-13-10-16(2,17(3)15(12)18)19-14-9-6-5-8-11(13)14/h4-9,13H,1,10H2,2-3H3', 'InChI=1S/C14H15NO3/c1-14-7-10(11(8-16)13(17)15(14)2)9-5-3-4-6-12(9)18-14/h3-6,8,10,16H,7H2,1-2H3', 'InChI=1S/C16H17NO2/c1-4-7-12-13-10-16(2,17(3)15(12)18)19-14-9-6-5-8-11(13)14/h4-9,13H,1,10H2,2-3H3', 'InChI=1S/C14H15NO3/c1-14-7-10(11(8-16)13(17)15(14)2)9-5-3-4-6-12(9)18-14/h3-6,8,10,16H,7H2,1-2H3', 'InChI=1S/C16H17NO2/c1-4-7-12-13-10-16(2,17(3)15(12)18)19-14-9-6-5-8-11(13)14/h4-9,13H,1,10H2,2-3H3', 'InChI=1S/C14H15NO3/c1-14-7-10(11(8-16)13(17)15(14)2)9-5-3-4-6-12(9)18-14/h3-6,8,10,16H,7H2,1-2H3', 'InChI=1S/C16H17NO2/c1-4-7-12-13-10-16(2,17(3)15(12)18)19-14-9-6-5-8-11(13)14/h4-9,13H,1,10H2,2-3H3', 'InChI=1S/C14H15NO3/c1-14-7-10(11(8-16)13(17)15(14)2)9-5-3-4-6-12(9)18-14/h3-6,8,10,16H,7H2,1-2H3', 'InChI=1S/C16H17NO2/c1-4-7-12-13-10-16(2,17(3)15(12)18)19-14-9-6-5-8-11(13)14/h4-9,13H,1,10H2,2-3H3', 'InChI=1S/C14H15NO3/c1-14-7-10(11(8-16)13(17)15(14)2)9-5-3-4-6-12(9)18-14/h3-6,8,10,16H,7H2,1-2H3', 'InChI=1S/C15H13N3/c1-11-7-9-16-14(10-11)18-13-6-2-4-12-5-3-8-17-15(12)13/h2-10H,1H3,(H,16,18)', 'InChI=1S/C7H10N2/c1-6-3-4-9-7(5-6)8-2/h3-5H,1-2H3,(H,8,9)', 'InChI=1S/C19H15N3/c1-13-10-11-20-18(12-13)22-19-14-6-2-4-8-16(14)21-17-9-5-3-7-15(17)19/h2-12H,1H3,(H,20,21,22)', 'InChI=1S/C11H11N3/c1-9-5-7-13-11(8-9)14-10-4-2-3-6-12-10/h2-8H,1H3,(H,12,13,14)', 'InChI=1S/C9H11NO/c1-6-3-4-8(9(10)11)7(2)5-6/h3-5H,1-2H3,(H2,10,11)', 'InChI=1S/C9H11NO/c1-6-4-3-5-7(2)8(6)9(10)11/h3-5H,1-2H3,(H2,10,11)', 'InChI=1S/C9H11NO/c1-7-5-3-4-6-8(7)9(11)10-2/h3-6H,1-2H3,(H,10,11)', 'InChI=1S/C12H16N2O/c1-9-4-2-3-5-11(9)12(15)14-10-6-7-13-8-10/h2-5,10,13H,6-8H2,1H3,(H,14,15)', 'InChI=1S/C9H11NO/c1-6-3-4-8(9(10)11)7(2)5-6/h3-5H,1-2H3,(H2,10,11)', 'InChI=1S/C9H12/c1-7-4-8(2)6-9(3)5-7/h4-6H,1-3H3', 'InChI=1S/C9H11NO/c1-6-4-3-5-7(2)8(6)9(10)11/h3-5H,1-2H3,(H2,10,11)', 'InChI=1S/C9H12/c1-7-4-8(2)6-9(3)5-7/h4-6H,1-3H3', 'InChI=1S/C9H11NO/c1-6-3-4-8(9(10)11)7(2)5-6/h3-5H,1-2H3,(H2,10,11)', 'InChI=1S/C9H11NO/c1-6-3-4-8(9(10)11)7(2)5-6/h3-5H,1-2H3,(H2,10,11)', 'InChI=1S/C9H11NO/c1-7-3-5-8(6-4-7)9(11)10-2/h3-6H,1-2H3,(H,10,11)', 'InChI=1S/C9H11NO/c1-6-3-4-8(9(10)11)7(2)5-6/h3-5H,1-2H3,(H2,10,11)', 'InChI=1S/C9H12/c1-7-4-8(2)6-9(3)5-7/h4-6H,1-3H3', 'InChI=1S/C9H11NO/c1-6-3-4-8(9(10)11)7(2)5-6/h3-5H,1-2H3,(H2,10,11)', 'InChI=1S/C9H12/c1-7-4-8(2)6-9(3)5-7/h4-6H,1-3H3', 'InChI=1S/C5H8N2S/c1-4-7-5(2-6)3-8-4/h3H,2,6H2,1H3', 'InChI=1S/C6H9NS/c1-3-6-4-8-5(2)7-6/h4H,3H2,1-2H3', 'InChI=1S/C6H9NS/c1-3-6-4-8-5(2)7-6/h4H,3H2,1-2H3', 'InChI=1S/C6H7NO2S/c1-4-7-5(3-10-4)2-6(8)9/h3H,2H2,1H3,(H,8,9)', 'InChI=1S/C11H11NS/c1-9-12-11(8-13-9)7-10-5-3-2-4-6-10/h2-6,8H,7H2,1H3']

    #assert hash_list == hash_list_check
    #assert all_inchi == all_inchi_check

def test_convert_to_networkx_canonical():

    bond_frequencies = get_bond_frequencies('../datasets/database1000/frequencies_11-20.txt')
    bond_frequencies = bond_frequencies_to_np(bond_frequencies)

    fragment_database = get_fragment_database('../datasets/database1000/fragments_11-20.sdf')

    fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    bond_frequencies = convert_bond_freq_np_to_dict(fragment_database_graph, bond_frequencies)

    f = FragmentMolecule()

    f.add_fragment(0, [0])
    f.add_fragment(31, [0,1,3,5])
    f.add_fragment(0, [0])
    f.add_fragment(11, [0,2])

    f.add_bond(0,1,0,3)
    f.add_bond(1,2,1,0)
    f.add_bond(1,3,0,0)

    print('line722')
    for i in f._graph.fragments:
        print(f._graph.fragments[i].attachment_points)
    print('line725')

    print('line727', f._graph.attachment_point_list)

    mol = convert_fragment_molecule_to_mol(f, fragment_database)

    lines = molecule_to_sdf(mol)
    with open('test2.sdf', 'w') as outfile:
        for line in lines:
            outfile.write(line)
        outfile.write('$$$$\n')

    print(fragment_database_graph.fragments[31].get_canonical_mapping())

    print('line740', f._graph.attachment_point_list)

    for i in f._graph.fragments:
        print(f._graph.fragments[i].attachment_points)


def test_get_unique_molecule_list_sort():

    f = FragmentMolecule()
    f.add_fragment(0, [0])
    f._graph._build_probability = 1

    f2 = FragmentMolecule()
    f2.add_fragment(1, [0])
    f2._graph._build_probability = 2

    f3 = FragmentMolecule()
    f3.add_fragment(2, [0])
    f3._graph._build_probability = 3

    sorted_list = get_unique_molecule_list([f, f2, f3])

    answers = [3,2,1]

    for idx, i in enumerate(sorted_list):
        print(i._graph.build_probability)
        assert i._graph.build_probability == answers[idx]

def test_get_unique_molecule_list_sort2():

    f = FragmentMolecule()
    f.add_fragment(0, [0])
    f._graph._build_probability = 1

    f2 = FragmentMolecule()
    f2.add_fragment(0, [0])
    f2._graph._build_probability = 2

    f3 = FragmentMolecule()
    f3.add_fragment(2, [0])
    f3._graph._build_probability = 3

    sorted_list = get_unique_molecule_list([f, f2, f3])

    answers = [3,3]

    for idx, i in enumerate(sorted_list):
        print(i._graph.build_probability)
        assert i._graph.build_probability == answers[idx]


def test_get_unique_molecule_list_sort3():

    f = FragmentMolecule()
    f.add_fragment(0, [0], {0:0})
    f.add_fragment(0, [0], {0:0})
    f.add_bond(0,1,0,0)
    f._graph._build_probability = 1

    f2 = FragmentMolecule()
    f2.add_fragment(0, [0], {0:0})
    f2.add_fragment(0, [1], {1:0})
    f2.add_bond(0,1,0,1)
    f2._graph._build_probability = 2

    f3 = FragmentMolecule()
    f3.add_fragment(2, [0], {0:0})
    f3._graph._build_probability = 3

    sorted_list = get_unique_molecule_list([f, f2, f3])

    answers = [3,3]

    for idx, i in enumerate(sorted_list):
        print(i._graph.build_probability)
        assert i._graph.build_probability == answers[idx]


def test_convert_bond_freq_np_to_dict():

    bond_frequencies = get_bond_frequencies('../datasets/database1000/frequencies1.txt')
    bond_frequencies_np = bond_frequencies_to_np(bond_frequencies)

    fragment_database = get_fragment_database('../datasets/database1000/fragments1.sdf')
  
    fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    bond_frequencies = convert_bond_freq_np_to_dict(fragment_database_graph, bond_frequencies_np)

    assert bond_frequencies == {(0, 0): {(1, 0): 1}, (1, 0): {(0, 0): 1}, (1, 2): {(2, 1): 1}, (2, 1): {(1, 2): 1}, (2, 2): {(3, 2): 1}, (3, 2): {(2, 2): 1, (4, 0): 1, (5, 3): 1, (6, 1): 1}, (4, 0): {(3, 2): 1, (5, 1): 1}, (5, 1): {(4, 0): 1}, (5, 3): {(3, 2): 1}, (5, 9): {(8, 0): 1}, (5, 11): {(8, 0): 1}, (6, 1): {(3, 2): 1, (7, 3): 1}, (7, 3): {(6, 1): 1}, (8, 0): {(5, 9): 1, (5, 11): 1}}

    bonds = bond_frequencies_np[0]
    freq = bond_frequencies_np[1]

    for key, val in zip(bonds, freq):
        i = key[0]
        j = key[1]
        k = key[2]
        l = key[3]
        assert bond_frequencies[(i,k)][(j,l)] == val

def test_read_write_bond_frequencies_dict():

    bond_frequencies = get_bond_frequencies('../datasets/database1000/frequencies1.txt')
    bond_frequencies_np = bond_frequencies_to_np(bond_frequencies)

    fragment_database = get_fragment_database('../datasets/database1000/fragments1.sdf')
  
    fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    bond_frequencies = convert_bond_freq_np_to_dict(fragment_database_graph, bond_frequencies_np)

    write_bond_frequencies_dict(bond_frequencies, 'bond_frequencies_dict.txt')

    bond_frequencies_read = read_bond_frequencies_dict('bond_frequencies_dict.txt')

    assert bond_frequencies == bond_frequencies_read


def test_prepare_fragment():

    bond_frequencies = get_bond_frequencies('../datasets/database1000/frequencies11-20.txt')
    bond_frequencies_np = bond_frequencies_to_np(bond_frequencies)

    fragment_database = get_fragment_database('../datasets/database1000/fragments11-20.sdf')
  
    fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    bond_frequencies = convert_bond_freq_np_to_dict(fragment_database_graph, bond_frequencies_np)    

    parent_file = '../datasets/database1000/benzene.sdf'
    parent_fragment_file_list = ['../datasets/database1000/ch4.sdf', '../datasets/database1000/nh3.sdf']

    remove_hydrogens = [11, 6]
    remove_hydrogens_parent_fragment = [4, 3]

    parent_mapping_1 = [0, 0, 1, 0]

    parent, bond_frequencies, fragment_database, fragment_database_graph = prepare_parent(bond_frequencies, fragment_database, fragment_database_graph, parent_file, parent_fragment_file_list, parent_mapping_1,remove_hydrogens, remove_hydrogens_parent_fragment)

    print(len(fragment_database))

    assert bond_frequencies == {(0, 0): {(1, 0): 4, (7, 0): 3, (31, 1): 2, (15, 0): 1, (18, 0): 1, (23, 0): 1, (23, 2): 1, (23, 3): 1, (25, 0): 1, (28, 4): 1, (28, 5): 1, (30, 0): 1, (31, 3): 1, (34, 4): 1}, (1, 0): {(0, 0): 4, (2, 2): 1, (5, 11): 1, (10, 4): 1, (4, 1): 1}, (2, 2): {(1, 0): 1}, (2, 4): {(3, 0): 1}, (3, 0): {(4, 1): 3, (2, 4): 1, (5, 0): 1, (30, 4): 1}, (4, 1): {(3, 0): 3, (4, 1): 3, (12, 0): 3, (10, 4): 3, (14, 3): 2, (33, 0): 1, (1, 0): 1, (32, 0): 1, (25, 1): 1, (24, 11): 1, (11, 0): 1, (16, 5): 1, (33, 2): 1, (9, 1): 1, (15, 0): 1, (13, 1): 1, (8, 0): 1, (7, 0): 1, (17, 5): 1, (34, 1): 1}, (5, 0): {(3, 0): 1}, (5, 4): {(6, 0): 1}, (5, 11): {(1, 0): 1}, (6, 0): {(5, 4): 1, (16, 3): 1}, (7, 0): {(0, 0): 3, (10, 4): 2, (4, 1): 1}, (8, 0): {(4, 1): 1}, (8, 5): {(9, 1): 1}, (8, 7): {(9, 1): 1}, (9, 1): {(8, 5): 1, (8, 7): 1, (4, 1): 1, (18, 5): 1}, (10, 4): {(4, 1): 3, (7, 0): 2, (11, 0): 1, (21, 7): 1, (22, 0): 1, (1, 0): 1, (27, 2): 1, (10, 4): 1}, (11, 0): {(10, 4): 1, (4, 1): 1, (31, 0): 1}, (11, 2): {(12, 0): 2, (21, 7): 1}, (12, 0): {(4, 1): 3, (14, 3): 3, (11, 2): 2, (13, 1): 1}, (13, 1): {(4, 1): 1, (12, 0): 1}, (14, 3): {(12, 0): 3, (4, 1): 2, (35, 0): 1}, (15, 0): {(0, 0): 1, (4, 1): 1}, (15, 4): {(16, 0): 1}, (16, 0): {(15, 4): 1}, (16, 3): {(6, 0): 1}, (16, 5): {(4, 1): 1}, (17, 5): {(4, 1): 1, (19, 4): 1}, (18, 0): {(0, 0): 1}, (18, 5): {(9, 1): 1}, (18, 8): {(19, 4): 1}, (19, 4): {(18, 8): 1, (20, 0): 1, (32, 11): 1, (17, 5): 1}, (20, 0): {(19, 4): 1}, (21, 7): {(11, 2): 1, (10, 4): 1}, (22, 0): {(10, 4): 1}, (23, 0): {(0, 0): 1}, (23, 2): {(0, 0): 1}, (23, 3): {(0, 0): 1}, (23, 4): {(24, 7): 1}, (24, 7): {(23, 4): 1}, (24, 11): {(4, 1): 1}, (25, 0): {(0, 0): 1, (26, 2): 1}, (25, 1): {(4, 1): 1}, (26, 2): {(25, 0): 1}, (27, 2): {(10, 4): 1, (28, 0): 1}, (28, 0): {(27, 2): 1, (29, 0): 1}, (28, 4): {(0, 0): 1}, (28, 5): {(0, 0): 1}, (29, 0): {(28, 0): 1}, (30, 0): {(0, 0): 1}, (30, 4): {(3, 0): 1}, (31, 0): {(11, 0): 1}, (31, 1): {(0, 0): 2}, (31, 3): {(0, 0): 1}, (32, 0): {(4, 1): 1}, (32, 11): {(19, 4): 1}, (33, 0): {(4, 1): 1}, (33, 2): {(4, 1): 1}, (34, 1): {(4, 1): 1}, (34, 4): {(0, 0): 1}, (35, 0): {(36, 0): 3, (14, 3): 1}, (36, 0): {(35, 0): 3}, (37, 0): {(1, 0): 4, (7, 0): 3, (31, 1): 2, (15, 0): 1, (18, 0): 1, (23, 0): 1, (23, 2): 1, (23, 3): 1, (25, 0): 1, (28, 4): 1, (28, 5): 1, (30, 0): 1, (31, 3): 1, (34, 4): 1}, (37, 1): {(8, 5): 1, (8, 7): 1, (4, 1): 1, (18, 5): 1}}

    assert parent.get_frag_id(0) == 37
    assert parent.get_canonical_mapping(0) == {0: 0, 1: 0, 2: 2, 3: 3, 4: 3, 5: 2, 7: 7, 8: 8, 9: 8, 10: 7}
    assert parent.list_free_valence_points() == [[0,1]]
    print(len(fragment_database))
    print(len(fragment_database_graph))

    assert fragment_database_graph.fragments[37].attachment_points == [0,1]
    assert fragment_database_graph.fragments[37].get_attribute('frag_id') == 37 
    assert fragment_database_graph.fragments[37].get_canonical_mapping() == {0: 0, 1: 0, 2: 2, 3: 3, 4: 3, 5: 2, 7: 7, 8: 8, 9: 8, 10: 7}

    # add mol not working in fragment database

test_prepare_fragment()