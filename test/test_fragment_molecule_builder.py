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

    ch3 = FragmentMolecule()

    ch3.add_fragment(0, [0])

    output_mol_list = extend_molecule_list_depth([ch3], bond_frequencies, fragment_database_graph, depth=1)

    for x in output_mol_list:
        mol = convert_fragment_molecule_to_mol(x, fragment_database)
        inchi = molecule_to_inchi(mol)
        assert str(x) == '0-1'
        assert inchi == 'InChI=1S/C4H5NO/c1-4-2-3-5-6-4/h2-3H,1H3'

    output_mol_list = extend_molecule_list_depth([ch3], bond_frequencies, fragment_database_graph, depth=2)

    for x in output_mol_list:
        mol = convert_fragment_molecule_to_mol(x, fragment_database)
        inchi = molecule_to_inchi(mol)
        assert str(x) == '0-1-2'
        assert inchi == 'InChI=1S/C5H6N2O2/c1-4-2-5(6-3-8)7-9-4/h2-3H,1H3,(H,6,7,8)'

    output_mol_list = extend_molecule_list_depth([ch3], bond_frequencies, fragment_database_graph, depth=3)

    for x in output_mol_list:
        mol = convert_fragment_molecule_to_mol(x, fragment_database)
        inchi = molecule_to_inchi(mol)
        assert str(x) == '0-1-2-3'
        assert inchi == 'InChI=1S/C6H8N2O2/c1-4-3-6(8-10-4)7-5(2)9/h3H,1-2H3,(H,7,8,9)'

    output_mol_list = extend_molecule_list_depth([ch3], bond_frequencies, fragment_database_graph, depth=4)

    answers = ['InChI=1S/C7H9N3O3/c1-4-2-6(10-13-4)9-7(12)3-5(8)11/h2H,3H2,1H3,(H2,8,11)(H,9,10,12)', 'InChI=1S/C6H8N2O4S/c1-4-2-5(8-12-4)7-6(9)3-13(10)11/h2,13H,3H2,1H3,(H,7,8,9)', 'InChI=1S/C14H13N3O2/c1-10-8-13(16-19-10)15-14(18)9-17-7-6-11-4-2-3-5-12(11)17/h2-8H,9H2,1H3,(H,15,16,18)', 'InChI=1S/C12H12N2O2/c1-9-7-11(14-16-9)13-12(15)8-10-5-3-2-4-6-10/h2-7H,8H2,1H3,(H,13,14,15)']

    for idx, x in enumerate(output_mol_list):
        mol = convert_fragment_molecule_to_mol(x, fragment_database)
        inchi = molecule_to_inchi(mol)
        print(inchi)
        assert inchi == answers[idx]

    output_mol_list = extend_molecule_list_depth([ch3], bond_frequencies, fragment_database_graph, depth=5)

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

    ch3 = FragmentMolecule()

    ch3.add_fragment(0, [0])

    output_mol_list = extend_molecule_list_depth([ch3], bond_frequencies, fragment_database_graph, depth=1)

    answers = ['0-0', '0-1']

    for j in range(len(output_mol_list)):
        assert str(output_mol_list[j]) == answers[j]

def test_extend_molecule_list_depth_rzt():

    bond_frequencies = get_bond_frequencies('../datasets/simple/frequencies_rzt.txt')
    print(bond_frequencies)
    bond_frequencies = bond_frequencies_to_np(bond_frequencies)

    fragment_database = get_fragment_database('../datasets/simple/fragments_rzt.sdf')
    fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    ch3 = FragmentMolecule()

    ch3.add_fragment(0, [0])

    output_mol_list = extend_molecule_list_depth([ch3], bond_frequencies, fragment_database_graph, depth=1)

    for j in output_mol_list:
        mol = convert_fragment_molecule_to_mol(j, fragment_database)
        inchi = molecule_to_inchi(mol)
        assert inchi == 'InChI=1S/C7H8/c1-7-5-3-2-4-6-7/h2-6H,1H3'

    output_mol_list = extend_molecule_list_depth([ch3], bond_frequencies, fragment_database_graph, depth=2)

    answers = ['InChI=1S/C8H10/c1-7-3-5-8(2)6-4-7/h3-6H,1-2H3', 'InChI=1S/C13H12/c1-11-7-9-13(10-8-11)12-5-3-2-4-6-12/h2-10H,1H3']

    for idx, x in enumerate(output_mol_list):
        mol = convert_fragment_molecule_to_mol(x, fragment_database)
        inchi = molecule_to_inchi(mol)
        assert inchi == answers[idx]

    output_mol_list = extend_molecule_list_depth([ch3], bond_frequencies, fragment_database_graph, depth=3)

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

    output_mol_list = extend_molecule_list_depth([f], bond_frequencies, fragment_database_graph, depth=1)

    assert len(output_mol_list) == 1
    assert str(output_mol_list[0]) == '0-0'
    assert output_mol_list[0].list_bonds() == [(0, 1, 1, 1)]

    output_mol_list = extend_molecule_list_depth([f], bond_frequencies, fragment_database_graph, depth=2)

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

    output_mol_list = extend_molecule_list_depth([f], bond_frequencies, fragment_database_graph, depth=1)

    for i in output_mol_list:
        print(i, i.list_bonds())    

    assert len(output_mol_list) == 1
    assert str(output_mol_list[0]) == '0-0'
    assert output_mol_list[0].list_bonds() == [(0, 1, 2, 2)]

    output_mol_list = extend_molecule_list_depth([f], bond_frequencies, fragment_database_graph, depth=2)

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

    output_mol_list = extend_molecule_list_depth([f], bond_frequencies, fragment_database_graph, depth=1)

    assert len(output_mol_list) == 2

    mols = ['0-0', '0-0']
    bonds = [[(0, 1, 1, 2)], [(0, 1, 2, 1)]]

    for idx, x in enumerate(output_mol_list):
        print(x, x.list_bonds())
        assert str(x) == mols[idx]
        assert x.list_bonds() == bonds[idx]

    output_mol_list = extend_molecule_list_depth([f], bond_frequencies, fragment_database_graph, depth=2)

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

    output_mol_list = extend_molecule_list_depth([f], bond_frequencies, fragment_database_graph, depth=1)

    #assert len(output_mol_list) == 2

    mols = ['0-0', '0-0','0-0']
    bonds = [[(0, 1, 1, 1)], [(0, 1, 1, 2)], [(0, 1, 2, 1)]]

    for idx, x in enumerate(output_mol_list):
        print(x, x.list_bonds())
        assert str(x) == mols[idx]
        assert x.list_bonds() == bonds[idx]

    output_mol_list = extend_molecule_list_depth([f], bond_frequencies, fragment_database_graph, depth=2)

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

    output_mol_list = extend_molecule_list_depth([f], bond_frequencies, fragment_database_graph, depth=1)

    mols = ['0-0', '0-0']
    bonds = [[(0, 1, 1, 2)], [(0, 1, 2, 2)]]

    for idx, x in enumerate(output_mol_list):
        print(x, x.list_bonds())
        assert str(x) == mols[idx]
        assert x.list_bonds() == bonds[idx]

    output_mol_list = extend_molecule_list_depth([f], bond_frequencies, fragment_database_graph, depth=2)

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

    output_mol_list = extend_molecule_list_depth([f], bond_frequencies, fragment_database_graph, depth=1)

    mols = ['0-0', '0-0', '0-0']
    bonds = [[(0, 1, 1, 2)], [(0, 1, 2, 2)], [(0, 1, 3, 2)]]

    for idx, x in enumerate(output_mol_list):
        print(x, x.list_bonds())
        assert str(x) == mols[idx]
        assert x.list_bonds() == bonds[idx]

    output_mol_list = extend_molecule_list_depth([f], bond_frequencies, fragment_database_graph, depth=2)

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

    parent = FragmentMolecule()

    parent.add_fragment(0, [0])

    output_mol_list = extend_molecule_list_depth([parent], bond_frequencies, fragment_database_graph, 2)

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

    output_mol_list = extend_molecule_list_depth([parent], bond_frequencies, fragment_database_graph, depth=2, threshold=0.025)

    inchi_list = ['InChI=1S/C2H6O/c1-3-2/h1-2H3', 'InChI=1S/C10H9NO/c1-12-9-4-5-10-8(7-9)3-2-6-11-10/h2-7H,1H3', 'InChI=1S/C14H11NO/c1-16-12-6-7-14-11(9-12)8-10-4-2-3-5-13(10)15-14/h2-9H,1H3', 'InChI=1S/C7H8O/c1-8-7-5-3-2-4-6-7/h2-6H,1H3', 'InChI=1S/C2H6O/c1-3-2/h1-2H3', 'InChI=1S/C2H7N/c1-3-2/h3H,1-2H3', 'InChI=1S/C2H7N/c1-3-2/h3H,1-2H3', 'InChI=1S/C9H9N3/c1-8-11-10-7-12(8)9-5-3-2-4-6-9/h2-7H,1H3', 'InChI=1S/C18H17N5O/c1-23-15(24)18(22-16(23)19)14-5-3-2-4-11(14)6-17(18)7-12-9-20-21-10-13(12)8-17/h2-5,9-10H,6-8H2,1H3,(H2,19,22)', 'InChI=1S/C24H20N4O/c1-28-15-25-24(22(28)29)21-9-17(16-5-3-2-4-6-16)7-8-18(21)10-23(24)11-19-13-26-27-14-20(19)12-23/h2-9,13-15H,10-12H2,1H3', 'InChI=1S/C3H7NO/c1-3(5)4-2/h1-2H3,(H,4,5)', 'InChI=1S/C6H8N2/c1-5-2-3-8-6(7)4-5/h2-4H,1H3,(H2,7,8)', 'InChI=1S/C8H9NO/c1-6-4-2-3-5-7(6)8(9)10/h2-5H,1H3,(H2,9,10)', 'InChI=1S/C8H10/c1-7-4-3-5-8(2)6-7/h3-6H,1-2H3', 'InChI=1S/C8H10/c1-7-4-3-5-8(2)6-7/h3-6H,1-2H3', 'InChI=1S/C5H7NS/c1-4-3-7-5(2)6-4/h3H,1-2H3']

    print(len(output_mol_list))

    all_inchi = []

    for idx, x in enumerate(output_mol_list):

        mol = convert_fragment_molecule_to_mol(x, fragment_database)
        inchi = molecule_to_inchi(mol)
        all_inchi.append(inchi)
        print(inchi, x.get_build_probability())

        assert inchi == inchi_list[idx]

test_extend_molecule_list_database11_20_threshold()