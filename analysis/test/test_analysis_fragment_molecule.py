import filecmp
import subprocess

from pymolgen.analysis.analysis_fragment_molecule import *


def test1():

    subprocess.run('python ../analysis_fragment_molecule.py -i input/inchi1.inchi -o output/inchi1_analysis.txt -a ../../datasets/fragments/fragments_30_50k_co_10_l5_5_sorted_filter_copy.sdf -p input/phenylisoxazole.sdf -r 20 21 -rf ../../datasets/fragments/fragment_database_30_50k_co_10_l5_5_sorted_filter_copy.txt', check=True, shell=True)


def test2():

    subprocess.run('python ../analysis_fragment_molecule.py -i input/inchi10.inchi -o output/inchi10_analysis.txt -a ../../datasets/fragments/fragments_30_50k_co_10_l5_5_sorted_filter_copy.sdf -p input/phenylisoxazole.sdf -r 20 21 -rf ../../datasets/fragments/fragment_database_30_50k_co_10_l5_5_sorted_filter_copy.txt', check=True, shell=True)

    assert filecmp.cmp('input/inchi10_analysis.txt', 'output/inchi10_analysis.txt') is True


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


def test_calculate_build_probability():

    fragment_database = get_fragment_database('../../datasets/fragments/fragments_30_50k_co_10_l5_5_sorted_filter_copy.sdf')

    fragment_database_graph = read_fragment_database_graph('../../datasets/fragments/fragment_database_30_50k_co_10_l5_5_sorted_filter_copy.txt')

    bond_frequencies = read_bond_frequencies_dict('../../datasets/fragments/bond_frequencies_30_50k_co_10_l5_5_sorted_filter_copy.txt')

    parent, bond_frequencies, fragment_database, fragment_database_graph = prepare_parent(bond_frequencies=bond_frequencies, fragment_database=fragment_database, fragment_database_graph=fragment_database_graph, parent_file='input/phenylisoxazole.sdf', parent_fragment_file_list=['input/benzene.sdf', 'input/benzene.sdf'], parent_mapping_1=[15,0,16,0], remove_hydrogens=[20,21], remove_hydrogens_parent_fragment=[6,6])

    with open('input/calculate10.txt') as f:

        for line in f:

            string_representation = line.strip()

            print(string_representation)

            fragment_molecule = generate_fragment_molecule_from_string(string_representation, fragment_database_graph)

            build_probability = calculate_build_probability(bond_frequencies, fragment_database_graph, fragment_molecule, root=16310)

            print(fragment_molecule._graph.build_probability, build_probability)


def test_calculate_build_probability2():

    fragment_database = get_fragment_database('../../datasets/database1000/fragments10.sdf')

    fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    bond_frequencies = get_bond_frequencies('../../datasets/database1000/frequencies10.txt')
    bond_frequencies = bond_frequencies_to_np(bond_frequencies)
    bond_frequencies = convert_bond_freq_np_to_dict(fragment_database_graph, bond_frequencies)

    string_representation = '2-3:(0,1,1,2):1.0'

    print(string_representation)

    fragment_molecule = generate_fragment_molecule_from_string(string_representation, fragment_database_graph)

    build_probability = calculate_build_probability(bond_frequencies, fragment_database_graph, fragment_molecule, root=2)

    print(fragment_molecule._graph.build_probability, build_probability)


def test_calculate_build_probability3():

    fragment_database = get_fragment_database('../../datasets/database1000/fragments10.sdf')

    fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    bond_frequencies = get_bond_frequencies('../../datasets/database1000/frequencies10.txt')
    bond_frequencies = bond_frequencies_to_np(bond_frequencies)
    bond_frequencies = convert_bond_freq_np_to_dict(fragment_database_graph, bond_frequencies)

    with open('input/builder-depth3.txt') as infile:

        for n, line in enumerate(infile):

            string_representation = line.strip()

            fragment_molecule = generate_fragment_molecule_from_string(string_representation, fragment_database_graph)

            build_probability = calculate_build_probability(bond_frequencies, fragment_database_graph, fragment_molecule, root=0)

            print('RESULT', n+1, string_representation, fragment_molecule._graph.build_probability, build_probability)

            #assert (fragment_molecule._graph.build_probability - build_probability) ** 2 < 0.00001

test_calculate_build_probability3()

def test_bfs():

    import networkx as nx

    # Create a directed graph (can also be undirected)
    G = nx.Graph()

    # Adding edges to create a graph structure
    G.add_edges_from([(0, 1), (0, 2), (1, 3), (1, 4), (2, 5), (2, 6), (5, 7)])

    # Root node for BFS
    root = 0

    # Perform BFS traversal from root
    def bfs_with_edges(graph, root):
        # Perform BFS using nx.bfs_tree, which returns nodes in BFS order
        bfs_nodes = list(nx.bfs_tree(graph, root))
        
        visited = set()

        for node in bfs_nodes:

            visited.add(node)

            print(f"Node {node}:")
            # List all edges originating from the current node
            for neighbor in graph.neighbors(node):
                if neighbor not in visited:
                    print(f"  Edge: ({node}, {neighbor})")

    # Run the BFS with edge listing
    bfs_with_edges(G, root)


def test_get_ordered_bonds():

    f = FragmentMolecule()

    f.add_fragment(0, [0])
    f.add_fragment(10, [1,2,3])
    f.add_fragment(20, [5, 6])
    f.add_fragment(40, [7, 8])
    f.add_fragment(30, [9])
    f.add_fragment(50, [4])

    f.add_bond(0, 1, 0, 1)
    f.add_bond(1, 2, 2, 5)
    f.add_bond(2, 3, 6, 7)
    f.add_bond(3, 4, 8, 9)
    f.add_bond(1, 5, 3, 4)

    bonds = get_ordered_bonds(f, root=0)

    print(bonds)

    assert bonds == [(0, 1), (1, 2), (1, 5), (2, 3), (3, 4)]

