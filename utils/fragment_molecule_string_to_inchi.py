import argparse

from pymolgen.fragment_molecule import *
from pymolgen.fragment_molecule_builder import *

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Build Molecules using the FragmentMolecule class')

    # required arguments
    parser.add_argument('-a','--fragments_sdf', help='SDF file of fragments',required=True)
    parser.add_argument('-i','--input', help='Input file name with FragmentMolecule strings',required=True)
    parser.add_argument('-o','--output', help='Output file name in inchi format',required=True)
    parser.add_argument('-p','--parent_file', help='Parent Structure File in SDF format',required=True)
    parser.add_argument('-r','--remove_hydrogens', type=int, nargs='+', help='Space-separated hydrogen atoms that will be created as attachment points, numbered from 0',required=True)

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

            string_representation = line.strip()

            build_probability = string_representation.split(';')[2]

            fragment_molecule = generate_fragment_molecule_from_string(string_representation, fragment_database_graph)

            mol = convert_fragment_molecule_to_mol(fragment_molecule, fragment_database)

            inchi = molecule_to_inchi(mol)

            outfile.write(f'{inchi} {build_probability}\n')


