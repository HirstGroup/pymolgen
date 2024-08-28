import argparse

from pymolgen.fragment_molecule import *

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Build Molecules using the FragmentMolecule class')

    # required arguments
    parser.add_argument('-a','--fragments_sdf', help='SDF file of fragments',required=True)
    parser.add_argument('-i','--input', help='Input file name with FragmentMolecule strings',required=True)
    parser.add_argument('-o','--output', help='Output file name in inchi format',required=True)

    args = parser.parse_args()

    if args.input == args.output:
    	raise Exception('Same input and output')

    fragment_database = get_fragment_database(args.fragments_sdf)

    fragment_database_graph = convert_fragment_database_to_graph(fragment_database)

    with open(args.input) as infile, open(args.output, 'w') as outfile:

    	for line in infile:

    		string_representation = line.strip()

    		fragment_molecule = generate_fragment_molecule_from_string(string_representation, fragment_database_graph)

    		mol = convert_fragment_molecule_to_mol(fragment_molecule, fragment_database)

    		inchi = molecule_to_inchi(mol)

    		outfile.write(f'{inchi}\n')


