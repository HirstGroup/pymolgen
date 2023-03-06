import argparse

from pymolgen.fragment_mol import *

parser = argparse.ArgumentParser(description='Fragment molecules in database')
parser.add_argument('-i','--input', help='Input file',required=True)
parser.add_argument('-o','--output', help='Output file',required=True)
parser.add_argument('-f', '--fragment', help='Fragment file')

args = parser.parse_args()

dataset = SDFDatasetLarge(args.input)

frag = molecule_from_sdf(args.fragment)
frag_len = len(frag)

outfile = open(args.output, 'w')

for i in range(len(dataset)):
    mol = dataset.load_molecule(i)

    fragments, pairs, bonds = get_fragments_dataset(mol)

    for j in fragments:

        if len(j) != frag_len:
            continue

        gm = isomorphism.GraphMatcher(j, frag, node_match=node_compare_element)

        if gm.is_isomorphic():
            lines = molecule_to_sdf(j)
            for line in lines:
                outfile.write(line)
            outfile.write('$$$$\n')