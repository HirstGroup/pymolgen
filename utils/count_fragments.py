import argparse

from rdkit import Chem

from pymolgen.canonicalise_tautomer import canonicalise_tautomer
from pymolgen.molecule_formats import molecule_from_smiles
from pymolgen.fragment_mol import get_fragments_dataset

parser = argparse.ArgumentParser(description='Count number of fragments in generated molecules')

parser.add_argument('-i','--input', help='Input file, inchi format', required=True)
parser.add_argument('-o','--output', help='Output file', required=True)

args = parser.parse_args()

with open(args.input) as infile, open(args.output, 'w') as outfile:

	for line in infile:

		inchi = line.split()[0]

		rdmol = Chem.MolFromInchi(inchi)

		smi = Chem.MolToSmiles(rdmol)

		smi = canonicalise_tautomer(smi)

		mol = molecule_from_smiles(smi)

		fragments, pairs, bonds = get_fragments_dataset(mol, carbonyl=True, fluorine=True)

		print(len(fragments))



