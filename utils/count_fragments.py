import argparse

import pandas as pd

from rdkit import Chem

from pymolgen.canonicalise_tautomer import canonicalise_tautomer
from pymolgen.molecule_formats import molecule_from_smiles
from pymolgen.fragment_mol import get_fragments_dataset

parser = argparse.ArgumentParser(description='Count number of fragments in generated molecules')

parser.add_argument('-i','--input', help='Input csv File Name, separated by semicolon', required=True)
parser.add_argument('-o','--output', help='Output file name, csv format', required=True)

args = parser.parse_args()


def count_fragments(row):

		inchi = row['inchi']

		rdmol = Chem.MolFromInchi(inchi)

		smi = Chem.MolToSmiles(rdmol)

		smi = canonicalise_tautomer(smi)

		mol = molecule_from_smiles(smi)

		fragments, pairs, bonds = get_fragments_dataset(mol, carbonyl=True, fluorine=False)
		fragments_fluorine, pairs, bonds = get_fragments_dataset(mol, carbonyl=True, fluorine=True)

		patt = Chem.MolFromSmarts("[#9]")

		n_fluorine = len(rdmol.GetSubstructMatches(patt))

		print(n_fluorine)

		n_fluorine_notcf3 = n_fluorine - len(fragments) + len(fragments_fluorine)

		return pd.Series([len(fragments_fluorine), n_fluorine_notcf3])


df = pd.read_csv(args.input, sep=';')

df[['n_fragments', 'n_fluorine_notcf3']] = df.apply(count_fragments, axis=1)

print(df)



