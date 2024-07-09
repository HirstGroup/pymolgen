import argparse

from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')

parser = argparse.ArgumentParser(description='Remove salts from molecules')
parser.add_argument('-i','--input', help='Input file in inchi format',required=True)
parser.add_argument('-o','--output', help='Output file, inchi format',required=True)
parser.add_argument('-a1','--aux1', help='Auxiliary file 1 for changed inchis, inchi format',required=True)
parser.add_argument('-a2','--aux2', help='Auxiliary file 2 for failed inchis, inchi format',required=True)

args = parser.parse_args()

if args.input == args.output:
	sys.exit('Same input as output')

infile = open(args.input)
outfile = open(args.output, 'w')

# write changed inchis in auxfile1
auxfile1 = open(args.aux1, 'w')

# write failed inchis in auxfile2
auxfile2 = open(args.aux2, 'w')

remover = SaltRemover()

for line in infile:

	inchi = line.strip().split()[0]

	try:

		mol = Chem.MolFromInchi(inchi)

		mol_stripped = remover.StripMol( mol )

		inchi_stripped = Chem.MolToInchi(mol_stripped)

		if inchi != inchi_stripped:
			auxfile1.write(f'{inchi}\n{inchi_stripped}\n')

		outfile.write(f'{inchi_stripped}\n')

	except:

		auxfile2.write(f'{inchi}\n')