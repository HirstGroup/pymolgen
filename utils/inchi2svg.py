#!/usr/bin/env python

import argparse
import os
import sys

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

from functools import partial

print = partial(print, flush=True)

# Canonicalise tautomer:
def canonicalise_tautomer(smi, method='RDKit'):
    """
    Convert a SMILES to a canonical tautomer, 
    using RDKit, OpenBabel or via InChIs.
    """

    if method == 'RDKit':
        enumerator = rdMolStandardize.TautomerEnumerator()
        # enumerator.SetReassignStereo = True
        # enumerator.SetRemoveBondStereo = False
        # enumerator.SetRemoveSp3Stereo = False
        mol = enumerator.Canonicalize(Chem.MolFromSmiles(smi))
        canon_smi = Chem.MolToSmiles(mol, canonical=True, 
                                     isomericSmiles=True)

        # If no change in SMILES, revert to original SMILES to retain
        # stereochemistry which may have been lost in 
        # TautomerEnumerator():
        canon_smi_2D = Chem.MolToSmiles(mol, canonical=True, 
                                        isomericSmiles=False)
        smi_2D = Chem.MolToSmiles(Chem.MolFromSmiles(smi),
                                  canonical=True, isomericSmiles=False)
        if canon_smi != smi and canon_smi_2D == smi_2D:
            canon_smi = smi

    elif method == 'InChI':
        # Convert to and from InChI to standarise tautomer:
        mol_smi = Chem.MolFromSmiles(smi)
        inchi = Chem.MolToInchi(mol_smi)
        mol_inchi = Chem.MolFromInchi(inchi)
        canon_smi = Chem.MolToSmiles(mol_inchi, canonical=True, 
                                     isomericSmiles=True)

    elif method == 'OpenBabel':
        # Could use otautomer to get obabel canonical tautomers
        print('Warning: Method not yet implemented', file=sys.stderr)
        canon_smi = smi

    elif method == 'OpenEye':
        # Could use openeye to get a canonical tautomer
        print('Warning: Method not yet implemented', file=sys.stderr)
        canon_smi = smi

    return canon_smi


parser = argparse.ArgumentParser(description='Convert inchi to svg with canonical tautomer')
parser.add_argument('-i', '--input', help='Inchi file',required=True)

args = parser.parse_args()

name = os.path.splitext(args.input)[0]

os.system(f'obabel {args.input} -O {name}.smi')

with open(name + '.smi') as infile, open(name + '.can', 'w') as smifile:

	for line in infile:
		smi = line.strip().split()[0]

		try:
			smi_can = canonicalise_tautomer(smi)
		except:
			smi_can = smi

		smifile.write(f'{smi_can}\n')

os.system(f'obabel {name}.can -O {name}.svg')

