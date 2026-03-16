#!/usr/bin/env python3
import sys,os

# Import RDKit Tools
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def sub_search(smi, sub_smi):

    m = Chem.MolFromSmiles(smi)
    m = Chem.AddHs(m)

    sub = Chem.MolFromSmiles(sub_smi)

    patt_match = m.HasSubstructMatch(sub)

    if patt_match is True:
        return True

    return False

input = sys.argv[1]
sub_smi = sys.argv[2]

with open(input) as infile:
    for line in infile:
        smi = line.strip()

        if sub_search(smi, sub_smi):
            print(smi)

