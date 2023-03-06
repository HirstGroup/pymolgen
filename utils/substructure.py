import os
import sys

from rdkit import Chem

input = sys.argv[1]
output = sys.argv[2]

patt = Chem.MolFromSmarts('O=C1Nc2c(C1)cccc2')
matches = []

suppl = Chem.SDMolSupplier(input)

for mol in suppl:
	if mol is not None:
		if mol.HasSubstructMatch(patt):
			matches.append(mol)

print(len(matches))

w = Chem.SDWriter(output)

for m in matches:
	w.write(m)
