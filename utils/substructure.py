import os
import sys

from rdkit import Chem

input = sys.argv[1]
output = sys.argv[2]
pattern = sys.argv[3]

smi_list = []

with open(input) as f:
	for line in f:
		smi_list.append(line.strip('\n'))

with open(output, 'w') as f:
	for idx, x in enumerate(smi_list):
		m = Chem.MolFromSmiles(x)
		patt = Chem.MolFromSmarts(pattern)
		if m.HasSubstructMatch(patt):
			print(idx, x)
			f.write(f'{x} {idx}\n')


"""
patt = Chem.MolFromSmarts(smi)
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
"""