import os
import sys

from rdkit import Chem

import argparse

 
parser = argparse.ArgumentParser(description='This script splits SDF file')
parser.add_argument('-i','--input', help='Input file name', required=True)
parser.add_argument('-p','--pattern_list', nargs='+', help='Smarts pattern list space separated, use individual quotation marks for each pattern', required=True)
parser.add_argument('-o1','--output_true', help='Output file name of molecules that contain pattern', required=True)
parser.add_argument('-o2','--output_false', help='Output file name of molecules that do not contain pattern', required=True)
args = parser.parse_args()

pattern_list = args.pattern_list


pattern_list_objects = [Chem.MolFromSmarts(x) for x in pattern_list]

smi_list = []

with open(args.input) as f:
	for line in f:
		smi_list.append(line.strip('\n').split()[0])

with open(args.output_true, 'w') as f, open(args.output_false, 'w') as f2:
	for idx, x in enumerate(smi_list):
		m = Chem.MolFromSmiles(x)
		match = False
		for pattern_i in pattern_list_objects:
			if m.HasSubstructMatch(pattern_i):
				print(idx, x)
				f.write(f'{x} {idx}\n')
				match = True
		if match is False:
			f2.write(f'{x} {idx}\n')