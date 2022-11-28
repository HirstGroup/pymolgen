import sys,os

import pandas as pd

input = sys.argv[1]
output = sys.argv[2]

print(input, output)

df = pd.read_csv(input, sep=';')

outfile = open(output, 'w')

inchi_list = df['inchi'].tolist()

for i in inchi_list:
	outfile.write('%s\n' %i)

