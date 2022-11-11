import sys,os

from pymolgen.utils.smi_funcs import canonicalise_tautomer

infile = open(sys.argv[1])
outfile = open(sys.argv[2], 'w')

for line in infile:
	smi = line.strip('\n')

	smi, warnings = canonicalise_tautomer(smi)

	outfile.write('%s\n' %smi)