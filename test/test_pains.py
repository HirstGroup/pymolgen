from pymolgen.properties_pymolgen import *

def test_pains():

	pains_database = gen_pains_database()

	pains_smi = []

	with open('test_pains.csv') as infile:
		next(infile)
		for line in infile:
			smi = line.split(';')[1].strip('\n')
			pains_smi.append(smi)

	for smi in pains_smi:
		print(pains_filter_rdkit(smi, pains_database))

test_pains()