import filecmp
import os
import sys

sys.path.append('../../')
from utils.classify_meta_para import *

def test1():

	os.system('python ../classify_meta_para.py -i input/inchi10.inchi -fi inchi -o output/inchi10_classify.inchi')

	assert filecmp.cmp('output/inchi10_classify.inchi', 'input/inchi10_classify.inchi') is True


def test2():

	inchi = 'InChI=1S/C11H11NO/c1-8-11(9(2)13-12-8)10-6-4-3-5-7-10/h3-7H,1-2H3'

	mol_type = classify(inchi)

	assert mol_type == 'NONE'

	inchi = 'InChI=1S/C24H22N2O/c1-16-11-12-25-23(13-16)22-10-9-20(24-17(2)26-27-18(24)3)15-21(22)14-19-7-5-4-6-8-19/h4-13,15H,14H2,1-3H3'

	mol_type = classify(inchi)

	assert mol_type == 'parameta'

	inchi = 'InChI=1S/C12H13NO/c1-8-5-4-6-11(7-8)12-9(2)13-14-10(12)3/h4-7H,1-3H3'

	mol_type = classify(inchi)

	assert mol_type == 'meta'

	inchi = 'InChI=1S/C12H13NO/c1-8-4-6-11(7-5-8)12-9(2)13-14-10(12)3/h4-7H,1-3H3'

	mol_type = classify(inchi)

	assert mol_type == 'para'


def test3():
	# test csv input

	df = pd.read_csv('input/calculator10.csv', sep=';')

	df['mol_type'] = df.apply(classify_row, axis=1)

	print(df)


def test4():
	# test csv input command line

	os.system('python ../classify_meta_para.py -i input/calculator10.csv -fi csv -o output/calculator10_classify.csv')

	assert filecmp.cmp('input/calculator10_classify.csv', 'output/calculator10_classify.csv') is True



