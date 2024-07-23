import os
import sys

sys.path.append('../../')
from utils.classify_meta_para import *
from utils.substructure import *

def test1():

	os.system("python ../substructure.py -i input/inchi10.inchi -fi inchi -o1 output/inchi10_substructure_true.inchi -o2 output/inchi10_substructure_false.inchi -p '[#6]-[#6]-1-[#8]-[#7]=[#6](-[#6])-[#6]-1-[#6]-1=[#6]-[#6]=[#6]-[#6]=[#6]-1'")


def test2():

	pattern = '[cH]1[cH]c(-*)[cH][cH]c1-c2c(-[CH3])onc2-[CH3]'

	os.system(f"python ../substructure.py -i input/inchi10.inchi -fi inchi -o1 output/inchi10_substructure_true.inchi -o2 output/inchi10_substructure_false.inchi -p '{pattern}'")


def test3():

	inchi = 'InChI=1S/C11H11NO/c1-8-11(9(2)13-12-8)10-6-4-3-5-7-10/h3-7H,1-2H3'

	mol_type = classify(inchi)

	assert mol_type == 'NONE'

	inchi = 'InChI=1S/C24H22N2O/c1-16-11-12-25-23(13-16)22-10-9-20(24-17(2)26-27-18(24)3)15-21(22)14-19-7-5-4-6-8-19/h4-13,15H,14H2,1-3H3'

	mol_type = classify(inchi)

	assert mol_type == 'parameta'


def test_load_patterns():

	pattern_list = ['smarts', '[cH]1[cH]c(-*)[cH][cH]c1-c2c(-[CH3])onc2-[CH3]']

	pattern_list_objects = load_patterns(pattern_list)


def test_read_input_inchi():

	mol_list = read_input('input/inchi10.inchi', 'inchi')

	inchi_list = []

	with open('input/inchi10.inchi') as f:

		for line in f:

			inchi_list.append(line.strip().split()[0])

	assert len(mol_list) == len(inchi_list)

	for mol, inchi in zip(mol_list, inchi_list):

		assert Chem.MolToInchi(mol) == inchi


def test_read_input_smi():

	mol_list = read_input('input/inchi10.smi', 'smi')

	smi_list = []

	with open('input/inchi10.smi') as f:

		for line in f:

			smi_list.append(line.strip().split()[0])

	assert len(mol_list) == len(smi_list)

	for mol, smi in zip(mol_list, smi_list):

		assert Chem.MolToInchi(mol) == Chem.MolToInchi(Chem.MolFromSmiles(smi))


def test_find_substructure():

	mol_list = read_input('input/inchi10.inchi', 'inchi')
	pattern_list_objects = load_patterns(['smarts', '[cH]1c(-*)[cH][cH][cH]c1-c2c(-[CH3])onc2-[CH3]'])

	output_false, output_true = find_substructure(mol_list, pattern_list_objects)

	print(output_false)
	print(output_true)

	assert output_false == ['InChI=1S/C11H11NO/c1-8-11(9(2)13-12-8)10-6-4-3-5-7-10/h3-7H,1-2H3', 'InChI=1S/C24H22N2O/c1-16-11-12-25-23(13-16)22-10-9-20(24-17(2)26-27-18(24)3)15-21(22)14-19-7-5-4-6-8-19/h4-13,15H,14H2,1-3H3', 'InChI=1S/C18H19N3O2/c1-11-18(12(2)23-21-11)13-4-9-16(19-3)17(10-13)20-14-5-7-15(22)8-6-14/h4-10,19-20,22H,1-3H3', 'InChI=1S/C24H24N4O/c1-15-24(16(2)29-28-15)19-10-9-18(11-17-7-8-17)20(12-19)13-26-23-14-25-21-5-3-4-6-22(21)27-23/h3-6,9-10,12,14,17H,7-8,11,13H2,1-2H3,(H,26,27)', 'InChI=1S/C12H13NO/c1-8-4-6-11(7-5-8)12-9(2)13-14-10(12)3/h4-7H,1-3H3', 'InChI=1S/C13H15NO/c1-8-5-6-12(7-9(8)2)13-10(3)14-15-11(13)4/h5-7H,1-4H3', 'InChI=1S/C15H19NO/c1-5-12-7-8-14(9-13(12)6-2)15-10(3)16-17-11(15)4/h7-9H,5-6H2,1-4H3', 'InChI=1S/C13H15NO/c1-4-11-5-7-12(8-6-11)13-9(2)14-15-10(13)3/h5-8H,4H2,1-3H3']

	assert output_true == ['InChI=1S/C12H13NO/c1-8-5-4-6-11(7-8)12-9(2)13-14-10(12)3/h4-7H,1-3H3', 'InChI=1S/C13H15NO/c1-4-11-6-5-7-12(8-11)13-9(2)14-15-10(13)3/h5-8H,4H2,1-3H3']


def test_find_substructure2():

	#pattern = '[#6]-[#6]-1-[#8]-[#7]=[#6](-[#6])-[#6]-1-[#6]-1=[#6]-[#6]=[#6]-[#6]=[#6]-1'
	pattern = '[#6]-[#6]-[#8]-[#7]=[#6](-[#6])-[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]'
	#pattern = 'CC1ON=C(C)C1C1=CC=CC=C1'
	#pattern = 'CC1ON=C(C)C1C1~CC~CC~C1'

	mol_list = read_input('input/inchi10.inchi', 'inchi')
	pattern_list_objects = load_patterns([pattern], 'smarts')

	output_false, output_true = find_substructure(mol_list, pattern_list_objects)

	print(output_false)
	print(output_true)

test_find_substructure2()
