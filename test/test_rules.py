import filecmp
import os
import sys

dir_path = os.path.dirname(os.path.realpath(__file__))

def test1():

	os.system('rm outputs/*')

	os.system('cp inputs/all_chembl_30_10k_calculator_ok_100.csv outputs/')

	os.chdir('outputs')

	os.system(f'python {dir_path}/../rules.py -i all_chembl_30_10k_calculator_ok_100.csv -o all_chembl_30_10k_calculator_ok_100_rules.csv -r all_chembl_30_10k_calculator_ok_100_rules.smi --all')

	os.chdir('../')

	assert filecmp.cmp('inputs/all_chembl_30_10k_calculator_ok_100_rules.csv', 'outputs/all_chembl_30_10k_calculator_ok_100_rules.csv') is True
	assert filecmp.cmp('inputs/all_chembl_30_10k_calculator_ok_100_rules.smi', 'outputs/all_chembl_30_10k_calculator_ok_100_rules.smi') is True

