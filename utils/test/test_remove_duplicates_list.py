import filecmp
import os
import shutil
import subprocess
import sys

def test1():

	subprocess.run('python ../remove_duplicates_list.py -i input/test_remove_duplicates_list_input2.txt input/test_remove_duplicates_list_input3.txt -o output/test_remove_duplicates_list_output1.txt --read_count', check=True, shell=True)

	assert filecmp.cmp('input/test_remove_duplicates_list_output1.txt', 'output/test_remove_duplicates_list_output1.txt') is True


def test_delete():

	shutil.copy('input/test_remove_duplicates_list_input2.txt', 'output')
	shutil.copy('input/test_remove_duplicates_list_input3.txt', 'output')

	subprocess.run('python ../remove_duplicates_list.py -i output/test_remove_duplicates_list_input2.txt output/test_remove_duplicates_list_input3.txt -o output/test_remove_duplicates_list_output1.txt --delete --read_count', check=True, shell=True)

	assert filecmp.cmp('input/test_remove_duplicates_list_output1.txt', 'output/test_remove_duplicates_list_output1.txt') is True

	assert os.path.isfile('output/test_remove_duplicates_list_input1.txt') is False
	assert os.path.isfile('output/test_remove_duplicates_list_input2.txt') is False


def test_delete_nocount():

	shutil.copy('input/test_remove_duplicates_list_input2.txt', 'output')
	shutil.copy('input/test_remove_duplicates_list_input3.txt', 'output')

	subprocess.run('python ../remove_duplicates_list.py -i output/test_remove_duplicates_list_input2.txt output/test_remove_duplicates_list_input3.txt -o output/test_remove_duplicates_list_output2.txt --delete', check=True, shell=True)

	assert filecmp.cmp('input/test_remove_duplicates_list_output2.txt', 'output/test_remove_duplicates_list_output2.txt') is True

	assert os.path.isfile('output/test_remove_duplicates_list_input1.txt') is False
	assert os.path.isfile('output/test_remove_duplicates_list_input2.txt') is False

