import filecmp
import os
import shutil
import subprocess

def test1():

	subprocess.run('python ../split_depth.py -i input/test_split_depth_input.txt -o output/test_split_depth1', check=True, shell=True)

	for i in range(3,13):

		assert filecmp.cmp(f'input/test_split_depth1_d{i}.inchi', f'output/test_split_depth1_d{i}.inchi') is True

def test_delete():

	shutil.copy('input/test_split_depth_input.txt', 'output/')

	subprocess.run('python ../split_depth.py -i output/test_split_depth_input.txt -o output/test_split_depth1 --delete', check=True, shell=True)

	for i in range(3,13):

		assert filecmp.cmp(f'input/test_split_depth1_d{i}.inchi', f'output/test_split_depth1_d{i}.inchi') is True

	assert os.path.isfile('output/test_split_depth_input.txt') is False	