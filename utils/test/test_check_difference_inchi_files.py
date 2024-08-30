import filecmp
import os
import sys


def test1():

	os.system('python ../check_difference_inchi_files.py -i1 input/test_check_difference_inchi_files1.inchi -i2 input/test_check_difference_inchi_files2.inchi -o1 output/test_check_difference_inchi_files1.inchi -o2 output/test_check_difference_inchi_files2.inchi')

	with open('output/test_check_difference_inchi_files1.inchi') as f1, open('output/test_check_difference_inchi_files1.inchi') as f2:
		lines1 = f1.readlines()
		lines2 = f2.readlines()

		assert sorted(lines1) == sorted(lines2)

	with open('output/test_check_difference_inchi_files2.inchi') as f1, open('output/test_check_difference_inchi_files2.inchi') as f2:
		lines1 = f1.readlines()
		lines2 = f2.readlines()

		assert sorted(lines1) == sorted(lines2)


def test2():

	os.system('python ../check_difference_inchi_files.py -i1 input/test_check_difference_inchi_files1.inchi -i2 input/test_check_difference_inchi_files1.inchi -o1 output/test_check_difference_inchi_files1.inchi -o2 output/test_check_difference_inchi_files2.inchi')

	with open('output/test_check_difference_inchi_files1.inchi') as f:

		assert len(f.readlines()) == 0

	with open('output/test_check_difference_inchi_files2.inchi') as f:

		assert len(f.readlines()) == 0


def test3():

	os.system('python ../check_difference_inchi_files.py -i1 input/test_check_difference_inchi_files1.inchi -i2 input/test_check_difference_inchi_files3.inchi -o1 output/test_check_difference_inchi_files1.inchi -o2 output/test_check_difference_inchi_files2.inchi')

	with open('output/test_check_difference_inchi_files1.inchi') as f:

		assert len(f.readlines()) == 0

	with open('output/test_check_difference_inchi_files2.inchi') as f:

		assert len(f.readlines()) == 10

