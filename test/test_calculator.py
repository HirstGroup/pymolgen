import filecmp
import os
import sys

def test():

	os.system('python ../calculator.py -i inputs/test_calculator.inchi -o outputs/test_calculator.txt')

	assert filecmp.cmp('outputs/test_calculator.txt', 'models/test_calculator.txt')

