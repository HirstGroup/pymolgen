import filecmp
import os
import sys

def test():

	os.system('python ../inchi2svg.py -i test_inchi2svg_input.inchi')

	assert filecmp.cmp('test_inchi2svg_input.can', 'test_inchi2svg_model.txt')

test()

