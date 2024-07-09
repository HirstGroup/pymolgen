import filecmp
import os
import sys


def test():
	# test doesn't work for build_probabilities but works for everything else

	os.system('python ../calculator.py -i inputs/test_calculator.inchi -o outputs/test_calculator.txt -l outputs/test_calculator.log -bp')

	assert filecmp.cmp('outputs/test_calculator.txt', 'models/test_calculator.txt')


def test2():
	
	os.system('python ../calculator.py -i inputs/test_calculator_classify.inchi -o outputs/test_calculator_classify.txt --classify -l outputs/test_calculator_classify.log')

	assert filecmp.cmp('outputs/test_calculator_classify.txt', 'models/test_calculator_classify.txt')
	assert filecmp.cmp('outputs/test_calculator_classify.log', 'models/test_calculator_classify.log')

test2()