import filecmp
import os
import sys

def test_probabilities():

	os.system('python ../remove_duplicates_list.py -i test_remove_duplicates_list_input.txt -o test_remove_duplicates_list_output.txt')

	assert filecmp.cmp('test_remove_duplicates_list_model.txt', 'test_remove_duplicates_list_output.txt') is True

