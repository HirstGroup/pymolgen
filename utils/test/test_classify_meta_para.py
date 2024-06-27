import filecmp
import os
import sys

def test1():

	os.system(f'python ../classify_meta_para.py -i input/inchi10.inchi -o output/inchi10_classify.inchi')

	assert filecmp.cmp('output/inchi10_classify.inchi', 'input/inchi10_classify.inchi') is True

test1()

