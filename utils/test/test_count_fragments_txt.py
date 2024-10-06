import filecmp
import os
import shutil
import subprocess
import sys

def test1():

	subprocess.run('python ../count_fragments_txt.py -a input/phenylisoxazole_random_0.txt -i input/phenylisoxazole_random_0.inchi -o output/phenylisoxazole_random_0_count.inchi', check=True, shell=True)


def test_delete():

	shutil.copy('input/phenylisoxazole_random_0.inchi', 'output/')

	subprocess.run('python ../count_fragments_txt.py -a input/phenylisoxazole_random_0.txt -i output/phenylisoxazole_random_0.inchi -o output/phenylisoxazole_random_0_count.inchi --delete', check=True, shell=True)

	assert os.path.isfile('output/phenylisoxazole_random_0.inchi') is False