import os
import subprocess
import sys

def test1():

	subprocess.run('python ../find_missing_inchis.py -i input/inchi5.inchi -a input/inchi10.inchi -o output/test1.inchi', check=True, shell=True)


def test2():

	subprocess.run('python ../find_missing_inchis.py -i input/inchi10.inchi -a input/inchi5.inchi -o output/test2.inchi', check=True, shell=True)


def test3():

	subprocess.run('python ../find_missing_inchis.py -i input/inchi10.inchi -a input/inchi10_1.inchi  input/inchi10_2.inchi -o output/test2.inchi', check=True, shell=True)