import subprocess

from pymolgen.analysis.analysis_fragment_molecule import *


def test1():

	subprocess.run('python ../analysis_fragment_molecule.py -i input/inchi10.inchi -o output/inchi10_analysis.txt -a ../../datasets/fragments/fragments_30_50k_co_10_l5_5_sorted_filter_copy.sdf -p input/phenylisoxazole.sdf -r 20 21 -rf ../../datasets/fragments/fragments_30_50k_co_10_l5_5_sorted_filter_copy.txt', check=True, shell=True)


test1()