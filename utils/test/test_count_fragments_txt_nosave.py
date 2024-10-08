import subprocess
import sys

sys.path.append('../')
from count_fragments_txt_nosave import *


def test1():

	depth_count = count_fragments_txt('input/phenylisoxazole_random_0.txt')

	print(depth_count)

	assert depth_count == {3: 3, 4: 3, 5: 2, 6: 1, 7: 1, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0}

test1()
