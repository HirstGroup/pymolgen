import os
import sys

sys.path.append('../../')
from utils.split_molecule_list import *


def test1():

    fragment_molecule_list = []

    with open('input/test_split_molecule_list.txt') as infile:

        for line in infile:

            fragment_molecule_list.append(line.strip())

    fragment_molecule_list_list = split_molecule_list(fragment_molecule_list, 5)    

    print(fragment_molecule_list_list)

    for i in fragment_molecule_list_list:

        total = 0

        for j in i:

            total += int(j.split(':')[2])

        print(total)


def test2():

    fragment_molecule_list = []

    with open('input/test_split_molecule_list.txt') as infile:

        for line in infile:

            fragment_molecule_list.append(line.strip())

    fragment_molecule_list_list = divide_into_n_lists(fragment_molecule_list, n=5)

    print(fragment_molecule_list_list)

    for i in fragment_molecule_list_list:

        total = 0

        for j in i:

            total += int(j.split(':')[2])

        print(total)

    assert fragment_molecule_list_list == [['10:10:10', '2:2:2', '1:1:1', '1:1:1', '1:1:1', '1:1:1', '1:1:1'], ['9:9:9', '2:2:2', '1:1:1', '1:1:1', '1:1:1', '1:1:1', '1:1:1', '1:1:1'], ['8:8:8', '2:2:2', '1:1:1', '1:1:1', '1:1:1', '1:1:1', '1:1:1', '1:1:1', '1:1:1'], ['7:7:7', '2:2:2', '2:2:2', '1:1:1', '1:1:1', '1:1:1', '1:1:1', '1:1:1', '1:1:1'], ['6:6:6', '5:5:5', '1:1:1', '1:1:1', '1:1:1', '1:1:1', '1:1:1', '1:1:1']]