import os
import sys

import pandas as pd


def test1():

	os.system('python ../count_fragments.py -i input/all_chembl_30_10k_calculator_ok_100_rules.csv -o output/all_chembl_30_10k_calculator_ok_100_rules_count.csv')


def test2():

	def function_row(row):

		return pd.Series([0, 1])

	df = pd.read_csv('input/all_chembl_30_10k_calculator_ok_100_rules.csv', sep=';')

	df[['row1', 'row2']] = df.apply(function_row, axis=1)


test1()




