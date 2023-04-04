#!/usr/bin/python

import sys,os
import argparse
import numpy as np
import pandas as pd

#   CHIRAL_THRESHOLD:   Maximum number of Chiral Centers (<= 2)
CHIRAL_THRESHOLD = 2
#   PSA_THRESHOLD:       Polar surface area (<= 140)
PSA_THRESHOLD = 140
#   PFI_THRESHOLD:       Property Forecast Index (LOGP + # of aromatic rings < 8)
PFI_THRESHOLD = 8
#   ROTBOND_THRESHOLD:  Number of rotatable bonds (<= 7)
ROTBOND_THRESHOLD = 7
#   WEIGHT_THRESHOLD:    Maximum MW in Daltons (<= 500)
WEIGHT_THRESHOLD = 500
#   H_DON_THRESHOLD:     Maximum number of hydrogen donors (<= 5)
H_DON_THRESHOLD = 5
#   H_ACC_THRESHOLD:     Maximum number of hydrogen acceptors (<= 10)
H_ACC_THRESHOLD = 10
#   LOGP_THRESHOLD:      Water/Octanol Partition Coefficient (0.5-5.0)
LOGP_THRESHOLD_UP = 5
LOGP_THRESHOLD_LOW = 0.5

parser = argparse.ArgumentParser(description='Recalculate filter_pass')
parser.add_argument('-i','--input', help='Input file name', required=True)
parser.add_argument('-o','--output', help='Output file name', required=True)

args = parser.parse_args()

input = args.input
output = args.output

print(input, output)

df = pd.read_csv(input, sep=';') #, dtype={'n_chiral': int})

def filter(row):

	if np.isnan(row['n_chiral']) or row['n_chiral'] > CHIRAL_THRESHOLD:
		return False
	if np.isnan(row['psa']) or row['psa'] > PSA_THRESHOLD:
		return False
	if np.isnan(row['pfi']) or row['pfi'] > PFI_THRESHOLD:
		return False
	if np.isnan(row['n_rot_bonds']) or row['n_rot_bonds'] > ROTBOND_THRESHOLD:
		return False
	if np.isnan(row['mw']) or row['mw'] > WEIGHT_THRESHOLD:
		return False
	if np.isnan(row['h_don']) or row['h_don'] > H_DON_THRESHOLD:
		return False
	if np.isnan(row['h_acc']) or row['h_acc'] > H_ACC_THRESHOLD:
		return False
	if np.isnan(row['logp']) or row['logp'] > LOGP_THRESHOLD_UP:
		return False
	if np.isnan(row['logp']) or row['logp'] < LOGP_THRESHOLD_LOW:
		return False

	return True

if 'filter_pass' in df.columns:
	df['filter_pass_original'] = df['filter_pass']	

df['filter_pass'] = df.apply (lambda row: filter(row), axis=1)

df.sort_values(['filter_pass','mpo'], ascending=[False, True], inplace=True)

df.to_csv(output, sep=';', index=False)