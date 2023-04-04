#!/usr/bin/python

import sys,os
import argparse
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

input = sys.argv[1]
output = sys.argv[2]

df = pd.read_csv(input, sep=';') #, dtype={'n_chiral': int})

def filter(row):
	if row['n_chiral'] > CHIRAL_THRESHOLD:
		return 0
	if row['psa'] > PSA_THRESHOLD:
		return 0
	if row['pfi'] > PFI_THRESHOLD:
		return 0
	if row['n_rot_bonds'] > ROTBOND_THRESHOLD:
		return 0
	if row['mw'] > WEIGHT_THRESHOLD:
		return 0
	if row['h_don'] > H_DON_THRESHOLD:
		return 0
	if row['h_acc'] > H_ACC_THRESHOLD:
		return 0
	if row['logp'] > LOGP_THRESHOLD_UP:
		return 0
	if row['logp'] < LOGP_THRESHOLD_LOW:
		return 0
	return 1

if 'filter_pass' in df.columns:
	df['filter_pass_original'] = df['filter_pass']	

df['filter_pass'] = df.apply (lambda row: filter(row), axis=1)

df.sort_values(['filter_pass','mpo'], ascending=[False, True], inplace=True)

total_pass = df['filter_pass'].sum()
print('Total_pass = ', total_pass)

df.to_csv(output, sep=';', index=False)