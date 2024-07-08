import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import sys

#chunksize = 100000
#tfr = pd.read_csv('all_rules_pass_classify.dat', sep=';', chunksize=chunksize, iterator=True)
#df = pd.concat(tfr, ignore_index=True)

#df = pd.read_csv('all_rules_pass_classify_10k.dat', sep=';')
df = pd.read_csv('all_chembl_30_10k_calculator_ok.csv', sep=';')

print(df)

"""
# 1. Plot MW histogram

fig1, ax1 = plt.subplots()

ax1.hist(df['mw'], density =True, bins=100)  # density=False would make counts
#plt.title('MW')
ax1.set_ylabel('Frequency')
ax1.set_xlabel('MW')

fig1.savefig('mw.pdf') #, dpi=600)

print(list(df.columns.values))


# 2. Plot n_rot_bonds histogram

fig2, ax2 = plt.subplots()

max = df['n_rot_bonds'].max()
min = df['n_rot_bonds'].min()
n_bins = int(max - min + 1)

ax2.hist(df['n_rot_bonds'], density=True, bins=n_bins)  # density=False would make counts
#plt.title('MW')
ax2.set_ylabel('Frequency')
ax2.set_xlabel('n_rot_bonds')

fig2.savefig('n_rot_bonds.pdf') #, dpi=600)

# 3. Plot h_acc histogram

fig2, ax2 = plt.subplots()

max = df['h_acc'].max()
min = df['h_acc'].min()
n_bins = int(max - min + 1)

ax2.hist(df['h_acc'], density=True, bins=n_bins)  # density=False would make counts
#plt.title('MW')
ax2.set_ylabel('Frequency')
ax2.set_xlabel('H-bond acceptors')

fig2.savefig('h_acc.pdf') #, dpi=600)
"""

for prop in ['n_rot_bonds', 'n_chiral', 'h_acc', 'h_don', 'n_aromatic']:

	print('PROP', prop)

	fig, ax = plt.subplots()

	max = df[prop].max()
	min = df[prop].min()
	n_bins = int(max - min)

	a = ax.hist(df[prop], density=True, bins=n_bins)  # density=False would make counts
	pd.DataFrame({'x_upper':a[1][1:], 'y': a[0]}).to_csv(f'{prop}.csv')
	#plt.title('MW')
	ax.set_ylabel('Frequency')
	ax.set_xlabel(prop)

	fig.savefig(f'{prop}.pdf') #, dpi=600)	

#, 'psa', 'logp', 'n_aromatic', 'pfi', 'pIC50_pred', 'mpo', 'filter_pass', 'index', 'rules_filter', 'mol_type']
# ['inchi', 'smi', 'mw', 'n_rot_bonds', 'n_chiral', 'h_acc', 'h_don', 'psa', 'logp', 'n_aromatic', 'pfi', 'pIC50_pred', 'mpo', 'filter_pass', 'index', 'rules_filter', 'mol_type']

for prop in ['mw', 'psa', 'logp', 'pfi', 'pIC50_pred', 'mpo']:

	fig, ax = plt.subplots()

	a = ax.hist(df[prop], density=True, bins=100)  # density=False would make counts
	pd.DataFrame({'x_upper':a[1][1:], 'y': a[0]}).to_csv(f'{prop}.csv')
	#plt.title('MW')
	ax.set_ylabel('Frequency')
	ax.set_xlabel(prop)

	fig.savefig(f'{prop}.pdf') #, dpi=600)	 