import pandas as pd
import numpy as np
import sys,os 
import argparse

parser = argparse.ArgumentParser(description='Pymolgen molecular generator from fragments')
parser.add_argument('-i','--input', help='Input file of inchi',required=True)
parser.add_argument('-o','--output', help='Output file',required=True)
parser.add_argument('--mpo', type=float, help='MPO threshold',required=False)

args = parser.parse_args()

print(args.input, args.output)

df = pd.read_csv(args.input, sep=';')

df = df.loc[df['filter_pass'] == True]
df = df.loc[df['rules_filter'] == True]

df.sort_values('mpo', inplace=True)

if args.mpo is not None:
	df = df.loc[df['mpo'] <= args.mpo]

df.to_csv(args.output, sep=';', index=False)
