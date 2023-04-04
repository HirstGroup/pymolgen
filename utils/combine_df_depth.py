#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys,os 
import argparse

parser = argparse.ArgumentParser(description='Combine data in dataframe format to single file')
parser.add_argument('-i','--input', help='Input file containing filenames and depth space separated',required=True)
parser.add_argument('-o','--output', help='Output file',required=True)
parser.add_argument('--sort', help='Name of column to sort by',required=False)
parser.add_argument('--reverse', action='store_true', default=False, help='Reverse sort', required=False)

args = parser.parse_args()

print(args.input, args.output)

file_list = []
depth_list = []

with open(args.input) as infile:
	for line in infile:
		file_list.append(line.split()[0])
		depth_list.append(line.split()[1])

print(file_list)
print(depth_list)

df = pd.read_csv(file_list[0], sep=';')

df['depth'] = depth_list[0]

for i in range(len(file_list)):

	df2 = pd.read_csv(file_list[i], sep=';')

	df2['depth'] = depth_list[i]

	df = pd.concat([df, df2])

if args.sort is not None:
	df.sort_values(args.sort, ascending=not args.reverse, inplace=True)

df.to_csv(args.output, sep=';', index=False)

print(df)
