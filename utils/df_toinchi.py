import argparse
import sys, os
import pandas as pd


parser = argparse.ArgumentParser(description='Extract inchi column from csv file and make into new file')

parser.add_argument('-i','--input', help='Input csv file',required=True)
parser.add_argument('-o','--output', help='Output file with inchis',required=True)

args = parser.parse_args()

if args.input == args.output:
	sys.exit('Same input and output')

df = pd.read_csv(args.input, sep=';')

outfile = open(args.output, 'w')

for i in df['inchi']:
	outfile.write(f'{i}\n')

