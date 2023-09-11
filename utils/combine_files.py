#!/usr/bin/env python

import sys,os 
import argparse

parser = argparse.ArgumentParser(description='Combine files containing header line')
parser.add_argument('-i','--input', nargs='+', help='List of input file names',required=True)
parser.add_argument('-o','--output', help='Output file',required=True)

args = parser.parse_args()

os.system(f'head -n1 {args.input[0]} > {args.output}')

for i in args.input:
	os.system(f'tail -n+2 {i} >> {args.output}')