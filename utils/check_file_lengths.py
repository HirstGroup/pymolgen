#!/usr/bin/python

import sys,os
import argparse

 
parser = argparse.ArgumentParser(description='This script checks that two files have the same number of lines')
parser.add_argument('-a','--afile', help='File a to check', required=True)
parser.add_argument('-b','--bfile', help='File b to check', required=True)
parser.add_argument('-d','--difference', help='Difference in number of lines, file b - file a', default=0, type=int, required=False)

args = parser.parse_args()

def line_count(file_path):
    return int(os.popen(f'wc -l {file_path}').read().split()[0])

an = line_count(args.afile)
bn = line_count(args.bfile)

if bn - an - args.difference == 0:
	result = 'OK'
else:
	result = 'NOK'

print(args.afile, args.bfile, result)