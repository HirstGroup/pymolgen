#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys,os 
import argparse

parser = argparse.ArgumentParser(description='Combine data in dataframe format to single file')
parser.add_argument('-i','--input', nargs='+', help='List of input file names',required=True)
parser.add_argument('-o','--output', help='Output file',required=True)

args = parser.parse_args()

print(args.input, args.output)

df = pd.read_csv(args.input[0], sep=';')

for i in args.input[1:]:

    print(i)

    df2 = pd.read_csv(i, sep=';')

    df = pd.concat([df, df2])

df.to_csv(args.output, sep=';', index=False)
