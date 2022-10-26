import pandas as pd
import numpy as np
import sys,os 
import argparse

parser = argparse.ArgumentParser(description='Combine calculator results into a single table')
parser.add_argument('-i','--input', nargs='+', help='Space-separated list of input files',required=True)
parser.add_argument('-o','--output', help='Output files',required=True)
args = parser.parse_args()

df = pd.read_csv(args.input[0], sep=';')

print(len(df))

for i in args.input[1:]:

    print(i, len(df))

    df2 = pd.read_csv(i, sep=';')

    print(len(df2))
    
    df = pd.concat([df, df2], ignore_index=True)

print(len(df))

df.to_csv(args.output, sep=';', index=False)


