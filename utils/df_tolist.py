import sys,os

import pandas as pd

infile = sys.argv[1]

df = pd.read_csv(infile, sep=';')

print(df.head())

