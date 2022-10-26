import pandas as pd
import numpy as np
import sys,os 

input = sys.argv[1]
output = sys.argv[2]

df = pd.read_csv(input, sep=';')

print(df.head())

df = df.loc[df['rules_filter'] == True]

df.sort_values('mpo', inplace=True)

df = df.loc[df['mpo'] <= -6.5]

df.to_csv(output, sep=';', index=False)
