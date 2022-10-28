import pandas as pd
import numpy as np
import sys,os 

input = sys.argv[1]
output = sys.argv[2]

print(input, output)

df = pd.read_csv(input, sep=';')

df = df.loc[df['rules_filter'] == True]

df.sort_values('mpo', inplace=True)

df.to_csv(output, sep=';', index=False)
