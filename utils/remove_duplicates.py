import pandas as pd
import numpy as np
import sys,os 

df = pd.read_csv('results.csv')

inchikey = []

with open('results.inchikey') as infile:
    for line in infile:
        inchikey.append(line.strip('\n'))

df['inchikey'] = inchikey

df['Count'] = 1

df2 = df.groupby(['inchikey'])['Count'].count().reset_index()

df = df.drop_duplicates(subset='inchikey', keep="first")

df['Count'] = df2['Count'].tolist()

df.sort_values(by=['MPO'], inplace=True)

df.to_csv('results-single.csv', index=False)

print(df)
