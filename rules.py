import sys,os
import argparse
import pandas as pd
import subprocess

parser = argparse.ArgumentParser(description='Pymolgen molecular generator from fragments')
parser.add_argument('-i','--input', help='Input File Name',required=True)
parser.add_argument('-o','--output', help='Output File Name',required=True)
parser.add_argument('-r','--rules_file', help='Rules File Name',required=True)

args = parser.parse_args()

df = pd.read_csv(args.input, sep=';')

df['index'] = df.index

df2 = df.loc[df['filter_pass'] == True]

df2 = df2[['smi', 'index']]

with open(args.rules_file, 'w') as outfile:
    for i, row in df2.iterrows():
        outfile.write('%s %s\n' %(row['smi'], row['index']))

home = os.path.expanduser('~/')

result = subprocess.run([home + 'Lilly-Medchem-Rules/Lilly_Medchem_Rules.rb %s' %args.rules_file], shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8')

print(result)

rules_pass_list = set()

for line in result.split('\n'):
    if not line.strip():
        continue
    i_mol = int(line.split()[1])

    rules_pass_list.add(i_mol)

print(rules_pass_list)

def rules_filter(row):

    if row['index'] in rules_pass_list:
        return True

    return False

df['rules_filter'] = df.apply(rules_filter, axis=1)

print(df)

df.to_csv(args.output, index=False, sep=';')