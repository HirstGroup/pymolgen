import sys,os
import argparse
import pandas as pd

# Import Openeye Modules
from openeye import oechem
from openeye import oemolprop as mp

from pymolgen.properties_pymolgen import gen_pains_database

parser = argparse.ArgumentParser(description='Pymolgen molecular generator from fragments')
parser.add_argument('-i','--input', help='Input File Name',required=True)
parser.add_argument('-o','--output', help='Output File Name',required=True)

args = parser.parse_args()

def pains_filter(row):

    smi = row['smi']

    if row['filter_pass'] is False:
        return None

    oemol = oechem.OEGraphMol()
    oechem.OESmilesToMol(oemol, smi)

    oechem.OEAddExplicitHydrogens(oemol)    

    for fragment in pains_database:
        fragment_search = oechem.OESubSearch(fragment)
        oechem.OEPrepareSearch(oemol, fragment_search)
        if fragment_search.SingleMatch(oemol):
            return False
            break

    return True

try:
    pains_database = gen_pains_database()
except:
    raise Exception("Could not generate pains database")

df = pd.read_csv(args.input, sep=';')

print(df)

df['pains_filter'] = df.apply(pains_filter, axis=1)

print(df['pains_filter'])

df.to_csv(args.output, index=False, sep=';')