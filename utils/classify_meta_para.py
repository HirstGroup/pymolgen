#!/usr/bin/env python

import argparse
import pandas as pd

from rdkit import Chem

para = Chem.MolFromSmarts("[cH]1[cH]c(-*)[cH][cH]c1-c2c(-[#6])onc2-[#6]")
meta = Chem.MolFromSmarts("[cH]1c(-*)[cH][cH][cH]c1-c2c(-[#6])onc2-[#6]")
parameta = Chem.MolFromSmarts("[cH]1c(-*)c(-*)[cH][cH]c1-c2c(-[#6])onc2-[#6]")


def classify1(mol):
    """
    Classify RDKit molecule into meta, para or meta/para phenylisoxazole

    Parameters
    ----------
    mol : RDKit molecule object

    Returns
    -------
    mol_type : str
        classification of molecule into meta, para or parameta phenylisoxazole

    """

    mol_type = 'NONE'

    if mol.HasSubstructMatch(para):
        mol_type = 'para'
    elif mol.HasSubstructMatch(meta):
        mol_type = 'meta'
    elif mol.HasSubstructMatch(parameta):
        mol_type = 'parameta'

    return mol_type


def classify(inchi):
    """
    Classify inchi into meta, para or meta/para phenylisoxazole

    Parameters
    ----------
    inchi : str
        Inchi string

    Returns
    -------
    mol_type : str
        classification of molecule into meta, para or parameta phenylisoxazole, NONE or error

    """

    mol_type = 'NONE'

    try:

        mol = Chem.inchi.MolFromInchi(inchi)

        if mol.HasSubstructMatch(para):
            mol_type = 'para'
        elif mol.HasSubstructMatch(meta):
            mol_type = 'meta'
        elif mol.HasSubstructMatch(parameta):
            mol_type = 'parameta'

    except:

        mol_type = 'error'

    return mol_type


def classify_row(row):
    """
    Run classify function for pandas row

    Parameters
    ----------
    row : pandas row
    """

    #mol = Chem.inchi.MolFromInchi(row['inchi'])

    mol_type = classify(row['inchi'])

    return mol_type


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Classify molecules as meta, para or meta/para phenylisoxazole')
    
    parser.add_argument('-i','--input', help='Input file',required=True)
    parser.add_argument('-fi','--format_input', help='Format of input file, either inchi or csv (separated by ;)',required=True)
    parser.add_argument('-o','--output', help='Output file',required=True)

    args = parser.parse_args()

    outfile = open(args.output, 'w')

    if args.format_input == 'inchi':

        infile = open(args.input)

        for line in infile:

            inchi = line.split()[0]

            mol_type = classify(inchi)

            outfile.write(f'{inchi} {mol_type}\n')

    elif args.format_input == 'csv':

        df = pd.read_csv(args.input, sep=';')

        df['mol_type'] = df.apply(classify_row, axis=1)

        df.to_csv(args.output, sep=';', index=False)

    print('Normal termination')

