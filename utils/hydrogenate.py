import sys,os
import numpy as np

import networkx

from pymolgen.fragment_mol import *
from pymolgen.fragment_builder import *

from functools import partial
print = partial(print, flush=True)

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def findfragment(fragments_sdf_in, fragments_txt_in, frequencies_txt_in, frag_frequencies_txt_in, output):

    fragment_database_mol = get_fragment_database('%s.sdf' %(fragments_sdf_in) )

    frequencies = get_bond_frequencies('%s.txt' %(frequencies_txt_in) )

    frag_frequencies = get_frag_frequencies('%s.txt' %(frag_frequencies_txt_in))

    frag_mapping = get_frag_mapping('%s.txt' %(fragments_txt_in) )

    fragment_database = []

    for i in fragment_database_mol:
        fragment_database.append(i.graph)

    frequencies = update_bond_frequencies(frequencies, frag_mapping)

    with open(output,'w') as outfile:

        for i in range(len(fragment_database_mol)):

            mol = fragment_database_mol[i]

            mol.hydrogenate()

            lines = molecule_to_sdf(mol)

            for line in lines:
                outfile.write(line)

            outfile.write('$$$$\n')


if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Combine fragmented molecules')
    parser.add_argument('-i','--in_sub', help='Input subscript',required=True)
    parser.add_argument('-o','--output', help='Output SDF file',required=True)

    args = parser.parse_args()

    in_sub = args.in_sub

    fragments_sdf_in = 'fragments%s' %in_sub
    fragments_txt_in = 'fragments%s' %in_sub
    frequencies_txt_in = 'frequencies%s' %in_sub
    frag_frequencies_txt_in = 'frag_frequencies%s' %in_sub

    findfragment(fragments_sdf_in, fragments_txt_in, frequencies_txt_in, frag_frequencies_txt_in, args.output)

    print('Normal termination')