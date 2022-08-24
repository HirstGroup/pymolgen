import sys,os
import argparse

from fragment_builder import build_molecule
from pymolgen.molecule_formats import molecule_to_smiles

def oracle(mol):

    mw = mol.molecular_weight()

    if mw < 500:
        return True

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Pymolgen molecular generator from fragments')
    parser.add_argument('-a','--fragments_sdf', help='SDF file of fragments',required=True)
    parser.add_argument('-f','--fragments_txt', help='List of fragments in TXT file',required=True)
    parser.add_argument('-d','--frequencies_txt', help='Bond frequencies dictionary in txt file',required=True)
    parser.add_argument('-p','--parent_file', help='Parent Structure File in SDF format',required=True)
    parser.add_argument('-x','--parent_fragment_file_list', nargs='+', help='Parent Fragment Structure File list space-separated to search fragment database in SDF format',required=True)
    parser.add_argument('--parent_mapping_1', nargs='+', type=int, help='Parent Fragment i dict list space-separated to search fragment database in SDF format',required=True)
    parser.add_argument('--dict', nargs='+', type=int, help='Parent Fragment i dict list space-separated to search fragment database in SDF format',required=True)
    parser.add_argument('-r','--remove_hydrogens', type=int, nargs='+', help='Space-separated hydrogen atoms that will be created as attachment points, numbered from 0',required=False)
    parser.add_argument('-R','--remove_hydrogens_parent_fragment', type=int, nargs='+', help='Space-separated hydrogen atoms that will be created as attachment points for the parent fragment in database, numbered from 0',required=True)
    parser.add_argument('-s','--seed', type=int, help='Seed for random number generator',required=False)
    parser.add_argument('-o','--outfile_name', help='Output File Name',required=True)
    parser.add_argument('-n','--n_mol', type=int, help='Number of molecules to generate',required=True)
    parser.add_argument('--unique', action='store_true', help='Generate unique set of molecules', required=False)
    parser.add_argument('--rules', action='store_true', help='Use rules to filter', required=False)
    parser.add_argument('--rules_file', help='Rules file name for rules to filter', required=False)
    parser.add_argument('--filters', action='store_true', help='Use filters', required=False)
    parser.add_argument('--restart', action='store_true', help='Restart generation from previous run')
    parser.add_argument('--verbose', action='store_true', help='Verbose output', required=False)
    parser.add_argument('--mw_check', action='store_true', help='MW filter in every fragment addition')
    parser.add_argument('--no_numpy', action='store_true', help='Do not use numpy for fragment bond frequencies')
    parser.add_argument('--batch_size', type=int, help='Batch size for rules')

    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    if args.unique:
        print('Unique not fully working since does not take symmetry into account')

    use_numpy = not args.no_numpy

    n = 0

    for mol_list in build_molecule(fragments_sdf=args.fragments_sdf, fragments_txt=args.fragments_txt, frequencies_txt=args.frequencies_txt, parent_file=args.parent_file, parent_fragment_file_list=args.parent_fragment_file_list, parent_mapping_1=args.parent_mapping_1, parent_fragment_i_dict=args.dict, remove_hydrogens=args.remove_hydrogens, remove_hydrogens_parent_fragment=args.remove_hydrogens_parent_fragment,outfile_name=args.outfile_name, unique=args.unique, rules=args.rules, rules_file=args.rules_file, filters=args.filters, restart=args.restart, verbose=args.verbose, use_numpy=use_numpy, batch_size=args.batch_size):

        for mol in mol_list:

            smi = molecule_to_smiles(mol)
            mw = mol.molecular_weight()

            if oracle(mol):

                n += 1

                print(n, mw, smi)

                if n == args.n_mol:
                    sys.exit('Normal termination')

            

