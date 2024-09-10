import sys,os
import argparse
import numpy as np
import random
import time

from fragment_builder import build_molecule
from pymolgen.molecule_formats import molecule_to_smiles, molecule_to_inchi, molecule_to_sdf
from pymolgen.properties_pymolgen import *

from rdkit import Chem
from PP_ML_models.predictive_models.ml_model_gcnn import Ensemble_Model_DC

# Import Openeye Modules
from openeye import oechem
from openeye import oemolprop as mp

from functools import partial
print = partial(print, flush=True)

from os.path import expanduser
home = expanduser("~")

# Add path so the predictive_models and properties modules can be found
head_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(head_path)
sys.path.append(home + '/PP_ML_models')

def oracle(mol):

    mw = mol.molecular_weight()

    if mw < 500:
        return True

# Calculate logP:
def oeLogP(smi):
    mol = oechem.OEGraphMol()
    if not oechem.OESmilesToMol(mol, smi):
        print('ERROR: {}'.format(smi))
    else:
        logp = mp.OEGetXLogP(mol, atomxlogps=None)
    return logp

def oeLogP_oemol(oemol):
    mol = oechem.OEGraphMol()
    if not oechem.OESmilesToMol(mol, smi):
        print('ERROR: {}'.format(smi))
    else:
        logp = mp.OEGetXLogP(mol, atomxlogps=None)
    return logp

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Pymolgen molecular generator from fragments and evaluation of molecular properties')

    # required arguments
    parser.add_argument('-a','--fragments_sdf', help='SDF file of fragments', required=True)
    parser.add_argument('-d','--frequencies_txt', help='Bond frequencies dictionary in txt file', required=True)
    parser.add_argument('-f','--fragments_txt', help='List of fragments in TXT file', required=True)
    parser.add_argument('-l','--log', help='Log file name', required=True)
    parser.add_argument('-p','--parent_file', help='Parent Structure File in SDF format', required=True)
    parser.add_argument('-R','--remove_hydrogens_parent_fragment', type=int, nargs='+', help='Space-separated hydrogen atoms that will be created as attachment points for the parent fragment in database, numbered from 0', required=True)
    parser.add_argument('-x','--parent_fragment_file_list', nargs='+', help='Parent Fragment Structure File list space-separated to search fragment database in SDF format', required=True)
    parser.add_argument('--mpo', type=float, help='MPO threshold', required=True)
    parser.add_argument('--parent_mapping_1', nargs='+', type=int, help='Parent Fragment i dict list space-separated to search fragment database in SDF format', required=True)

    # optional arguments
    parser.add_argument('-n','--n_mol', type=int, help='Number of molecules to generate', required=True)
    parser.add_argument('-o','--outfile_name', help='Output File Name', required=True)
    parser.add_argument('-r','--remove_hydrogens', type=int, nargs='+', help='Space-separated hydrogen atoms that will be created as attachment points, numbered from 0', required=False)
    parser.add_argument('-s','--seed', type=int, help='Seed for random number generator', required=False)
    parser.add_argument('--batch_size', type=int, help='Batch size for rules', required=False)
    parser.add_argument('--filters', action='store_true', help='Use filters', required=False)
    parser.add_argument('--mw_check', action='store_true', help='MW filter in every fragment addition', required=False)
    parser.add_argument('--no_numpy', action='store_true', help='Do not use numpy for fragment bond frequencies', required=False)
    parser.add_argument('--no_time', action='store_true', help='Do not print time interval in csv file', required=False)
    parser.add_argument('--rules', action='store_true', help='Use rules to filter', required=False)
    parser.add_argument('--rules_file', help='Rules file name for rules to filter', required=False)
    parser.add_argument('--restart', action='store_true', help='Restart generation from previous run', required=False)
    parser.add_argument('--unique', action='store_true', help='Generate unique set of molecules', required=False)
    parser.add_argument('--verbose', action='store_true', help='Verbose output', required=False)

    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)
        if args.n_mol > 1000:
            sys.exit('Cannot run with seed and n_mol > 1000')

    use_numpy = not args.no_numpy

    pIC50_pred_model = Ensemble_Model_DC(home + '/PP_ML_models/pIC50.pk')
    print(pIC50_pred_model.info)
    print(pIC50_pred_model.version)
    # Run prediction model once to initialise:
    _ = pIC50_pred_model.predict('C')[0]

    with open(args.log, 'w') as out:
        if args.no_time:
            out.write('inchi;smi;pIC50_pred;mpo;pfi;psa;logp;n_aromatic;n_rot_bonds;n_chiral;h_acc;mw;h_don\n')
        else:
            out.write('inchi;smi;pIC50_pred;mpo;pfi;psa;logp;n_aromatic;n_rot_bonds;n_chiral;h_acc;mw;h_don;interval_time\n')

    with open(args.outfile_name, 'w') as outfile:
        print(f'Writing to {args.outfile_name}')

    start_time = time.time()
    current_time = start_time

    n = 0

    for mol_list in build_molecule(fragments_sdf=args.fragments_sdf, fragments_txt=args.fragments_txt, frequencies_txt=args.frequencies_txt, parent_file=args.parent_file, parent_fragment_file_list=args.parent_fragment_file_list, parent_mapping_1=args.parent_mapping_1, remove_hydrogens=args.remove_hydrogens, remove_hydrogens_parent_fragment=args.remove_hydrogens_parent_fragment, unique=args.unique, rules=args.rules, rules_file=args.rules_file, filters=args.filters, restart=args.restart, verbose=args.verbose, mw_check=args.mw_check, use_numpy=use_numpy, batch_size=args.batch_size):

        for mol in mol_list:

            smi = molecule_to_smiles(mol)
            mw = mol.molecular_weight()

            try:
                pIC50_pred = pIC50_pred_model.predict(smi)[0]

                oemol = oechem.OEGraphMol()
                oechem.OESmilesToMol(oemol, smi)
                oechem.OEAddExplicitHydrogens(oemol)

                logp = mp.OEGetXLogP(oemol, atomxlogps=None)

                n_aromatic = Chem.rdMolDescriptors.CalcNumAromaticRings(Chem.MolFromSmiles(smi))
            except:
                print('Could not calculate properties for', smi)
                continue

            pfi = n_aromatic + logp
            mpo = (-pIC50_pred)*(1/(1 + np.exp(pfi - 8)))

            inchi = molecule_to_inchi(mol)

            if mpo < args.mpo:

                n += 1

                previous_time = current_time
                current_time = time.time() - start_time
                interval_time = current_time - previous_time

                print(f'EVALUATOR {n} {mw:.2f} {inchi} {interval_time:.2f}')

                psa = mp.OEGet2dPSA(oemol,atomPSA = None)

                n_rot_bonds = num_rot_bond(oemol)

                n_chiral = num_chiral_centres(oemol)

                h_acc = num_lipinsky_acceptors(oemol)

                h_don = num_lipinsky_donors(oemol)

                with open(args.log, 'a') as out:
                    if args.no_time:
                        out.write(f'{inchi};{smi};{pIC50_pred};{mpo};{pfi};{psa};{logp};{n_aromatic};{n_rot_bonds};{n_chiral};{h_acc};{mw:.2f};{h_don}\n')
                    else:
                        out.write(f'{inchi};{smi};{pIC50_pred};{mpo};{pfi};{psa};{logp};{n_aromatic};{n_rot_bonds};{n_chiral};{h_acc};{mw:.2f};{h_don};{interval_time:.2f}\n')


                lines = molecule_to_sdf(mol)

                with open(args.outfile_name, 'a') as outfile:
                    for line in lines:
                        outfile.write(line)
                    outfile.write('$$$$\n')

                if n == args.n_mol:
                    sys.exit('Normal termination')

            

