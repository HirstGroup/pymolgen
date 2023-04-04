#!/usr/bin/env python

import argparse
import numpy as np
import random
import os
import sys
import time

from pymolgen.canonicalise_tautomer import canonicalise_tautomer
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

#   CHIRAL_THRESHOLD:   Maximum number of Chiral Centers (<= 2)
CHIRAL_THRESHOLD = 2
#   PSA_THRESHOLD:       Polar surface area (<= 140)
PSA_THRESHOLD = 140
#   PFI_THRESHOLD:       Property Forecast Index (LOGP + # of aromatic rings < 8)
PFI_THRESHOLD = 8
#   ROTBOND_THRESHOLD:  Number of rotatable bonds (<= 7)
ROTBOND_THRESHOLD = 7


#   H_DON_THRESHOLD:     Maximum number of hydrogen donors (<= 5)
H_DON_THRESHOLD = 5
#   H_ACC_THRESHOLD:     Maximum number of hydrogen acceptors (<= 10)
H_ACC_THRESHOLD = 10
#   LOGP_THRESHOLD:      Water/Octanol Partition Coefficient (0.5-5.0)
LOGP_THRESHOLD_UP = 5
LOGP_THRESHOLD_LOW = 0.5

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

    parser = argparse.ArgumentParser(description='Pymolgen molecular generator from fragments')
    parser.add_argument('-i','--input', help='Input file of inchi',required=True)
    parser.add_argument('-o','--output', help='Output file',required=True)
    parser.add_argument('-bp', '--build_probability', action='store_true', default=False, help='Read build probability from input and save into output', required=False)
    parser.add_argument('--mw_threshold', type=float, help='MW threshold', default = 500.0, required=False)

    args = parser.parse_args()

    WEIGHT_THRESHOLD = args.mw_threshold

    pIC50_pred_model = Ensemble_Model_DC(home + '/PP_ML_models/pIC50.pk', tauto=False)
    print(pIC50_pred_model.info)
    print(pIC50_pred_model.version)
    # Run prediction model once to initialise:
    _ = pIC50_pred_model.predict('C')[0]

    with open(args.output  , 'w') as outfile:
        outfile.write('inchi;smi;mw;n_rot_bonds;n_chiral;h_acc;h_don;psa;logp;n_aromatic;pfi;pIC50_pred;mpo;filter_pass')
        if args.build_probability is True:
            outfile.write(';build_probability')
        outfile.write('\n')
        print('Writing to', args.output)

    infile = open(args.input)

    for line in infile:

        inchi = line.split()[0].strip('\n')

        mw, n_rot_bonds, n_chiral, h_acc, h_don, psa, logp, n_aromatic, pfi, pIC50_pred, mpo = '', '', '', '', '', '', '', '', '', '', ''

        filter_pass = None

        if args.build_probability is True:
            build_probability = round(np.log10(float(line.strip().split()[1])), 2)

        try:
            rdmod, smi, oemol = None, None, None

            rdmol = Chem.MolFromInchi(inchi)

            smi = Chem.MolToSmiles(rdmol)

            smi = canonicalise_tautomer(smi)

            filter_pass = True

            mw = round(Chem.Descriptors.MolWt(rdmol), 2)

            oemol = oechem.OEGraphMol()
            oechem.OESmilesToMol(oemol, smi)
            oechem.OEAddExplicitHydrogens(oemol)

            n_rot_bonds = num_rot_bond(oemol)

            if n_rot_bonds > ROTBOND_THRESHOLD:
                filter_pass = False

            n_chiral = num_chiral_centres(oemol)

            if n_chiral > CHIRAL_THRESHOLD:
                filter_pass = False

            h_acc = num_lipinsky_acceptors(oemol)

            if h_acc > H_ACC_THRESHOLD:
                filter_pass = False

            h_don = num_lipinsky_donors(oemol)

            if h_don > H_DON_THRESHOLD:
                filter_pass = False

            psa = round(mp.OEGet2dPSA(oemol,atomPSA = None), 2)

            if psa > PSA_THRESHOLD:
                filter_pass = False

            logp = round(mp.OEGetXLogP(oemol, atomxlogps=None), 2)

            if logp > LOGP_THRESHOLD_UP:
                filter_pass = False

            if logp < LOGP_THRESHOLD_LOW:
                filter_pass = False

            n_aromatic = Chem.rdMolDescriptors.CalcNumAromaticRings(Chem.MolFromSmiles(smi))

            pfi = round(n_aromatic + logp, 2)

            pIC50_pred = round(pIC50_pred_model.predict(smi)[0], 2)

            mpo = round((-pIC50_pred)*(1/(1 + np.exp(pfi - 8))), 2)

        except:
            print('Could not calculate properties for', inchi)

        with open(args.output, 'a') as out:
            out.write(f'{inchi};{smi};{mw};{n_rot_bonds};{n_chiral};{h_acc};{h_don};{psa};{logp};{n_aromatic};{pfi};{pIC50_pred};{mpo};{filter_pass}')
            if args.build_probability is True:
                out.write(f';{build_probability}')
            out.write('\n')

    sys.exit('Normal termination')

            

