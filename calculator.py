import sys,os
import argparse
import numpy as np
import random
import time

from pymolgen.molecule_formats import molecule_to_smiles, molecule_to_inchi, molecule_to_sdf
from pymolgen.properties_pymolgen import *

from rdkit import Chem
from PP_ML_models.predictive_models.ml_model_gcnn_ens import Ensemble_Model_DC

# Import Openeye Modules
from openeye import oechem
from openeye import oemolprop as mp

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

    parser = argparse.ArgumentParser(description='Pymolgen molecular generator from fragments')
    parser.add_argument('-i','--input', help='Input file of inchi',required=True)
    parser.add_argument('-o','--output', help='Output file',required=True)

    args = parser.parse_args()

    pIC50_pred_model = Ensemble_Model_DC(home + '/PP_ML_models/pIC50.pk')
    print(pIC50_pred_model.info)
    print(pIC50_pred_model.version)
    # Run prediction model once to initialise:
    _ = pIC50_pred_model.predict('C')[0]

    with open(args.output  , 'w') as outfile:
        print('Writing to', args.output)

    infile = open(args.input)

    for line in infile:

        inchi = line.strip('\n')
        print(inchi)

        try: 
            rdmol = Chem.MolFromInchi(inchi)
            smi = Chem.MolToSmiles(rdmol)

            pIC50_pred = pIC50_pred_model.predict(smi)[0]

            oemol = oechem.OEGraphMol()
            oechem.OESmilesToMol(oemol, smi)
            oechem.OEAddExplicitHydrogens(oemol)

            logp = mp.OEGetXLogP(oemol, atomxlogps=None)

            n_aromatic = Chem.rdMolDescriptors.CalcNumAromaticRings(Chem.MolFromSmiles(smi))


        except:
            print('Could not calculate properties for', inchi)
            continue

        pfi = n_aromatic + logp
        mpo = (-pIC50_pred)*(1/(1 + np.exp(pfi - 8)))

        n_rot_bonds = num_rot_bond(oemol)

        n_chiral = num_chiral_centres(oemol)

        h_acc = num_lipinsky_acceptors(oemol)

        h_don = num_lipinsky_donors(oemol)

        with open(args.output, 'a') as out:
            out.write(f'{inchi};{smi};{pIC50_pred};{mpo};{pfi};{psa};{logp};{n_aromatic};{n_rot_bonds};{n_chiral};{h_acc};{mw:.2f};{h_don}\n')

    sys.exit('Normal termination')

            

