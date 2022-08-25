import sys,os
import pandas as pd
import numpy as np
import sys
from datetime import datetime
from rdkit import Chem
from PP_ML_models.predictive_models.ml_model_gcnn_ens import Ensemble_Model_DC

import openeye.oechem as oe
from openeye import oemolprop as mp

from os.path import expanduser
home = expanduser("~")

# Add path so the predictive_models and properties modules can be found
head_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(head_path)
sys.path.append(home + '/PP_ML_models')

# Calculate logP:
def oeLogP(smi):
    mol = oe.OEGraphMol()
    if not oe.OESmilesToMol(mol, smi):
        print('ERROR: {}'.format(smi))
    else:
        logp = mp.OEGetXLogP(mol, atomxlogps=None)
    return logp


pIC50_pred_model = Ensemble_Model_DC(home + '/PP_ML_models/pIC50.pk')
print(pIC50_pred_model.info)
print(pIC50_pred_model.version)
# Run prediction model once to initialise:
_ = pIC50_pred_model.predict('C')[0]

# Number of compounds to generate predictions for in one go:
n_compounds = 4000

# Section of csv file to generate predictions for:
start_row = int(sys.argv[1])
end_row = int(sys.argv[2])

df = pd.read_csv('pymolgen.csv').iloc[start_row:end_row]['smi']

if (end_row is None) or (end_row > start_row + df.shape[0]):
    end_row = start_row + df.shape[0]

out_file = 'pymolgen_predictions_'+str(start_row)+'-'+str(end_row)+'.csv'

out = open(out_file, 'w')
out.write('SMILES,pIC50_pred,MPO,logP,PFI\n')

start_time = datetime.now()

for i in range(0, end_row-start_row, n_compounds):
    smi = df.iloc[i:i+n_compounds].to_list()
    pIC50_pred = pIC50_pred_model.predict(smi)[0]

    logp = np.array([oeLogP(s) for s in smi])
    n_aromatic = np.array([Chem.rdMolDescriptors.CalcNumAromaticRings(Chem.MolFromSmiles(s)) for s in smi])
    pfi = n_aromatic + logp
    mpo = (-pIC50_pred)*(1/(1 + np.exp(pfi - 8)))

    for s, p, m, f, l, n in zip(smi, pIC50_pred, mpo, pfi, logp, n_aromatic):
        out.write('{},{},{},{},{},{}\n'.format(s, p, m, f, l, n))
    out.flush()

end_time = datetime.now()

out.close()

total_time = (end_time - start_time).total_seconds()
print('Total time: {} s (~{} s/molecule)'.format(total_time, total_time/(end_row - start_row)))