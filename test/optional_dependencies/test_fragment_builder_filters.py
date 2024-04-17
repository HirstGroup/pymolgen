import sys,os
import random

from pymolgen.molecule_formats import *
from pymolgen.molecule_visualization import *
from pymolgen.molecule import *
from pymolgen.fragment_mol import *
from pymolgen.fragment_builder import *
from rdkit import Chem

def test_filters_additive_mol():

    from openeye import oechem
    from pymolgen.newmol import filters_additive, gen_pains_database, filters_final
    pains_database = gen_pains_database()
    smi = 'COCCN1C(=O)C2(CN(C3C4=C(C=CC(C#N)=C4)OC(C)(C)C3O)C(=O)O2)C2=CC=C(C3=C(C)ON=C3C)C=C21'

    oemol = oechem.OEGraphMol()
    oechem.OESmilesToMol(oemol, smi)

    oechem.OEAddExplicitHydrogens(oemol)

    print(filters_additive(oemol, smi))
    print(filters_final(oemol, smi, pains_database))

def test_filters_final_mol():

    mol = molecule_from_sdf('../datasets/sdf/mol-1.sdf')
    from pymolgen.newmol import filters_additive, gen_pains_database, filters_final_mol
    pains_database = gen_pains_database()
    print(filters_final_mol(pains_database, mol))