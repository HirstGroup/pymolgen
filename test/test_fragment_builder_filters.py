import filecmp
import os
import random
import sys

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


def test_fragment_builder():
    # complete fragment_builder test including filters option

    cpu = 1

    batch_size = 10

    random.seed(100)

    os.system('rm outputs/*')

    fragment_builder(fragments_sdf='../datasets/database1000/fragments1000.sdf', fragments_txt='../datasets/database1000/fragments1000.txt', frequencies_txt='../datasets/database1000/frequencies1000.txt', parent_file='../datasets/database1000/phenylisoxazole.sdf', parent_fragment_file_list=['../datasets/database1000/benzene.sdf','../datasets/database1000/benzene.sdf'], parent_mapping_1=[16,0,15,0], remove_hydrogens=[20,21], remove_hydrogens_parent_fragment=[11,11],outfile_name='outputs/fragment_builder.sdf', n_mol=10, unique=True, rules=False, rules_file=None, filters=True, fragments_used_file='outputs/fragments_used.txt', restart=False, verbose=False, mw_check=True, use_numpy=True, batch_size=batch_size, cpu=cpu, candidate_file='outputs/candidates.txt', cap=True, intermediates=True)

    assert filecmp.cmp('models/candidates.txt', 'outputs/candidates.txt') is True
    assert filecmp.cmp('models/fragments_used.txt', 'outputs/fragments_used.txt') is True
    assert filecmp.cmp('models/fragment_builder.sdf', 'outputs/fragment_builder.sdf') is True