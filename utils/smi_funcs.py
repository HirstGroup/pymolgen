# Writen October 2020

"""
Functions for modifying SMILES and calculating descriptors.
"""

__version__ = '2021.3.1'


from rdkit import Chem
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

# Need openbabel/openeye for pH conversions:
from openbabel import openbabel as ob
import openeye as oe
from openeye import oechem

from io import StringIO
import sys
import os
import numpy as np
import re


# Canonicalise SMILES:
def canonicalise_smiles(smi, method='rdkit'):
    """
    Convert SMILES to canonical RDKit form.
    """
    canon_smi = Chem.MolToSmiles(Chem.MolFromSmiles(smi), 
                                 canonical=True,
                                 isomericSmiles=True)
    return canon_smi


# Canonicalise tautomer:
def canonicalise_tautomer(smi, method='rdkit'):
    """
    Convert a SMILES to a canonical tautomer, 
    using RDKit, OpenBabel or via InChIs.
    """

    # Access any rdkit errors from stderr:
    Chem.WrapLogs()
    sio = sys.stderr = StringIO()

    if method == 'rdkit':
        enumerator = rdMolStandardize.TautomerEnumerator()
        # enumerator.SetReassignStereo = True
        # enumerator.SetRemoveBondStereo = False
        # enumerator.SetRemoveSp3Stereo = False
        mol = enumerator.Canonicalize(Chem.MolFromSmiles(smi))
        canon_smi = Chem.MolToSmiles(mol, canonical=True, 
                                     isomericSmiles=True)

        # If no change in SMILES, revert to original SMILES to retain
        # stereochemistry which may have been lost in 
        # TautomerEnumerator():
        canon_smi_2D = Chem.MolToSmiles(mol, canonical=True, 
                                        isomericSmiles=False)
        smi_2D = Chem.MolToSmiles(Chem.MolFromSmiles(smi),
                                  canonical=True, isomericSmiles=False)
        if canon_smi != smi and canon_smi_2D == smi_2D:
            canon_smi = smi

    elif method == 'inchi':
        # Convert to and from InChI to standarise tautomer:
        mol_smi = Chem.MolFromSmiles(smi)
        inchi = Chem.MolToInchi(mol_smi)
        mol_inchi = Chem.MolFromInchi(inchi)
        canon_smi = Chem.MolToSmiles(mol_inchi, canonical=True, 
                                     isomericSmiles=True)

    elif method == 'obabel':
        # Could use otautomer to get obabel canonical tautomers
        print('Warning: Method not yet implemented', file=sys.stderr)
        canon_smi = smi

    elif method == 'openeye':
        # Could use openeye to get a canonical tautomer
        print('Warning: Method not yet implemented', file=sys.stderr)
        canon_smi = smi

    # Get any warnings in stderr:
    warnings = sio.getvalue()

    # Redirect errors back to stderr:
    sys.stderr = sys.__stderr__

    return canon_smi, warnings


# Correct SMILES using regex:
def correct_smiles(smi, smi_transforms, check_rdkit=True):
    """
    Modify SMILES using regex.
    """

    # Access any rdkit errors from stderr:
    Chem.WrapLogs()
    sio = sys.stderr = StringIO()

    for match_pat, repl_pat in smi_transforms:
        corr_smi = re.sub(match_pat, repl_pat, smi)

    # Check that new SMILES can be loaded into RDKit:
    if check_rdkit == True and Chem.MolFromSmiles(corr_smi) is None:
        warnings = sio.getvalue() + \
        'WARNING (correct_smiles): '
        'Error during tautomer correction, '
        'will use original SMILES.\n'
        corr_smi = smi
    else:
        warnings = sio.getvalue()

    # Redirect errors back to stderr:
    sys.stderr = sys.__stderr__

    return corr_smi, warnings


# Protonate/deprotonate to get SMILES at a given pH:
def adjust_for_ph(smi, 
                  ph=7.4, 
                  phmodel='OpenEye', 
                  phmodel_dir=None, 
                  verbose=0 # False
                 ):
    """
    Protonate/deprotonate SMILES according
    to a given pH using OpenBabel or OpenEye.
    """

    if phmodel == 'OpenEye' and ph != 7.4:
        print('WARNING:')
        phmodel = 'OpenBabel'

    if verbose:
        print('Adjusting SMILES to pH: {}, using {}'.format(ph, phmodel))

    if phmodel == 'OpenBabel':

        # Access any rdkit errors from stderr:
        Chem.WrapLogs()
        sio = sys.stderr = StringIO()

        # Convert smiles


        # Set BABEL_DATADIR environment variable if using a modified 
        # phmodel.txt file containing pH transformations not in the
        # current directory
        if phmodel_dir is not None:
            os.environ['BABEL_DATADIR'] = phmodel_dir

        # Use obabel from python to do pH correction:
        ob_conv = ob.OBConversion()
        ob_mol = ob.OBMol()
        ob_conv.SetInAndOutFormats("smi", "smi")
        ob_conv.ReadString(ob_mol, smi)
        ob_mol.AddHydrogens(False,  # only add polar H (i.e. not to C atoms)
                            True,   # correct for pH
                            ph)     # pH value
        ph_smi = ob_conv.WriteString(ob_mol,
                                     True)  # trimWhitespace

        # Check that pH adjusted SMILES can be read by RDKit, 
        # if not return original SMILES and a warning:
        if Chem.MolFromSmiles(ph_smi) is None:
            warnings = sio.getvalue() + \
                'WARNING (adjust_for_ph): ' \
                'Error during pH correction with obabel, ' \
                'will use original SMILES.\n'
            ph_smi = smi
        else:
            warnings = sio.getvalue()

        # Convert to canonical smi?

        # Redirect errors back to stderr:
        sys.stderr = sys.__stderr__

    elif phmodel == 'OpenEye':

        # Save any OE errors:
        warnos = oechem.oeosstream()
        oechem.OEThrow.SetOutputStream(warnos)

        # Need to add OE error handling, OEThrow?
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        oe.OESetNeutralpHModel(mol)
        ph_smi = oechem.OEMolToSmiles(mol)

        # Get warnings from output stream, gives a
        # byte-string, so also have to decode to Unicode:
        warnings = str(warnos.str(), 'UTF-8')

        # Redirect OE errors back to stderr:
        oechem.OEThrow.SetOutputStream(oe.oeerr)

    return ph_smi, warnings

# Calculate rdkit descriptors:
def calc_rdkit_descs(smi, desc_ls):
    """
    Calculate specific rdkit descriptors from SMILES.
    """

    # Access any rdkit errors from stderr:
    Chem.WrapLogs()
    sio = sys.stderr = StringIO()

    mol = Chem.MolFromSmiles(smi)
    # Check SMILES has been read by rdkit:
    if mol is None:
        raise Exception('ERROR (calc_rdkit_desc): '
            'Cannot read SMILES using RDKit.')

    # Calculate descriptors in desc_ls:
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(desc_ls)
    descs = np.zeros(len(desc_ls))
    descs = np.array(calc.CalcDescriptors(mol))

    # Calculate Ipc descriptor separately to ensure average value is
    # calculated, this is not default in RDKit-2020.09.1, but otherwise 
    # leads to very large values which can cause problems:
    if 'Ipc' in desc_ls:
        descs[desc_ls.index('Ipc')] = Chem.GraphDescriptors.Ipc(mol, 
                                                                avg=True)

    warnings = sio.getvalue()

    # Redirect errors back to stderr:
    sys.stderr = sys.__stderr__

    return descs, warnings