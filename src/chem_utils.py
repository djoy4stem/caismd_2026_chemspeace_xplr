"""
chem_utils.py — Workshop utility functions
===========================================
Low-level infrastructure helpers used by the workshop notebook.
These are NOT teaching cells — they handle format conversions and
correctness checks that would clutter the notebook without adding
pedagogical value.

Imported in Cell A1 as:
    from chem_utils import fragment_count, numpy_to_rdkit_fp
"""

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdmolops


def fragment_count(mol):
    """
    Return the number of disconnected fragments in a parsed Mol object.

    Uses RDKit's GetMolFrags() — accurate regardless of SMILES notation.
    A '.' in a SMILES string is NOT a reliable proxy for fragment count,
    because some single-fragment SMILES (e.g. certain ring closures or
    isotope/charge notations) can also contain dots.

    Parameters
    ----------
    mol : RDKit Mol or None

    Returns
    -------
    int : number of fragments (0 if mol is None)
    """
    if mol is None:
        return 0
    return len(rdmolops.GetMolFrags(mol))


def numpy_to_rdkit_fp(arr):
    """
    Convert a numpy uint8 bit array back to an RDKit ExplicitBitVect.

    This conversion is required by RDKit's BulkTanimotoSimilarity(),
    which does not accept numpy arrays directly.

    Parameters
    ----------
    arr : numpy.ndarray, shape (nbits,), dtype uint8
        Binary fingerprint as produced by mol_to_fp().

    Returns
    -------
    rdkit.DataStructs.ExplicitBitVect
    """
    fp = DataStructs.ExplicitBitVect(len(arr))
    fp.SetBitsFromList(np.where(arr == 1)[0].tolist())
    return fp
