"""
ECFP Feature Module
===================
Computes Extended Connectivity Fingerprints (ECFP4) and physicochemical descriptors.

References:
- Morgan Fingerprints: Rogers, D., & Hahn, M. (2010). Extended-connectivity fingerprints.
  J. Chem. Inf. Model., 50(5), 742-754. https://doi.org/10.1021/ci100050t
- RDKit: Landrum, G. (2006). RDKit: Open-source cheminformatics.
"""
import numpy as np
from typing import List, Dict, Optional

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, Lipinski
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

from ..config import ECFP_RADIUS, ECFP_BITS


def compute_ecfp(smiles: str, radius: int = ECFP_RADIUS, bits: int = ECFP_BITS) -> List[int]:
    """
    Compute ECFP (Morgan) fingerprint from SMILES string.
    
    Args:
        smiles: Canonical SMILES string.
        radius: Fingerprint radius (default 2 for ECFP4).
        bits: Fingerprint length (default 1024).
    
    Returns:
        List of bit values (0/1). Returns NaN-filled list on failure.
    """
    if not RDKIT_AVAILABLE:
        raise ImportError("RDKit is required for ECFP computation.")
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=bits)
            return list(fp)
    except Exception:
        pass
    return [np.nan] * bits


def compute_physchem(smiles: str) -> Dict[str, float]:
    """
    Compute physicochemical properties from SMILES.
    
    Returns:
        Dict with keys: logp, mw, tpsa, hbd, hba.
    """
    if not RDKIT_AVAILABLE:
        raise ImportError("RDKit is required for PhysChem computation.")
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return {
                'logp': Descriptors.MolLogP(mol),
                'mw': Descriptors.MolWt(mol),
                'tpsa': Descriptors.TPSA(mol),
                'hbd': Lipinski.NumHDonors(mol),
                'hba': Lipinski.NumHAcceptors(mol),
            }
    except Exception:
        pass
    return {k: np.nan for k in ['logp', 'mw', 'tpsa', 'hbd', 'hba']}


def canonicalize_smiles(smiles: str) -> Optional[str]:
    """
    Canonicalize a SMILES string using RDKit.
    
    Args:
        smiles: Input SMILES string.
    
    Returns:
        Canonical SMILES or None if parsing fails.
    """
    if not RDKIT_AVAILABLE:
        raise ImportError("RDKit is required for SMILES canonicalization.")
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        pass
    return None
