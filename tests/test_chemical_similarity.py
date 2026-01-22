#!/usr/bin/env python3
"""
Tests for Chemical Similarity Analysis.

Tests fingerprint generation, Tanimoto similarity, caching, and results.
"""

import pytest
import json
import tempfile
from pathlib import Path
from unittest.mock import patch

# Skip if RDKit not available
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit import DataStructs
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False



@pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not installed")
class TestFingerprints:
    
    def test_fingerprint_from_valid_smiles(self):
        """Test fingerprint generation from valid SMILES."""
        # Quercetin SMILES
        smiles = "OC1=CC(O)=C2C(=O)C(O)=C(OC2=C1)C3=CC(O)=C(O)C=C3"
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None
        
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        assert fp is not None
        assert len(fp) == 2048
    
    def test_fingerprint_from_invalid_smiles(self):
        """Test fingerprint generation from invalid SMILES returns None."""
        invalid_smiles = "NOT_A_VALID_SMILES"
        mol = Chem.MolFromSmiles(invalid_smiles)
        assert mol is None
    
    def test_fingerprint_hyperforin(self):
        """Test Hyperforin fingerprint generation."""
        smiles = "CC(C)C(=O)C12C(=O)C(=C(C(C1=O)(CC(C2(C)CCC=C(C)C)CC=C(C)C)CC=C(C)C)O)CC=C(C)C"
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None
        
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        assert fp is not None



@pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not installed")
class TestTanimotoSimilarity:
    """Test Tanimoto similarity calculations."""
    
    def test_identical_molecules_similarity_1(self):
        """Identical molecules should have Tanimoto = 1.0."""
        smiles = "CCO"  # Ethanol
        mol = Chem.MolFromSmiles(smiles)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        
        sim = DataStructs.TanimotoSimilarity(fp, fp)
        assert sim == 1.0
    
    def test_different_molecules_similarity_less_than_1(self):
        """Different molecules should have Tanimoto < 1.0."""
        smiles1 = "CCO"  # Ethanol
        smiles2 = "CCCCCCCC"  # Octane
        
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
        
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius=2, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius=2, nBits=2048)
        
        sim = DataStructs.TanimotoSimilarity(fp1, fp2)
        assert 0 <= sim < 1.0
    
    def test_similar_molecules_higher_similarity(self):
        """Similar molecules should have higher Tanimoto than dissimilar ones."""
        # Ethanol vs Methanol (similar - both alcohols)
        ethanol = Chem.MolFromSmiles("CCO")
        methanol = Chem.MolFromSmiles("CO")
        # Ethanol vs Benzene (dissimilar)
        benzene = Chem.MolFromSmiles("c1ccccc1")
        
        fp_eth = AllChem.GetMorganFingerprintAsBitVect(ethanol, radius=2, nBits=2048)
        fp_meth = AllChem.GetMorganFingerprintAsBitVect(methanol, radius=2, nBits=2048)
        fp_benz = AllChem.GetMorganFingerprintAsBitVect(benzene, radius=2, nBits=2048)
        
        sim_similar = DataStructs.TanimotoSimilarity(fp_eth, fp_meth)
        sim_dissimilar = DataStructs.TanimotoSimilarity(fp_eth, fp_benz)
        
        assert sim_similar > sim_dissimilar
    
    def test_tanimoto_range(self):
        """Tanimoto similarity should be in [0, 1]."""
        smiles_list = ["CCO", "c1ccccc1", "CC(=O)O", "CCN"]
        fps = []
        
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
            fps.append(fp)
        
        for i, fp1 in enumerate(fps):
            for j, fp2 in enumerate(fps):
                sim = DataStructs.TanimotoSimilarity(fp1, fp2)
                assert 0 <= sim <= 1.0



class TestSMILESCache:
    """Test SMILES caching functionality."""
    
    def test_cache_save_and_load(self):
        """Test saving and loading SMILES cache."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_file = Path(tmpdir) / "test_cache.json"
            
            # Save cache
            test_cache = {
                "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
                "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            }
            
            with open(cache_file, 'w') as f:
                json.dump(test_cache, f)
            
            # Load cache
            with open(cache_file, 'r') as f:
                loaded_cache = json.load(f)
            
            assert loaded_cache == test_cache
    
    def test_cache_empty_values(self):
        """Test cache handles empty strings (missed lookups)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_file = Path(tmpdir) / "test_cache.json"
            
            # Cache with missed lookup
            test_cache = {
                "ValidDrug": "CCO",
                "InvalidDrug": "",  # Empty = miss
            }
            
            with open(cache_file, 'w') as f:
                json.dump(test_cache, f)
            
            with open(cache_file, 'r') as f:
                loaded_cache = json.load(f)
            
            assert loaded_cache["ValidDrug"] == "CCO"
            assert loaded_cache["InvalidDrug"] == ""



class TestResultsValidation:
    """Test that results files exist and have correct structure."""
    
    def test_summary_file_exists(self):
        """Test chemical similarity summary file exists."""
        summary_file = Path("results/tables/chemical_similarity_summary.csv")
        assert summary_file.exists(), "Chemical similarity summary not found"
    
    def test_summary_has_required_columns(self):
        """Test summary file has required columns."""
        import pandas as pd
        summary_file = Path("results/tables/chemical_similarity_summary.csv")
        
        if summary_file.exists():
            df = pd.read_csv(summary_file)
            required_cols = [
                'compound', 
                'max_sim_DILI_positive', 
                'max_sim_DILI_negative',
                'structural_analog_to_hepatotoxins'
            ]
            for col in required_cols:
                assert col in df.columns, f"Missing column: {col}"
    
    def test_hyperforin_not_analog(self):
        """Test Hyperforin is not a structural analog to hepatotoxins."""
        import pandas as pd
        summary_file = Path("results/tables/chemical_similarity_summary.csv")
        
        if summary_file.exists():
            df = pd.read_csv(summary_file)
            hyp = df[df['compound'] == 'Hyperforin']
            
            if len(hyp) > 0:
                assert hyp.iloc[0]['structural_analog_to_hepatotoxins'] == 'No'
                assert hyp.iloc[0]['max_sim_DILI_positive'] < 0.4
    
    def test_quercetin_not_analog(self):
        """Test Quercetin is not a structural analog to hepatotoxins."""
        import pandas as pd
        summary_file = Path("results/tables/chemical_similarity_summary.csv")
        
        if summary_file.exists():
            df = pd.read_csv(summary_file)
            quer = df[df['compound'] == 'Quercetin']
            
            if len(quer) > 0:
                assert quer.iloc[0]['structural_analog_to_hepatotoxins'] == 'No'
                assert quer.iloc[0]['max_sim_DILI_positive'] < 0.4
    
    def test_dilirank_reference_set_exists(self):
        """Test DILIrank reference set file exists."""
        ref_file = Path("results/tables/dilirank_reference_set.csv")
        assert ref_file.exists(), "DILIrank reference set not found"
    
    def test_dilirank_reference_counts(self):
        """Test DILIrank reference set has expected counts."""
        import pandas as pd
        ref_file = Path("results/tables/dilirank_reference_set.csv")
        
        if ref_file.exists():
            df = pd.read_csv(ref_file)
            
            dili_pos = len(df[df['category'] == 'DILI_positive'])
            dili_neg = len(df[df['category'] == 'DILI_negative'])
            
            # Expected counts from FDA DILIrank 2.0
            assert dili_pos >= 500, f"Expected ~542 DILI+, got {dili_pos}"
            assert dili_neg >= 300, f"Expected ~365 DILI-, got {dili_neg}"
