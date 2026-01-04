"""
Unified Compound Identification System
======================================

Robust compound identification to prevent:
1. Name mismatches (cyclosporine vs ciclosporin)
2. Duplicate compounds (same drug, different names)
3. Wrong matching (similar names, different drugs)

Uses canonical identifiers: DrugBank ID > InChIKey > Normalized Name
"""
import pandas as pd
import numpy as np
from pathlib import Path
import hashlib
import re
import warnings
warnings.filterwarnings('ignore')

try:
    from rdkit import Chem
    from rdkit.Chem.inchi import MolFromInchi, MolToInchiKey
    RDKIT_AVAILABLE = True
except:
    RDKIT_AVAILABLE = False

ROOT = Path(r'v:\new\h-perforatum-network-tox')

# =============================================================================
# DRUG NAME SYNONYMS DATABASE
# =============================================================================

DRUG_SYNONYMS = {
    # International vs US naming
    'ciclosporin': 'cyclosporine',
    'ciclosporine': 'cyclosporine',
    'cyclosporin': 'cyclosporine',
    'cyclosporin a': 'cyclosporine',
    'paracetamol': 'acetaminophen',
    'rifampin': 'rifampicin',
    'epinephrine': 'adrenaline',
    'norepinephrine': 'noradrenaline',
    'furosemide': 'frusemide',
    'metolazone': 'metalazone',
    
    # Brand to generic
    'tylenol': 'acetaminophen',
    'advil': 'ibuprofen',
    'motrin': 'ibuprofen',
    'lipitor': 'atorvastatin',
    'zocor': 'simvastatin',
    'plavix': 'clopidogrel',
    'nexium': 'esomeprazole',
    'prilosec': 'omeprazole',
    'prozac': 'fluoxetine',
    'zoloft': 'sertraline',
    'lexapro': 'escitalopram',
    'cymbalta': 'duloxetine',
    'coumadin': 'warfarin',
    'glucophage': 'metformin',
    'synthroid': 'levothyroxine',
    'eltroxin': 'levothyroxine',
    'neoral': 'cyclosporine',
    'sandimmune': 'cyclosporine',
    'prograf': 'tacrolimus',
    'cellcept': 'mycophenolate',
    'imuran': 'azathioprine',
    'rheumatrex': 'methotrexate',
    'trexall': 'methotrexate',
    
    # Spelling variants
    'sulfasalazine': 'sulphasalazine',
    'sulfamethoxazole': 'sulphamethoxazole',
    'esomeprazole': 'esomeprazol',
    'omeprazole': 'omeprazol',
    'lansoprazole': 'lansoprazol',
    
    # Salt forms to base
    'metformin hydrochloride': 'metformin',
    'atorvastatin calcium': 'atorvastatin',
    'simvastatin sodium': 'simvastatin',
    'omeprazole magnesium': 'omeprazole',
    'esomeprazole magnesium': 'esomeprazole',
    'levothyroxine sodium': 'levothyroxine',
    'cyclosporine a': 'cyclosporine',
    
    # Common abbreviations
    'mtx': 'methotrexate',
    'aza': 'azathioprine',
    'csa': 'cyclosporine',
    'mpa': 'mycophenolic acid',
    '5-asa': 'mesalamine',
    '6-mp': 'mercaptopurine',
}

# =============================================================================
# COMPOUND IDENTIFIER CLASS
# =============================================================================

class CompoundIdentifier:
    """
    Unified compound identification system.
    
    Priority:
    1. DrugBank ID (most reliable)
    2. InChIKey (structure-based)
    3. Normalized name (after synonyms)
    4. Fuzzy match (last resort)
    """
    
    def __init__(self):
        self.synonyms = DRUG_SYNONYMS
        self.compound_registry = {}  # canonical_id -> compound info
        self.name_to_id = {}  # normalized name -> canonical_id
        self.smiles_to_id = {}  # canonical smiles -> canonical_id
        
    def normalize_name(self, name):
        """Normalize drug name to canonical form."""
        if pd.isna(name):
            return None
            
        # Lowercase and strip
        name = str(name).lower().strip()
        
        # Remove common suffixes
        name = re.sub(r'\s+(hydrochloride|hcl|sodium|potassium|calcium|magnesium|acetate|succinate)$', '', name)
        name = re.sub(r'\s+(monohydrate|dihydrate|trihydrate)$', '', name)
        
        # Apply synonym mapping
        if name in self.synonyms:
            name = self.synonyms[name]
            
        return name
    
    def get_canonical_smiles(self, smiles):
        """Get canonical SMILES from any SMILES string."""
        if not RDKIT_AVAILABLE or pd.isna(smiles):
            return None
            
        try:
            mol = Chem.MolFromSmiles(str(smiles))
            if mol:
                return Chem.MolToSmiles(mol, canonical=True)
        except:
            pass
        return None
    
    def get_inchikey(self, smiles):
        """Get InChIKey from SMILES (structure fingerprint)."""
        if not RDKIT_AVAILABLE or pd.isna(smiles):
            return None
            
        try:
            mol = Chem.MolFromSmiles(str(smiles))
            if mol:
                return Chem.MolToInchiKey(mol)
        except:
            pass
        return None
    
    def generate_compound_id(self, drugbank_id=None, smiles=None, name=None):
        """
        Generate canonical compound ID.
        Priority: DrugBank ID > InChIKey > Normalized name hash
        """
        # 1. DrugBank ID (best)
        if drugbank_id and not pd.isna(drugbank_id):
            return f"DB:{drugbank_id}"
        
        # 2. InChIKey from structure
        inchikey = self.get_inchikey(smiles)
        if inchikey:
            return f"IK:{inchikey[:14]}"  # First 14 chars = connectivity layer
        
        # 3. Normalized name hash
        norm_name = self.normalize_name(name)
        if norm_name:
            return f"NM:{norm_name}"
        
        return None
    
    def register_compound(self, drugbank_id=None, smiles=None, name=None, **metadata):
        """Register a compound in the registry."""
        compound_id = self.generate_compound_id(drugbank_id, smiles, name)
        if not compound_id:
            return None
        
        # Check for duplicates
        norm_name = self.normalize_name(name)
        canonical_smiles = self.get_canonical_smiles(smiles)
        
        # Register by name
        if norm_name and norm_name not in self.name_to_id:
            self.name_to_id[norm_name] = compound_id
        
        # Register by structure
        if canonical_smiles and canonical_smiles not in self.smiles_to_id:
            self.smiles_to_id[canonical_smiles] = compound_id
        
        # Store metadata
        if compound_id not in self.compound_registry:
            self.compound_registry[compound_id] = {
                'compound_id': compound_id,
                'drugbank_id': drugbank_id,
                'name': name,
                'normalized_name': norm_name,
                'smiles': smiles,
                'canonical_smiles': canonical_smiles,
                **metadata
            }
        
        return compound_id
    
    def find_compound(self, drugbank_id=None, smiles=None, name=None):
        """Find a compound in the registry by any identifier."""
        # Try DrugBank ID
        if drugbank_id and not pd.isna(drugbank_id):
            comp_id = f"DB:{drugbank_id}"
            if comp_id in self.compound_registry:
                return self.compound_registry[comp_id]
        
        # Try structure
        canonical_smiles = self.get_canonical_smiles(smiles)
        if canonical_smiles and canonical_smiles in self.smiles_to_id:
            comp_id = self.smiles_to_id[canonical_smiles]
            return self.compound_registry.get(comp_id)
        
        # Try name
        norm_name = self.normalize_name(name)
        if norm_name and norm_name in self.name_to_id:
            comp_id = self.name_to_id[norm_name]
            return self.compound_registry.get(comp_id)
        
        return None
    
    def detect_duplicates(self, df, name_col='drug_name', smiles_col='smiles'):
        """Detect duplicate compounds in a dataframe."""
        duplicates = []
        seen_ids = {}
        
        for idx, row in df.iterrows():
            comp_id = self.generate_compound_id(
                drugbank_id=row.get('drugbank_id'),
                smiles=row.get(smiles_col),
                name=row.get(name_col)
            )
            
            if comp_id in seen_ids:
                duplicates.append({
                    'idx1': seen_ids[comp_id],
                    'idx2': idx,
                    'compound_id': comp_id,
                    'name1': df.loc[seen_ids[comp_id], name_col],
                    'name2': row[name_col]
                })
            else:
                seen_ids[comp_id] = idx
        
        return duplicates

# =============================================================================
# APPLY TO OUR DATA
# =============================================================================

print("="*60)
print("Compound Identification System")
print("="*60)

# Initialize identifier
identifier = CompoundIdentifier()

# Load our data
print("\nLoading predictions...")
df = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'final_dili_predictions.csv')
print(f"Compounds: {len(df)}")

# Register all compounds
print("\nRegistering compounds...")
for idx, row in df.iterrows():
    identifier.register_compound(
        drugbank_id=row.get('drugbank_id'),
        smiles=row.get('smiles'),
        name=row.get('drug_name')
    )

print(f"Registered: {len(identifier.compound_registry)}")

# Check for duplicates
print("\nChecking for duplicates...")
duplicates = identifier.detect_duplicates(df)
if duplicates:
    print(f"Found {len(duplicates)} potential duplicates:")
    for dup in duplicates[:10]:
        print(f"  '{dup['name1']}' vs '{dup['name2']}' -> {dup['compound_id']}")
else:
    print("No duplicates found")

# Add canonical IDs to dataframe
df['compound_id'] = df.apply(
    lambda row: identifier.generate_compound_id(
        drugbank_id=row.get('drugbank_id'),
        smiles=row.get('smiles'),
        name=row.get('drug_name')
    ), axis=1
)

df['normalized_name'] = df['drug_name'].apply(identifier.normalize_name)

# =============================================================================
# TEST CROSS-DATASET MATCHING
# =============================================================================

print("\n" + "="*60)
print("Cross-Dataset Matching Test")
print("="*60)

# Test on known variants
test_names = [
    'ciclosporin',
    'cyclosporine',
    'paracetamol',
    'acetaminophen',
    'rifampin',
    'rifampicin',
    'metformin hydrochloride',
    'metformin'
]

print("\nName normalization test:")
for name in test_names:
    normalized = identifier.normalize_name(name)
    print(f"  {name:30s} -> {normalized}")

# =============================================================================
# SAVE ENHANCED DATA
# =============================================================================

output = ROOT / 'research' / 'ecnp-generalization' / 'results' / 'predictions_with_compound_ids.csv'
df.to_csv(output, index=False)
print(f"\nSaved: {output}")

print("\n" + "="*60)
print("COMPOUND IDENTIFICATION SUMMARY")
print("="*60)
print(f"""
Compounds processed: {len(df)}
Unique IDs generated: {len(identifier.compound_registry)}
Duplicates found: {len(duplicates)}

ID Types:
  DrugBank-based (DB:): {sum(1 for c in df['compound_id'] if c and c.startswith('DB:'))}
  Structure-based (IK:): {sum(1 for c in df['compound_id'] if c and c.startswith('IK:'))}
  Name-based (NM:): {sum(1 for c in df['compound_id'] if c and c.startswith('NM:'))}

Name synonyms in database: {len(DRUG_SYNONYMS)}
""")
