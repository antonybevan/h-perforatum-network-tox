"""
ECNP Interpretable Report Generator

For a PRESCRIBABLE safety assessment system, every prediction must be:
1. Quantified (Z-score)
2. Explained (which targets drive the signal)
3. Traceable (pathway connections)
4. Comparable (vs reference compounds)

This module generates human-readable reports for any compound.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import List, Dict, Optional
from ecnp_optimized import ECNPOptimized, ECNPConfig, ECNPStatus


@dataclass
class InfluenceReport:
    """Structured report for compound influence assessment."""
    compound_name: str
    targets: List[str]
    k: int
    z_score: float
    risk_tier: str  # HIGH, MODERATE, LOW, MINIMAL
    confidence: str  # HIGH, MEDIUM, LOW
    
    # Key metrics
    total_influence: float
    mean_influence: float
    ptni: float  # per-target normalized influence
    
    # Top contributors
    top_targets: List[Dict]  # [{gene, influence, pct_contribution, percentile}]
    
    # Pathway analysis
    hub_fraction: float  # fraction of targets in top 10%
    influence_concentration: float  # top 3 as fraction of total
    
    # DILI pathway connections  
    dili_connections: List[Dict]  # [{target, dili_gene, influence}]
    
    # Comparisons
    vs_hyperforin: Optional[float]  # ratio to reference
    vs_quercetin: Optional[float]
    
    # Status
    status: str
    warnings: List[str]


class ECNPReportGenerator:
    """Generate interpretable reports for compound safety assessment."""
    
    # Risk tier thresholds (based on manuscript Z-scores)
    RISK_TIERS = {
        'HIGH': 8.0,      # Z >= 8 (Hyperforin-like)
        'MODERATE': 4.0,  # Z >= 4 (Quercetin-like)
        'LOW': 2.0,       # Z >= 2
        'MINIMAL': 0.0    # Z < 2
    }
    
    # Reference Z-scores from Monte Carlo
    REF_HYPERFORIN = 10.27
    REF_QUERCETIN = 4.42
    
    def __init__(self, ecnp: ECNPOptimized = None):
        if ecnp is None:
            ecnp = ECNPOptimized()
        self.ecnp = ecnp
        
        # Load DILI genes
        dili_df = pd.read_csv(self.ecnp.data_dir / "dili_900_lcc.csv")
        self.dili_genes = set(dili_df['gene_name'].tolist())
    
    def generate_report(self, targets: List[str], compound_name: str = "Unknown") -> InfluenceReport:
        """Generate comprehensive influence report for a compound."""
        
        # Run ECNP
        result = self.ecnp.compute(targets)
        
        # Handle failures
        if result['status'] != ECNPStatus.SUCCESS and result['status'] != ECNPStatus.WARNING_LARGE_K:
            return InfluenceReport(
                compound_name=compound_name,
                targets=targets,
                k=len(targets),
                z_score=float('nan'),
                risk_tier="UNKNOWN",
                confidence="LOW",
                total_influence=0,
                mean_influence=0,
                ptni=0,
                top_targets=[],
                hub_fraction=0,
                influence_concentration=0,
                dili_connections=[],
                vs_hyperforin=None,
                vs_quercetin=None,
                status=result['status'].value,
                warnings=[result.get('message', 'Computation failed')]
            )
        
        # Extract values
        z = result['Z']
        k = result['k']
        I_T = result['I_T']
        
        # Determine risk tier
        risk_tier = self._get_risk_tier(z)
        
        # Analyze targets
        top_targets, hub_fraction, concentration = self._analyze_targets(targets)
        
        # Find DILI connections
        dili_connections = self._find_dili_connections(targets)
        
        # Determine confidence
        warnings = []
        if result['status'] == ECNPStatus.WARNING_LARGE_K:
            warnings.append(f"Large target set (k={k}) may have inflated variance")
        if hub_fraction < 0.2:
            warnings.append("Few hub targets - signal may be diffuse")
        if concentration < 0.3:
            warnings.append("Low influence concentration - distributed perturbation")
        
        confidence = self._get_confidence(result, hub_fraction, warnings)
        
        return InfluenceReport(
            compound_name=compound_name,
            targets=targets,
            k=k,
            z_score=z,
            risk_tier=risk_tier,
            confidence=confidence,
            total_influence=I_T,
            mean_influence=I_T / k if k > 0 else 0,
            ptni=I_T / k if k > 0 else 0,
            top_targets=top_targets,
            hub_fraction=hub_fraction,
            influence_concentration=concentration,
            dili_connections=dili_connections,
            vs_hyperforin=z / self.REF_HYPERFORIN,
            vs_quercetin=z / self.REF_QUERCETIN,
            status=result['status'].value,
            warnings=warnings
        )
    
    def _get_risk_tier(self, z: float) -> str:
        """Assign risk tier based on Z-score."""
        if z >= self.RISK_TIERS['HIGH']:
            return "HIGH"
        elif z >= self.RISK_TIERS['MODERATE']:
            return "MODERATE"
        elif z >= self.RISK_TIERS['LOW']:
            return "LOW"
        else:
            return "MINIMAL"
    
    def _analyze_targets(self, targets: List[str]):
        """Analyze target distribution."""
        target_data = []
        for t in targets:
            if t in self.ecnp.node_to_idx:
                idx = self.ecnp.node_to_idx[t]
                m = self.ecnp.m_array[idx]
                pct = self.ecnp.percentiles_array[idx]
                target_data.append({'gene': t, 'influence': m, 'percentile': pct})
        
        if not target_data:
            return [], 0.0, 0.0
        
        # Sort by influence
        target_data.sort(key=lambda x: -x['influence'])
        
        total = sum(t['influence'] for t in target_data)
        
        # Add percentage contribution
        for t in target_data:
            t['pct_contribution'] = t['influence'] / total if total > 0 else 0
        
        # Hub fraction (top 10% influence percentile)
        hub_fraction = sum(1 for t in target_data if t['percentile'] >= 0.90) / len(target_data)
        
        # Concentration (top 3)
        top3_sum = sum(t['influence'] for t in target_data[:3])
        concentration = top3_sum / total if total > 0 else 0
        
        return target_data[:10], hub_fraction, concentration
    
    def _find_dili_connections(self, targets: List[str], top_n: int = 10):
        """Find strongest target -> DILI gene connections."""
        connections = []
        
        for t in targets:
            if t not in self.ecnp.node_to_idx:
                continue
            
            target_idx = self.ecnp.node_to_idx[t]
            
            for dili_gene in self.dili_genes:
                if dili_gene not in self.ecnp.node_to_idx:
                    continue
                
                dili_idx = self.ecnp.node_to_idx[dili_gene]
                influence = self.ecnp.M[dili_idx, target_idx]
                
                if influence > 0.01:  # Threshold for meaningful connection
                    connections.append({
                        'target': t,
                        'dili_gene': dili_gene,
                        'influence': influence
                    })
        
        # Sort and return top N
        connections.sort(key=lambda x: -x['influence'])
        return connections[:top_n]
    
    def _get_confidence(self, result, hub_fraction, warnings) -> str:
        """Determine confidence level."""
        if len(warnings) >= 2:
            return "LOW"
        elif len(warnings) == 1 or hub_fraction < 0.3:
            return "MEDIUM"
        else:
            return "HIGH"
    
    def print_report(self, report: InfluenceReport):
        """Print formatted report."""
        print("\n" + "=" * 70)
        print(f"ECNP INFLUENCE ASSESSMENT: {report.compound_name}")
        print("=" * 70)
        
        print(f"\n[SUMMARY]")
        print(f"  Risk Tier:    {report.risk_tier}")
        print(f"  Z-Score:      {report.z_score:.2f}")
        print(f"  Confidence:   {report.confidence}")
        print(f"  Targets:      {report.k}")
        
        print(f"\n[COMPARISON TO REFERENCE COMPOUNDS]")
        if report.vs_hyperforin:
            print(f"  vs Hyperforin (Z=10.27): {report.vs_hyperforin:.1%}")
        if report.vs_quercetin:
            print(f"  vs Quercetin (Z=4.42):   {report.vs_quercetin:.1%}")
        
        print(f"\n[TARGET ANALYSIS]")
        print(f"  Hub targets (top 10%):    {report.hub_fraction:.0%}")
        print(f"  Influence concentration:  Top 3 = {report.influence_concentration:.0%}")
        
        print(f"\n[TOP CONTRIBUTORS]")
        print(f"  {'Gene':<12} | {'Influence':>10} | {'% Total':>8} | {'Percentile':>10}")
        print(f"  {'-'*48}")
        for t in report.top_targets[:5]:
            print(f"  {t['gene']:<12} | {t['influence']:>10.4f} | {t['pct_contribution']:>7.1%} | {t['percentile']:>9.1%}")
        
        if report.dili_connections:
            print(f"\n[KEY DILI PATHWAY CONNECTIONS]")
            print(f"  {'Target':<12} -> {'DILI Gene':<12} | {'Influence':>10}")
            print(f"  {'-'*42}")
            for c in report.dili_connections[:5]:
                print(f"  {c['target']:<12} -> {c['dili_gene']:<12} | {c['influence']:>10.4f}")
        
        if report.warnings:
            print(f"\n[WARNINGS]")
            for w in report.warnings:
                print(f"  ! {w}")
        
        print(f"\n[INTERPRETATION]")
        if report.risk_tier == "HIGH":
            print("  This compound shows strong DILI-directed network influence,")
            print("  comparable to or exceeding Hyperforin. High-leverage targets")
            print("  suggest efficient perturbation of hepatotoxicity-relevant pathways.")
        elif report.risk_tier == "MODERATE":
            print("  This compound shows moderate DILI-directed network influence,")
            print("  comparable to Quercetin. The signal warrants attention but")
            print("  may reflect distributed rather than focused perturbation.")
        elif report.risk_tier == "LOW":
            print("  This compound shows weak DILI-directed network influence.")
            print("  Network topology does not suggest strong hepatotoxicity risk,")
            print("  though other mechanisms may still be relevant.")
        else:
            print("  This compound shows minimal DILI-directed network influence.")
            print("  Network-based prioritization does not flag this compound.")
        
        print("\n" + "=" * 70)


def demo():
    """Demonstrate report generation."""
    gen = ECNPReportGenerator()
    
    targets_df = pd.read_csv(gen.ecnp.data_dir / "targets_lcc.csv")
    hyperforin = targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'].tolist()
    quercetin = targets_df[targets_df['compound'] == 'Quercetin']['gene_symbol'].tolist()
    
    # Generate reports
    report_hyp = gen.generate_report(hyperforin, "Hyperforin")
    gen.print_report(report_hyp)
    
    report_que = gen.generate_report(quercetin, "Quercetin")
    gen.print_report(report_que)
    
    # Test with random targets
    print("\n\n" + "#" * 70)
    print("TESTING WITH RANDOM TARGETS")
    print("#" * 70)
    
    np.random.seed(42)
    random_targets = list(np.random.choice(gen.ecnp.node_list, 15, replace=False))
    report_random = gen.generate_report(random_targets, "Random-15")
    gen.print_report(report_random)


if __name__ == "__main__":
    demo()
