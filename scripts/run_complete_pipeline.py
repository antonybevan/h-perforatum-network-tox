"""
Master Pipeline: Complete Reproducibility from Source

This script orchestrates the entire analysis pipeline from raw data to final results.
Running this script will regenerate ALL processed files and analysis results from source.

Usage:
    python scripts/run_complete_pipeline.py [--skip-data] [--skip-validation]
    
    --skip-data        Skip data regeneration (use existing processed files)
    --skip-validation  Skip validation steps (only regenerate data)
"""

import subprocess
import sys
import argparse
import time

def run_step(step_name, command, critical=True):
    """Run a pipeline step."""
    print("\n" + "="*80)
    print(f"STEP: {step_name}")
    print("="*80)
    print(f"Command: {' '.join(command)}")
    print()
    
    start = time.time()
    result = subprocess.run(command, capture_output=False)
    elapsed = time.time() - start
    
    if result.returncode != 0:
        print(f"\n[ERR] ERROR: {step_name} failed!")
        if critical:
            print("Pipeline stopped due to critical error.")
            sys.exit(1)
        else:
            print("Continuing despite error (non-critical step)...")
    else:
        print(f"\n[OK] {step_name} complete! ({elapsed:.1f}s)")
    
    return result.returncode == 0

def main():
    parser = argparse.ArgumentParser(description='Run complete analysis pipeline')
    parser.add_argument('--skip-data', action='store_true',
                       help='Skip data regeneration steps')
    parser.add_argument('--skip-validation', action='store_true',
                       help='Skip validation steps')
    parser.add_argument('--quick', action='store_true',
                       help='Quick mode (skip non-essential steps)')
    
    args = parser.parse_args()
    
    print("="*80)
    print("COMPLETE ANALYSIS PIPELINE")
    print("H. perforatum Network Pharmacology")
    print("="*80)
    print("\nConfiguration:")
    print(f"  Skip data regeneration: {args.skip_data}")
    print(f"  Skip validation: {args.skip_validation}")
    print(f"  Quick mode: {args.quick}")
    
    # Track results
    results = {}
    
    # =========================================================================
    # PHASE 1: DATA REGENERATION
    # =========================================================================
    if not args.skip_data:
        print("\n" + "="*80)
        print("PHASE 1: DATA REGENERATION FROM SOURCE")
        print("="*80)
        
        # Step 1: Regenerate liver proteome from GTEx
        results['liver_proteome'] = run_step(
            "Regenerate Liver Proteome",
            ['python', 'scripts/regenerate_liver_proteome.py'],
            critical=True
        )
        
        # Step 2: Regenerate networks from STRING (optional - slow!)
        if not args.quick:
            results['network_900'] = run_step(
                "Extract STRING Network (≥900)",
                ['python', 'scripts/extract_string_network.py',
                 '--threshold', '900',
                 '--output', 'data/processed/network_900.parquet'],
                critical=False  # Non-critical (can use existing)
            )
            
            results['network_700'] = run_step(
                "Extract STRING Network (≥700)",
                ['python', 'scripts/extract_string_network.py',
                 '--threshold', '700',
                 '--output', 'data/processed/network_700.parquet'],
                critical=False
            )
        
        # Step 3: Regenerate DILI modules
        results['dili_modules'] = run_step(
            "Regenerate DILI Modules",
            ['python', 'scripts/regenerate_dili.py'],
            critical=True
        )
        
        # Step 4: Curate targets (optional - slow!)
        if not args.quick:
            results['targets'] = run_step(
                "Curate Targets from Raw",
                ['python', 'scripts/curate_targets.py'],
                critical=False
            )
    
    # =========================================================================
    # PHASE 2: VALIDATION & ANALYSIS
    # =========================================================================
    if not args.skip_validation:
        print("\n" + "="*80)
        print("PHASE 2: STATISTICAL VALIDATION")
        print("="*80)
        
        # Step 5: Primary validation (>=900)
        results['primary_validation'] = run_step(
            "Primary Validation (>=900, 1000 permutations)",
            ['python', 'scripts/run_full_validation.py'],
            critical=True
        )
        
        # Step 6: Robustness validation (>=700)
        results['robustness_validation'] = run_step(
            "Robustness Validation (>=700, 1000 permutations)",
            ['python', 'scripts/run_validation_700.py'],
            critical=True
        )
        
        # Step 7: Bootstrap sensitivity
        results['bootstrap'] = run_step(
            "Bootstrap Sensitivity Analysis (100 iterations)",
            ['python', 'scripts/run_bootstrap_sensitivity.py'],
            critical=True
        )
    
    # =========================================================================
    # PHASE 3: VERIFICATION
    # =========================================================================
    print("\n" + "="*80)
    print("PHASE 3: FINAL VERIFICATION")
    print("="*80)
    
    # Step 8: Run validation checks
    results['final_check'] = run_step(
        "Final Validation Checks",
        ['python', 'scripts/final_validation_check.py'],
        critical=True
    )
    
    # =========================================================================
    # SUMMARY
    # =========================================================================
    print("\n" + "="*80)
    print("PIPELINE COMPLETE")
    print("="*80)
    
    print("\nResults Summary:")
    for step, success in results.items():
        status = "[OK]" if success else "[FAIL]"
        print(f"  {status} {step}")
    
    failed = [k for k, v in results.items() if not v]
    if failed:
        print(f"\n[!] {len(failed)} steps failed: {', '.join(failed)}")
        print("Review errors above for details.")
        sys.exit(1)
    else:
        print("\n[OK] All steps completed successfully!")
        print("\nGenerated files:")
        print("  - data/processed/liver_proteome.csv")
        print("  - data/processed/dili_900_lcc.csv, dili_700_lcc.csv")
        print("  - results/final_statistics.csv")
        print("  - results/final_statistics_700.csv")
        print("  - results/bootstrap_sensitivity.csv")
        print("\nAnalysis is complete and reproducible!")

if __name__ == '__main__':
    main()
