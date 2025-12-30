#!/usr/bin/env python3
"""
Master Pipeline: Complete Reproducible Analysis

This script orchestrates the entire H. perforatum network toxicology analysis
from raw data to final results, ensuring 100% reproducibility.

Usage:
    python scripts/run_pipeline.py              # Run all steps
    python scripts/run_pipeline.py --step 3     # Run from step 3 onwards
    python scripts/run_pipeline.py --only 5     # Run only step 5
    python scripts/run_pipeline.py --validate   # Run validation only
    
Pipeline Steps:
    1. Extract STRING networks (700 & 900 thresholds)
    2. Regenerate targets from raw data
    3. Regenerate DILI genes from DisGeNET
    4. Create LCC-filtered data
    5. Validate data integrity (27 checks)
    6. Run Standard RWR permutation analysis
    7. Run Expression-weighted RWR permutation analysis
    8. Run Shortest Path permutation analysis
    9. Run Bootstrap sensitivity analysis
    10. Run Chemical Similarity control
    11. Generate DATA_FLOW.md documentation

Author: Network Toxicology Pipeline v2.0
"""

import subprocess
import sys
import time
import argparse
from pathlib import Path
from datetime import datetime

# Configuration
SCRIPTS_DIR = Path('scripts')
RESULTS_DIR = Path('results')
DATA_DIR = Path('data')

# Pipeline steps in execution order
PIPELINE_STEPS = [
    {
        'step': 1,
        'name': 'Extract STRING Networks',
        'script': None,  # Custom handler
        'description': 'Extract networks from STRING v12.0 at 700 and 900 thresholds',
        'outputs': ['data/processed/network_700.parquet', 'data/processed/network_900.parquet'],
        'skip_if_exists': True,
    },
    {
        'step': 2,
        'name': 'Regenerate Targets',
        'script': 'regenerate_targets.py',
        'description': 'Create targets.csv from raw data with documented filtering',
        'outputs': ['data/processed/targets.csv'],
    },
    {
        'step': 3,
        'name': 'Regenerate DILI Genes',
        'script': 'regenerate_dili.py',
        'description': 'Extract DILI genes from DisGeNET source',
        'outputs': ['data/raw/dili_genes_raw.csv', 'data/processed/dili_700_lcc.csv', 'data/processed/dili_900_lcc.csv'],
    },
    {
        'step': 4,
        'name': 'Create LCC-Filtered Data',
        'script': 'create_lcc_filtered_data.py',
        'description': 'Filter targets and networks to liver-expressed LCC',
        'outputs': ['data/processed/targets_lcc.csv', 'data/processed/network_700_liver_lcc.parquet'],
    },
    {
        'step': 5,
        'name': 'Validate Data Integrity',
        'script': 'validate_data_integrity.py',
        'description': 'Run 27 automated integrity checks',
        'outputs': [],
        'required': True,  # Must pass for pipeline to continue
    },
    {
        'step': 6,
        'name': 'Standard RWR Analysis',
        'script': 'run_standard_rwr_lcc_permutations.py',
        'description': 'Degree-matched permutation testing (1000 permutations)',
        'outputs': ['results/tables/standard_rwr_lcc_permutation_results.csv'],
    },
    {
        'step': 7,
        'name': 'Expression-Weighted RWR Analysis',
        'script': 'run_expression_weighted_rwr_permutations.py',
        'description': 'Liver expression-weighted RWR with permutation testing',
        'outputs': ['results/tables/expression_weighted_rwr_permutation_results.csv'],
    },
    {
        'step': 8,
        'name': 'Shortest Path Analysis',
        'script': 'run_shortest_path_permutations.py',
        'description': 'Network proximity to DILI genes with permutation testing',
        'outputs': ['results/tables/shortest_path_permutation_results.csv'],
    },
    {
        'step': 9,
        'name': 'Bootstrap Sensitivity Analysis',
        'script': 'run_bootstrap_sensitivity.py',
        'description': 'Test robustness to target count asymmetry',
        'outputs': ['results/bootstrap_sensitivity.csv', 'results/tables/bootstrap_summary.csv'],
    },
    {
        'step': 10,
        'name': 'Chemical Similarity Control',
        'script': 'run_chemical_similarity_control.py',
        'description': 'Tanimoto similarity vs FDA DILIrank hepatotoxins',
        'outputs': ['results/tables/chemical_similarity_summary.csv'],
    },
    {
        'step': 11,
        'name': 'Generate Documentation',
        'script': 'generate_dataflow.py',
        'description': 'Auto-generate DATA_FLOW.md with complete traceability',
        'outputs': ['docs/DATA_FLOW.md'],
    },
]


def print_header(text, char='='):
    """Print formatted header."""
    width = 70
    print()
    print(char * width)
    print(f" {text}")
    print(char * width)


def print_step(step_info, status='RUNNING'):
    """Print step information."""
    icons = {'RUNNING': '🔄', 'DONE': '✅', 'SKIP': '⏭️', 'FAIL': '❌'}
    icon = icons.get(status, '•')
    print(f"\n{icon} Step {step_info['step']}: {step_info['name']}")
    print(f"   {step_info['description']}")


def run_script(script_name, capture_output=False):
    """Run a Python script and return success status."""
    script_path = SCRIPTS_DIR / script_name
    
    if not script_path.exists():
        print(f"   ERROR: Script not found: {script_path}")
        return False
    
    cmd = [sys.executable, str(script_path)]
    
    try:
        if capture_output:
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"   STDERR: {result.stderr[:500]}")
                return False
        else:
            result = subprocess.run(cmd)
            if result.returncode != 0:
                return False
        return True
    except Exception as e:
        print(f"   Exception: {e}")
        return False


def run_string_extraction():
    """Run STRING network extraction for both thresholds."""
    for threshold in [700, 900]:
        output_file = DATA_DIR / 'processed' / f'network_{threshold}.parquet'
        
        if output_file.exists():
            print(f"   Network {threshold} exists, skipping...")
            continue
        
        print(f"   Extracting STRING network (≥{threshold})...")
        cmd = [
            sys.executable, str(SCRIPTS_DIR / 'extract_string_network.py'),
            '--threshold', str(threshold),
            '--output', str(output_file)
        ]
        result = subprocess.run(cmd)
        if result.returncode != 0:
            return False
    return True


def check_outputs_exist(outputs):
    """Check if all output files exist."""
    for output in outputs:
        if not Path(output).exists():
            return False
    return True


def run_pipeline(start_step=1, only_step=None, validate_only=False):
    """Run the complete pipeline."""
    
    print_header("H. PERFORATUM NETWORK TOXICOLOGY PIPELINE v2.0")
    print(f"\nStarted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Python: {sys.version.split()[0]}")
    print(f"Working directory: {Path.cwd()}")
    
    if validate_only:
        print("\n⚠️  VALIDATION MODE: Running data integrity checks only")
        start_step = 5
        only_step = 5
    
    if only_step:
        print(f"\n⚠️  SINGLE STEP MODE: Running step {only_step} only")
    elif start_step > 1:
        print(f"\n⚠️  RESUME MODE: Starting from step {start_step}")
    
    # Track results
    results = []
    start_time = time.time()
    
    for step_info in PIPELINE_STEPS:
        step_num = step_info['step']
        
        # Skip if before start_step
        if step_num < start_step:
            continue
        
        # Skip if only_step specified and not matching
        if only_step and step_num != only_step:
            continue
        
        print_step(step_info, 'RUNNING')
        step_start = time.time()
        
        # Check if we can skip (outputs exist)
        if step_info.get('skip_if_exists') and check_outputs_exist(step_info['outputs']):
            print("   Outputs exist, skipping...")
            print_step(step_info, 'SKIP')
            results.append({'step': step_num, 'status': 'SKIPPED', 'time': 0})
            continue
        
        # Run the step
        if step_num == 1:
            success = run_string_extraction()
        else:
            success = run_script(step_info['script'])
        
        elapsed = time.time() - step_start
        
        if success:
            print_step(step_info, 'DONE')
            print(f"   Completed in {elapsed:.1f}s")
            results.append({'step': step_num, 'status': 'SUCCESS', 'time': elapsed})
        else:
            print_step(step_info, 'FAIL')
            results.append({'step': step_num, 'status': 'FAILED', 'time': elapsed})
            
            # Stop if required step fails
            if step_info.get('required'):
                print("\n❌ PIPELINE HALTED: Required step failed")
                break
    
    # Summary
    total_time = time.time() - start_time
    print_header("PIPELINE SUMMARY")
    
    print("\nStep Results:")
    for r in results:
        status_icon = {'SUCCESS': '✅', 'FAILED': '❌', 'SKIPPED': '⏭️'}[r['status']]
        step_name = next(s['name'] for s in PIPELINE_STEPS if s['step'] == r['step'])
        print(f"  {status_icon} Step {r['step']}: {step_name} ({r['time']:.1f}s)")
    
    failed = [r for r in results if r['status'] == 'FAILED']
    
    print(f"\nTotal time: {total_time:.1f}s ({total_time/60:.1f} min)")
    print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    if failed:
        print(f"\n❌ PIPELINE FAILED: {len(failed)} step(s) failed")
        return 1
    else:
        print("\n✅ PIPELINE COMPLETED SUCCESSFULLY")
        return 0


def main():
    parser = argparse.ArgumentParser(
        description='H. perforatum Network Toxicology Master Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/run_pipeline.py              # Run complete pipeline
  python scripts/run_pipeline.py --step 6     # Resume from step 6
  python scripts/run_pipeline.py --only 5     # Run validation only
  python scripts/run_pipeline.py --validate   # Quick validation check
        """
    )
    
    parser.add_argument('--step', type=int, default=1,
                       help='Start from this step (1-11)')
    parser.add_argument('--only', type=int,
                       help='Run only this step')
    parser.add_argument('--validate', action='store_true',
                       help='Run data validation only (step 5)')
    parser.add_argument('--list', action='store_true',
                       help='List all pipeline steps')
    
    args = parser.parse_args()
    
    if args.list:
        print_header("PIPELINE STEPS")
        for step in PIPELINE_STEPS:
            print(f"\n  Step {step['step']}: {step['name']}")
            print(f"         {step['description']}")
            if step.get('script'):
                print(f"         Script: {step['script']}")
        print()
        return 0
    
    return run_pipeline(
        start_step=args.step,
        only_step=args.only,
        validate_only=args.validate
    )


if __name__ == '__main__':
    sys.exit(main())
