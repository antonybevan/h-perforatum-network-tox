#!/usr/bin/env python3
"""
H. perforatum Network Toxicology Pipeline

Orchestrates the analysis from raw data to final results.

Usage:
    python scripts/run_pipeline.py              # Run all steps
    python scripts/run_pipeline.py --step 3     # Run from step 3 onwards
    python scripts/run_pipeline.py --only 5     # Run only step 5
    python scripts/run_pipeline.py --validate   # Run validation only
"""

import subprocess
import sys
import time
import argparse
from pathlib import Path
from datetime import datetime

try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False

# Configuration
SCRIPTS_DIR = Path('scripts')
RESULTS_DIR = Path('results')
DATA_DIR = Path('data')

# Pipeline steps in execution order
PIPELINE_STEPS = [
    {
        'step': 1,
        'name': 'Extract STRING Networks',
        'script': None,
        'description': 'Extract networks from STRING v12.0 at 700 and 900 thresholds',
        'outputs': [DATA_DIR / 'processed' / 'network_700.parquet', DATA_DIR / 'processed' / 'network_900.parquet'],
        'skip_if_exists': True,
    },
    {
        'step': 2,
        'name': 'Regenerate Targets',
        'script': 'regenerate_targets.py',
        'description': 'Create targets.csv from raw data with documented filtering',
        'outputs': [DATA_DIR / 'processed' / 'targets.csv'],
    },
    {
        'step': 3,
        'name': 'Regenerate DILI Genes',
        'script': 'regenerate_dili.py',
        'description': 'Extract DILI genes from DisGeNET source',
        'outputs': [DATA_DIR / 'raw' / 'dili_genes_raw.csv', DATA_DIR / 'processed' / 'dili_700_lcc.csv', DATA_DIR / 'processed' / 'dili_900_lcc.csv'],
    },
    {
        'step': 4,
        'name': 'Create LCC-Filtered Data',
        'script': 'create_lcc_filtered_data.py',
        'description': 'Filter targets and networks to liver-expressed LCC',
        'outputs': [DATA_DIR / 'processed' / 'targets_lcc.csv', DATA_DIR / 'processed' / 'network_700_liver_lcc.parquet'],
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
        'outputs': [RESULTS_DIR / 'tables' / 'standard_rwr_lcc_permutation_results.csv'],
    },
    {
        'step': 7,
        'name': 'Expression-Weighted RWR Analysis',
        'script': 'run_expression_weighted_rwr_permutations.py',
        'description': 'Liver expression-weighted RWR with permutation testing',
        'outputs': [RESULTS_DIR / 'tables' / 'expression_weighted_rwr_permutation_results.csv'],
    },
    {
        'step': 8,
        'name': 'Shortest Path Analysis',
        'script': 'run_shortest_path_permutations.py',
        'description': 'Network proximity to DILI genes with permutation testing',
        'outputs': [RESULTS_DIR / 'tables' / 'shortest_path_permutation_results.csv'],
    },
    {
        'step': 9,
        'name': 'Bootstrap Sensitivity Analysis',
        'script': 'run_bootstrap_sensitivity.py',
        'description': 'Test robustness to target count asymmetry',
        'outputs': [RESULTS_DIR / 'bootstrap_sensitivity.csv', RESULTS_DIR / 'tables' / 'bootstrap_summary.csv'],
    },
    {
        'step': 10,
        'name': 'EWI Bootstrap Sensitivity Analysis',
        'script': 'run_ewi_bootstrap_sensitivity.py',
        'description': 'Test robustness of expression-weighted results',
        'outputs': [RESULTS_DIR / 'tables' / 'ewi_bootstrap_summary.csv'],
    },
    {
        'step': 11,
        'name': 'Chemical Similarity Control',
        'script': 'run_chemical_similarity_control.py',
        'description': 'Tanimoto similarity vs FDA DILIrank hepatotoxins',
        'outputs': [RESULTS_DIR / 'tables' / 'chemical_similarity_summary.csv'],
    },
    {
        'step': 12,
        'name': 'Generate Documentation',
        'script': 'generate_dataflow.py',
        'description': 'Auto-generate DATA_FLOW.md with complete traceability',
        'outputs': [Path('docs') / 'DATA_FLOW.md'],
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
    """Print step info."""
    icons = {'RUNNING': '...', 'DONE': 'OK', 'SKIP': 'SKIP', 'FAIL': 'FAIL'}
    icon = icons.get(status, '-')
    print(f"\n[{icon}] Step {step_info['step']}: {step_info['name']}")
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
    
    print_header("NETWORK TOXICOLOGY PIPELINE")
    print(f"\nStarted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Working directory: {Path.cwd()}")
    
    if validate_only:
        print("\nVAL MODE: Running data integrity checks only")
        start_step = 5
        only_step = 5
    
    if only_step:
        print(f"\nSINGLE STEP MODE: Running step {only_step} only")
    elif start_step > 1:
        print(f"\nRESUME MODE: Starting from step {start_step}")
    
    # Track results
    results = []
    start_time = time.time()
    
    # Filter steps to run
    steps_to_run = [s for s in PIPELINE_STEPS 
                    if s['step'] >= start_step and (only_step is None or s['step'] == only_step)]
    
    # Create progress bar if tqdm available
    if TQDM_AVAILABLE and len(steps_to_run) > 1:
        pbar = tqdm(
            steps_to_run,
            desc="Pipeline Progress",
            unit="step",
            ncols=80,
            bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]'
        )
    else:
        pbar = steps_to_run
    
    for step_info in pbar:
        step_num = step_info['step']
        
        # Update progress bar description
        if TQDM_AVAILABLE and hasattr(pbar, 'set_description'):
            pbar.set_description(f"Step {step_num}: {step_info['name'][:25]}")
        
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
            
            if step_info.get('required'):
                print("\nHALTED: Required step failed")
                break
    
    # Summary
    total_time = time.time() - start_time
    print_header("PIPELINE SUMMARY")
    
    print("\nResults:")
    for r in results:
        status_icon = {'SUCCESS': 'OK  ', 'FAILED': 'FAIL', 'SKIPPED': 'SKIP'}[r['status']]
        step_name = next(s['name'] for s in PIPELINE_STEPS if s['step'] == r['step'])
        print(f"  [{status_icon}] Step {r['step']}: {step_name} ({r['time']:.1f}s)")
    
    failed = [r for r in results if r['status'] == 'FAILED']
    
    print(f"\nTotal time: {total_time:.1f}s ({total_time/60:.1f} min)")
    print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    if failed:
        print(f"\nFAILED: {len(failed)} step(s) failed")
        return 1
    else:
        print("\nCOMPLETED SUCCESSFULLY")
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
