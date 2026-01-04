"""
ECNP Closed-Form Algorithm
==========================

Two-Layer Statistical Architecture for Network Toxicology.

Modules:
    core/          - Production algorithms (Layer 1 + Layer 2)
    precompute/    - One-time matrix computation
    validation/    - Test suite and stress tests
    analysis/      - Power analysis and calibration

Quick Start:
    from src.core.ecnp_optimized import ECNPOptimized
    from src.core.ecnp_permutation_test import StratifiedPermutationTest

Author: Research Session 2026-01-02
"""

__version__ = "1.0.0"
__all__ = ["core", "precompute", "validation", "analysis"]
