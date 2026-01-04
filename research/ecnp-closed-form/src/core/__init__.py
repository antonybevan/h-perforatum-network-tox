"""
Core ECNP algorithms for production use.

Layer 1 (ecnp_optimized): Fast ranking score (~1ms)
Layer 2 (ecnp_permutation_test): Valid p-values (~1-7ms for 10K perms)
"""

from .ecnp_optimized import ECNPOptimized, ECNPConfig
from .ecnp_permutation_test import StratifiedPermutationTest

__all__ = ["ECNPOptimized", "ECNPConfig", "StratifiedPermutationTest"]
