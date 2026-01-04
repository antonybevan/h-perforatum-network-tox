"""
Conformal Prediction Module
===========================
Implements Inductive Conformal Prediction (ICP) for classification.

References:
- Vovk, V., Gammerman, A., & Shafer, G. (2005). Algorithmic learning in a random world.
  Springer.
- Angelopoulos, A. N., & Bates, S. (2021). A Gentle Introduction to Conformal Prediction.
  arXiv:2107.07511. https://arxiv.org/abs/2107.07511
"""
import numpy as np
from typing import List, Tuple, Dict

from ..config import CONFORMAL_ALPHA


def compute_nonconformity(probs: np.ndarray, y_true_idx: np.ndarray) -> np.ndarray:
    """
    Compute non-conformity scores (1 - probability of true class).
    
    Args:
        probs: Probability matrix (N, K).
        y_true_idx: True class indices (N,).
    
    Returns:
        Non-conformity scores (N,).
    """
    return 1.0 - probs[np.arange(len(probs)), y_true_idx]


def calibrate_threshold(scores: np.ndarray, alpha: float = CONFORMAL_ALPHA) -> float:
    """
    Compute the conformal quantile threshold from calibration scores.
    
    Args:
        scores: Non-conformity scores from calibration set (N,).
        alpha: Miscoverage rate (default 0.1 for 90% coverage).
    
    Returns:
        Quantile threshold (q_hat).
    """
    n = len(scores)
    # Use the (1-alpha)(1+1/n) quantile for finite-sample correction
    q_val = np.quantile(scores, np.ceil((n + 1) * (1 - alpha)) / n, method='higher')
    return q_val


def get_prediction_set(probs: np.ndarray, q_hat: float) -> List[np.ndarray]:
    """
    Construct prediction sets based on calibrated threshold.
    
    Args:
        probs: Probability matrix (N, K).
        q_hat: Calibrated threshold.
    
    Returns:
        List of predicted class arrays.
    """
    threshold = 1.0 - q_hat
    sets = []
    for p_vec in probs:
        classes = np.where(p_vec >= threshold)[0]
        sets.append(classes)
    return sets


def evaluate_conformal(
    pred_sets: List[np.ndarray],
    y_true: np.ndarray,
    labels: np.ndarray = None
) -> Dict[str, float]:
    """
    Evaluate conformal prediction performance.
    
    Args:
        pred_sets: List of prediction sets.
        y_true: True labels (N,).
        labels: Optional array to compute label-conditional validity.
    
    Returns:
        Dict with: validity, mean_set_size, size_counts, ssc, lcv.
    """
    n = len(pred_sets)
    correct = np.array([y_true[i] in pred_sets[i] for i in range(n)])
    sizes = np.array([len(s) for s in pred_sets])
    
    results = {
        'validity': correct.mean(),
        'mean_set_size': sizes.mean(),
        'size_distribution': dict(zip(*np.unique(sizes, return_counts=True))),
    }
    
    # Size-Stratified Coverage (SSC)
    results['ssc'] = {}
    for sz in np.unique(sizes):
        mask = sizes == sz
        results['ssc'][sz] = correct[mask].mean()
    
    # Label-Conditional Validity (LCV)
    if labels is not None:
        results['lcv'] = {}
        for lbl in np.unique(labels):
            mask = y_true == lbl
            results['lcv'][int(lbl)] = correct[mask].mean()
    
    return results
