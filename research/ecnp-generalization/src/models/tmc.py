"""
Trusted Multi-View Classifier (TMC) with Evidential Deep Learning
===================================================================
Fuses chemical fingerprints and network features using EDL for uncertainty.

References:
- Sensoy, M., Kaplan, L., & Kandemir, M. (2018). Evidential deep learning to quantify
  classification uncertainty. NeurIPS. https://arxiv.org/abs/1806.01768
- Han, Z., Zhang, C., Fu, H., & Zhou, J. (2021). Trusted Multi-View Classification.
  ICLR. https://arxiv.org/abs/2102.02051
"""
import numpy as np
from scipy.optimize import minimize
from sklearn.preprocessing import StandardScaler
from typing import Tuple, Dict, Any

from ..config import RANDOM_SEED

np.random.seed(RANDOM_SEED)


def softplus(x: np.ndarray) -> np.ndarray:
    """Numerically stable softplus activation."""
    return np.log1p(np.exp(np.clip(x, -20, 20)))


def edl_loss_numpy(weights: np.ndarray, X: np.ndarray, y_onehot: np.ndarray, num_classes: int) -> float:
    """
    Evidential Deep Learning loss function using NumPy.
    
    Combines cross-entropy loss with a variance penalty to encourage
    evidence accumulation for the correct class.
    
    Args:
        weights: Flattened weight matrix (D*K,).
        X: Feature matrix (N, D).
        y_onehot: One-hot encoded labels (N, K).
        num_classes: Number of output classes.
    
    Returns:
        Scalar loss value.
    """
    D = X.shape[1]
    K = num_classes
    W = weights.reshape(D, K)
    logits = X @ W
    evidence = softplus(logits)
    alpha = evidence + 1  # Dirichlet concentration
    S = np.sum(alpha, axis=1, keepdims=True)
    
    # Mean Squared Error term
    err = np.sum((y_onehot - (alpha / S))**2, axis=1, keepdims=True)
    # Variance term (encourages peakedness)
    var = np.sum((alpha * (S - alpha)) / (S * S * (S + 1)), axis=1, keepdims=True)
    # L2 regularization
    reg = 0.01 * np.sum(W**2)
    
    return np.mean(err + var) + reg


def train_view_model(X: np.ndarray, y_onehot: np.ndarray) -> np.ndarray:
    """
    Train a single-view EDL model using L-BFGS-B optimization.
    
    Args:
        X: Feature matrix (N, D).
        y_onehot: One-hot encoded labels (N, K).
    
    Returns:
        Trained weight matrix (D, K).
    """
    D = X.shape[1]
    K = y_onehot.shape[1]
    initial_weights = np.random.randn(D * K) * 0.01
    
    res = minimize(
        fun=edl_loss_numpy,
        x0=initial_weights,
        args=(X, y_onehot, K),
        method='L-BFGS-B',
        options={'maxiter': 200, 'disp': False}
    )
    return res.x.reshape(D, K)


def predict_evidence(X: np.ndarray, W: np.ndarray) -> np.ndarray:
    """
    Predict evidence for each class.
    
    Args:
        X: Feature matrix (N, D).
        W: Weight matrix (D, K).
    
    Returns:
        Evidence matrix (N, K).
    """
    return softplus(X @ W)


class TrustedMultiViewClassifier:
    """
    Two-view classifier (Chemical + Network) with symbolic gating.
    
    The network view is gated by a "quality" signal (e.g., target count).
    """
    
    def __init__(self, quality_threshold: float = 0.693):
        self.W_chem = None
        self.W_net = None
        self.scaler_chem = StandardScaler()
        self.scaler_net = StandardScaler()
        self.quality_threshold = quality_threshold
    
    def fit(self, X_chem: np.ndarray, X_net: np.ndarray, y: np.ndarray):
        """Train both view models."""
        y_onehot = np.zeros((len(y), 2))
        y_onehot[np.arange(len(y)), y] = 1
        
        X_chem_scaled = self.scaler_chem.fit_transform(X_chem)
        X_net_scaled = self.scaler_net.fit_transform(X_net)
        
        self.W_chem = train_view_model(X_chem_scaled, y_onehot)
        self.W_net = train_view_model(X_net_scaled, y_onehot)
    
    def predict_proba(self, X_chem: np.ndarray, X_net: np.ndarray, quality: np.ndarray) -> np.ndarray:
        """
        Predict probabilities with symbolic gating.
        
        Args:
            X_chem: Chemical features (N, D_chem).
            X_net: Network features (N, D_net).
            quality: Quality signal (N,). Gate opens if quality > threshold.
        
        Returns:
            Probability matrix (N, 2).
        """
        X_chem_scaled = self.scaler_chem.transform(X_chem)
        X_net_scaled = self.scaler_net.transform(X_net)
        
        e_chem = predict_evidence(X_chem_scaled, self.W_chem)
        e_net = predict_evidence(X_net_scaled, self.W_net)
        
        # Symbolic Gate
        gate = (quality > self.quality_threshold).astype(float).reshape(-1, 1)
        e_net_gated = e_net * gate
        
        # Fusion
        e_fused = e_chem + e_net_gated
        alpha = e_fused + 1
        return alpha / np.sum(alpha, axis=1, keepdims=True)
