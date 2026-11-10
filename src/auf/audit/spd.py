"""SPD-manifold primitives for the Stage-1 audit.

Everything here is batched: inputs of shape (..., d, d) are symmetric positive
definite matrices, and reductions act on the trailing two axes.
"""
from __future__ import annotations

import numpy as np

__all__ = [
    "isotropize",
    "fractional_anisotropy",
    "airm_distance",
    "logfrob_distance",
]


def _eigh(sigma: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    w, v = np.linalg.eigh(sigma)
    if np.any(w <= 0):
        raise ValueError("matrix is not positive definite")
    return w, v


def isotropize(sigma: np.ndarray) -> np.ndarray:
    """Trace-matched isotropic surrogate ``(tr Σ / d) I``.

    This is the scalar baseline the published covariance models do not report.
    It keeps the magnitude of the uncertainty and discards its direction.
    """
    sigma = np.asarray(sigma, dtype=float)
    d = sigma.shape[-1]
    scale = np.trace(sigma, axis1=-2, axis2=-1) / d
    eye = np.eye(d)
    return scale[..., None, None] * eye


def fractional_anisotropy(sigma: np.ndarray) -> np.ndarray:
    """Fractional anisotropy in [0, 1): 0 is a sphere, 1 a degenerate line.

    Borrowed from diffusion-tensor imaging, where the same 3x3 SPD object has
    been studied for decades.
    """
    sigma = np.asarray(sigma, dtype=float)
    w, _ = _eigh(sigma)
    d = w.shape[-1]
    mean = w.mean(axis=-1, keepdims=True)
    num = np.sum((w - mean) ** 2, axis=-1)
    den = np.sum(w**2, axis=-1)
    return np.sqrt((d / (d - 1.0)) * num / den)


def logfrob_distance(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Log-Frobenius (log-Euclidean) distance ``||log A - log B||_F``."""
    return np.linalg.norm(_logm_spd(a) - _logm_spd(b), axis=(-2, -1))


def airm_distance(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Affine-invariant Riemannian distance ``||log(A^-1/2 B A^-1/2)||_F``."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a_inv_half = _powm_spd(a, -0.5)
    middle = a_inv_half @ b @ a_inv_half
    return np.linalg.norm(_logm_spd(middle), axis=(-2, -1))


def _logm_spd(sigma: np.ndarray) -> np.ndarray:
    w, v = _eigh(np.asarray(sigma, dtype=float))
    return (v * np.log(w)[..., None, :]) @ np.swapaxes(v, -1, -2)


def _powm_spd(sigma: np.ndarray, p: float) -> np.ndarray:
    w, v = _eigh(np.asarray(sigma, dtype=float))
    return (v * (w**p)[..., None, :]) @ np.swapaxes(v, -1, -2)
