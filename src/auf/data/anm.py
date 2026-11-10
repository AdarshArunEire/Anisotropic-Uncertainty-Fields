"""Anisotropic network model (ANM): a structure-only covariance baseline.

The cheapest physics that predicts *direction*: Cα atoms joined by identical
springs within a cutoff (Atilgan et al. 2001). One global parameter (spring
constant) sets magnitude, so with per-residue magnitude matched to the truth
it is a pure test of whether direction is predictable from structure alone,
with no learning. Any learned covariance model should be measured against
this, not only against the isotropic surrogate.
"""
from __future__ import annotations

import numpy as np

__all__ = ["anm_hessian", "anm_covariance"]


def anm_hessian(ca: np.ndarray, cutoff: float = 15.0, gamma: float = 1.0) -> np.ndarray:
    """3N x 3N ANM Hessian for Cα coordinates (N, 3), Å."""
    ca = np.asarray(ca, dtype=float)
    n = len(ca)
    diff = ca[None, :, :] - ca[:, None, :]  # r_j - r_i
    d2 = np.sum(diff**2, axis=-1)
    contact = (d2 < cutoff**2) & ~np.eye(n, dtype=bool)
    with np.errstate(divide="ignore", invalid="ignore"):
        block = -gamma * diff[..., :, None] * diff[..., None, :] / d2[..., None, None]
    block[~contact] = 0.0
    h = block.copy()
    idx = np.arange(n)
    h[idx, idx] = -block.sum(axis=1)
    return h.transpose(0, 2, 1, 3).reshape(3 * n, 3 * n)


def anm_covariance(ca: np.ndarray, cutoff: float = 15.0, n_zero: int = 6) -> np.ndarray:
    """Per-residue 3x3 ANM covariance blocks, (N, 3, 3), arbitrary units.

    Pseudo-inverse of the Hessian with the six rigid-body modes removed.
    """
    h = anm_hessian(ca, cutoff)
    w, v = np.linalg.eigh(h)
    w, v = w[n_zero:], v[:, n_zero:]
    if np.any(w <= 1e-8):
        raise ValueError("ANM network has extra zero modes; increase cutoff")
    n = len(ca)
    vr = v.reshape(n, 3, -1)
    return np.einsum("aik,ajk,k->aij", vr, vr, 1.0 / w)
