"""Directional metrics that separate *direction* from *amount* of anisotropy.

"Does Sigma beat isotropize(Sigma)?" conflates two things: whether the
predicted axes point the right way, and whether the predicted ellipsoid is
as elongated as the truth. A predictor with perfect axes but twice the
elongation can lose to the scalar. These metrics pull the two apart:

- ``direction_wins``: the strict test, Sigma vs its trace-matched surrogate.
- ``shape_power`` / ``optimal_shape_power``: the AIRM geodesic from the
  isotropic surrogate (alpha = 0) through Sigma (alpha = 1). The alpha that
  best matches the target says how much of the predicted anisotropy is real:
  0 means none (the scalar is optimal), 1 means all of it, > 1 means the
  predictor under-states it.
- ``principal_axis_cos``: |cos| between principal axes, 0.5 for random
  directions in 3D, 1 for perfect agreement. Pure direction, no magnitude.
- ``gaussian_loglik``: held-out average log-likelihood, the quantity a
  downstream consumer of the covariance actually pays for.
"""
from __future__ import annotations

import numpy as np

from auf.audit.spd import _powm_spd, airm_distance, isotropize

__all__ = [
    "unit_trace",
    "direction_wins",
    "shape_power",
    "optimal_shape_power",
    "principal_axis_cos",
    "gaussian_loglik",
]


def unit_trace(sigma: np.ndarray) -> np.ndarray:
    sigma = np.asarray(sigma, dtype=float)
    return sigma / np.trace(sigma, axis1=-2, axis2=-1)[..., None, None]


def direction_wins(pred: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Per residue: is pred closer to target than its isotropic surrogate?

    Both are compared at unit trace, so magnitude errors cannot decide it.
    """
    p, t = unit_trace(pred), unit_trace(target)
    return airm_distance(p, t) < airm_distance(isotropize(p), t)


def shape_power(sigma: np.ndarray, alpha: float) -> np.ndarray:
    """Unit-trace point at ``alpha`` on the geodesic from isotropic to Sigma."""
    return unit_trace(_powm_spd(unit_trace(sigma), alpha))


def optimal_shape_power(
    pred: np.ndarray,
    target: np.ndarray,
    alphas: np.ndarray = np.linspace(0.0, 2.0, 41),
) -> tuple[float, np.ndarray]:
    """alpha minimising median AIRM(shape_power(pred, alpha), target).

    Returns (best alpha, median distance for every alpha).
    """
    t = unit_trace(target)
    d = np.array([np.median(airm_distance(shape_power(pred, a), t)) for a in alphas])
    return float(alphas[np.argmin(d)]), d


def principal_axis_cos(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """|cos| between the largest-eigenvalue axes of a and b (random: 0.5)."""
    va = np.linalg.eigh(a)[1][..., -1]
    vb = np.linalg.eigh(b)[1][..., -1]
    return np.abs(np.sum(va * vb, axis=-1))


def gaussian_loglik(sigma: np.ndarray, sample_cov: np.ndarray) -> np.ndarray:
    """Average log N(x; mu, sigma) over data whose covariance about mu is
    ``sample_cov`` (ddof 0). Per residue, nats per frame."""
    d = sigma.shape[-1]
    _, logdet = np.linalg.slogdet(sigma)
    quad = np.trace(np.linalg.solve(sigma, sample_cov), axis1=-2, axis2=-1)
    return -0.5 * (quad + logdet + d * np.log(2 * np.pi))
