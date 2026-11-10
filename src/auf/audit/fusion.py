"""Downstream test: fusing sparse distance measurements with a prior covariance.

The use case where direction should matter most. An experiment (crosslinks,
FRET, NMR restraints, cryo-EM fitting) reports a few noisy inter-residue
distances; a prior covariance says how each residue is free to move. The
posterior structure is the prior conditioned on the data.

Distances are nonlinear in coordinates, so this is Gauss-Newton on the
linear-Gaussian model (an iterated Kalman update):

    d_ab = |x_a - x_b| + eps,  eps ~ N(0, s^2)
    x ~ N(mu, P),  P is 3N x 3N (block-diagonal for per-residue priors).
"""
from __future__ import annotations

import numpy as np

__all__ = ["block_diag_3", "fuse_distances", "posterior_blocks"]


def block_diag_3(blocks: np.ndarray) -> np.ndarray:
    """(N, 3, 3) per-residue blocks -> (3N, 3N) block-diagonal matrix."""
    n = len(blocks)
    out = np.zeros((n, 3, n, 3))
    idx = np.arange(n)
    out[idx, :, idx, :] = blocks
    return out.reshape(3 * n, 3 * n)


def _jacobian(x: np.ndarray, pairs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    a, b = pairs[:, 0], pairs[:, 1]
    diff = x[a] - x[b]
    dist = np.linalg.norm(diff, axis=-1)
    u = diff / dist[:, None]
    m, n = len(pairs), len(x)
    h = np.zeros((m, n, 3))
    h[np.arange(m), a] = u
    h[np.arange(m), b] = -u
    return dist, h.reshape(m, 3 * n)


def fuse_distances(
    mu: np.ndarray,
    prior: np.ndarray,
    pairs: np.ndarray,
    d_obs: np.ndarray,
    noise_sd: float,
    n_iter: int = 5,
) -> tuple[np.ndarray, np.ndarray]:
    """MAP structure and Laplace posterior covariance given distance data.

    Parameters
    ----------
    mu : (N, 3) prior mean.
    prior : (3N, 3N) prior covariance.
    pairs : (M, 2) residue index pairs measured.
    d_obs : (..., M) observed distances; a leading batch axis fuses several
        independent data sets against the same prior.
    noise_sd : measurement noise, Å.

    Returns
    -------
    x : (..., N, 3) posterior mean; cov : (..., 3N, 3N) posterior covariance.
    """
    d_obs = np.asarray(d_obs, float)
    batch = d_obs.shape[:-1]
    d_obs = d_obs.reshape(-1, d_obs.shape[-1])
    n = len(mu)
    r = noise_sd**2 * np.eye(len(pairs))
    xs, covs = [], []
    for d in d_obs:
        x = mu.copy()
        for _ in range(n_iter):
            dist, h = _jacobian(x, pairs)
            ph = prior @ h.T
            s = h @ ph + r
            # Gauss-Newton step for min |d - f(x)|^2_R + |x - mu|^2_P
            innov = d - dist + h @ (x - mu).reshape(-1)
            x = mu + (ph @ np.linalg.solve(s, innov)).reshape(n, 3)
        dist, h = _jacobian(x, pairs)
        ph = prior @ h.T
        cov = prior - ph @ np.linalg.solve(h @ ph + r, ph.T)
        xs.append(x)
        covs.append(cov)
    return (np.stack(xs).reshape(*batch, n, 3),
            np.stack(covs).reshape(*batch, 3 * n, 3 * n))


def posterior_blocks(cov: np.ndarray) -> np.ndarray:
    """(..., 3N, 3N) -> (..., N, 3, 3) per-residue diagonal blocks."""
    n = cov.shape[-1] // 3
    c = cov.reshape(*cov.shape[:-2], n, 3, n, 3)
    idx = np.arange(n)
    return np.moveaxis(c[..., idx, :, idx, :], 0, -3)
