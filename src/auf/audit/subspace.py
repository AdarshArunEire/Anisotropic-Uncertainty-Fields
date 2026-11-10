"""Do errors concentrate in a covariance's top directions?

Capture ratio C_k(S, e): the fraction of |e|^2 lying in the span of the top-k
eigenvectors of S, divided by the fraction a random k-dimensional subspace
would capture, k / d with d = 3N - 6 (rigid-body directions removed from
both). C = 1 is random; C = 2 means the top-k directions carry twice their
share. Defined in notes/phaseA_preregistration.md.
"""
from __future__ import annotations

import numpy as np

__all__ = ["rigid_basis", "project_out_rigid", "capture_ratio"]


def rigid_basis(x: np.ndarray) -> np.ndarray:
    """Orthonormal (3N, 6) basis of infinitesimal translations and rotations
    of coordinates x (N, 3) about their centroid."""
    x = np.asarray(x, float)
    n = len(x)
    r = x - x.mean(0)
    cols = []
    for k in range(3):
        t = np.zeros((n, 3))
        t[:, k] = 1.0
        cols.append(t.ravel())
    for k in range(3):
        w = np.zeros(3)
        w[k] = 1.0
        cols.append(np.cross(w, r).ravel())
    q, _ = np.linalg.qr(np.stack(cols, 1))
    return q


def project_out_rigid(x: np.ndarray, vec_or_cov: np.ndarray) -> np.ndarray:
    """Remove rigid-body components from a 3N vector (or a batch, last axis
    3N) or from a 3N x 3N covariance."""
    q = rigid_basis(x)
    a = np.asarray(vec_or_cov, float)
    if a.ndim == 2 and a.shape[0] == a.shape[1] == q.shape[0]:
        p = np.eye(q.shape[0]) - q @ q.T
        return p @ a @ p
    return a - (a @ q) @ q.T


def capture_ratio(cov: np.ndarray, err: np.ndarray, ks=(5, 10, 20),
                  rigid_ref: np.ndarray | None = None) -> dict[int, np.ndarray]:
    """C_k for each k. ``err`` is (3N,) or (M, 3N); returns arrays (M,) or ().

    If ``rigid_ref`` (N, 3) is given, rigid directions at that structure are
    removed from both cov and err first.
    """
    err = np.atleast_2d(np.asarray(err, float))
    if rigid_ref is not None:
        cov = project_out_rigid(rigid_ref, cov)
        err = project_out_rigid(rigid_ref, err)
        d = cov.shape[0] - 6
    else:
        d = cov.shape[0]
    w, v = np.linalg.eigh(cov)
    v = v[:, ::-1]
    total = np.sum(err**2, axis=1)
    proj = err @ v
    out = {}
    for k in ks:
        frac = np.sum(proj[:, :k] ** 2, axis=1) / total
        out[k] = np.squeeze(frac / (k / d))
    return out
