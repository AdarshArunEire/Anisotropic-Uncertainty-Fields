"""Choosing which residue-pair distances to measure (Phase B probe).

All choices use only the AlphaFold model, its pLDDT and a prior covariance;
never the experimental structure. See notes/phaseB_design_preregistration.md.
"""
from __future__ import annotations

import numpy as np

__all__ = ["candidate_pairs", "greedy_variance_reduction", "plddt_pairs", "random_pairs"]


def candidate_pairs(x: np.ndarray, seq: str | None = None, dmin: float = 5.0,
                    dmax: float = 25.0, min_sep: int = 3, crosslinkable: bool = False) -> np.ndarray:
    """(C, 2) residue pairs from model coordinates x (N, 3)."""
    n = len(x)
    d = np.linalg.norm(x[:, None] - x[None], axis=-1)
    i, j = np.triu_indices(n, min_sep)
    if crosslinkable:
        ok = np.array([c == "K" for c in seq])
        ok[0] = True  # N-terminal amine
        keep = ok[i] & ok[j] & (d[i, j] <= dmax)
    else:
        keep = (d[i, j] >= dmin) & (d[i, j] <= dmax)
    return np.stack([i[keep], j[keep]], 1)


def random_pairs(cands: np.ndarray, m: int, rng: np.random.Generator) -> np.ndarray:
    return cands[rng.choice(len(cands), size=min(m, len(cands)), replace=False)]


def plddt_pairs(cands: np.ndarray, plddt: np.ndarray, m: int) -> np.ndarray:
    """Greedy on (100 - pLDDT_i) + (100 - pLDDT_j); each residue used once."""
    score = (100 - plddt[cands[:, 0]]) + (100 - plddt[cands[:, 1]])
    used, out = set(), []
    for k in np.argsort(-score, kind="stable"):
        a, b = cands[k]
        if a in used or b in used:
            continue
        out.append(cands[k])
        used.update((a, b))
        if len(out) == m:
            break
    return np.array(out).reshape(-1, 2)


def greedy_variance_reduction(x: np.ndarray, prior: np.ndarray, cands: np.ndarray,
                              m: int, noise_sd: float) -> np.ndarray:
    """A-optimal greedy: repeatedly measure the pair that most reduces the
    trace of the (linearised) posterior covariance."""
    n = len(x)
    p = prior.copy()
    a, b = cands[:, 0], cands[:, 1]
    u = x[a] - x[b]
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    s2 = noise_sd**2
    chosen = []
    alive = np.ones(len(cands), bool)
    for _ in range(min(m, len(cands))):
        g = p.reshape(3 * n, n, 3)
        ph = np.einsum("dck,ck->dc", g[:, a, :], u) - np.einsum("dck,ck->dc", g[:, b, :], u)
        hph = _rows(ph, a, b, u)
        gain = np.sum(ph**2, axis=0) / (hph + s2)
        gain[~alive] = -np.inf
        k = int(np.argmax(gain))
        chosen.append(k)
        alive[k] = False
        v = ph[:, k]
        p = p - np.outer(v, v) / (hph[k] + s2)
    return cands[chosen]


def _rows(ph: np.ndarray, a: np.ndarray, b: np.ndarray, u: np.ndarray) -> np.ndarray:
    """hᵀ P h for every candidate, from the columns P h (3N, C)."""
    c = np.arange(len(a))
    pa = ph[(3 * a[:, None] + np.arange(3)).ravel(), np.repeat(c, 3)].reshape(-1, 3)
    pb = ph[(3 * b[:, None] + np.arange(3)).ravel(), np.repeat(c, 3)].reshape(-1, 3)
    return np.sum(u * (pa - pb), axis=1)
