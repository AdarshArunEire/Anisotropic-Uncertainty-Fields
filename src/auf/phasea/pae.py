"""PAE / pLDDT / ANM -> 3N x 3N Cα covariances (roadmap A3), plus helpers.

Conventions: coordinates (N, 3) in Å; covariances (3N, 3N) in Å², residue-major
(x1, y1, z1, x2, ...). ``pae[i, j]`` is AlphaFold's expected error of residue j
when the model is aligned on residue i.

Constructions (pre-registered in notes/phaseA_preregistration.md):

(a) ``pae_reldisp``: PAE²/3 as the per-coordinate variance of x_j − x_i.
    A matrix of pairwise "variances of differences" D is a variogram; the
    centred covariance with that variogram is C = −½ J D J (J = centring),
    clipped to PSD, then C ⊗ I₃.
(b) ``pae_domains``: PAE graph clusters (ChimeraX/Phenix convention: edges
    with weight 1/PAE for PAE < 5 Å, greedy modularity) as rigid bodies with
    translation sd σ_t and rotation sd σ_r, plus per-residue pLDDT noise.
(c) ``pae_lowrank``: top-r eigencomponents of (a), weight w, plus pLDDT noise.

Fitted parameters (σ_t, σ_r, r, w, the pLDDT variance table, global scales)
are fitted on the train split only, by ``experiments/phaseA/a3_fit_recipes.py``.
"""
from __future__ import annotations

import numpy as np

from auf.audit.subspace import rigid_basis
from auf.data.anm import anm_hessian

__all__ = [
    "kabsch_superpose", "error_vector", "sym_pae", "pae_reldisp_cov", "pae_domains",
    "rigid_domain_cov", "plddt_variance", "diag_cov", "lowrank", "anm_spectrum_cov",
    "random_rotation_cov", "nonrigid_basis", "whitened_quad", "lddt_ca", "rmsd",
]


def kabsch_superpose(mobile: np.ndarray, target: np.ndarray,
                     mask: np.ndarray | None = None) -> np.ndarray:
    """Rigidly move ``mobile`` (N, 3) onto ``target`` using residues in ``mask``."""
    m = np.ones(len(mobile), bool) if mask is None else np.asarray(mask, bool)
    a, b = mobile[m], target[m]
    ca, cb = a.mean(0), b.mean(0)
    u, _, vt = np.linalg.svd((a - ca).T @ (b - cb))
    d = np.sign(np.linalg.det(u @ vt))
    r = u @ np.diag([1.0, 1.0, d]) @ vt
    return (mobile - ca) @ r + cb


def rmsd(x: np.ndarray, ref: np.ndarray, superpose: bool = True) -> float:
    if superpose:
        x = kabsch_superpose(x, ref)
    return float(np.sqrt(np.mean(np.sum((x - ref) ** 2, axis=1))))


def error_vector(exp_xyz: np.ndarray, af_xyz: np.ndarray, plddt: np.ndarray,
                 mode: str = "core", core_cut: float = 70.0, min_core: int = 10):
    """AF model superposed on the experiment; returns (af_sup, e, used_core).

    e = x_exp − x_AF (3N) with the AF model's rigid-body directions removed.
    ``mode='core'`` superposes on pLDDT >= 70 residues (falls back to all if
    fewer than ``min_core``); ``mode='all'`` on all residues.
    """
    mask = plddt >= core_cut if mode == "core" else np.ones(len(plddt), bool)
    used_core = bool(mode == "core" and mask.sum() >= min_core)
    if not used_core:
        mask = np.ones(len(plddt), bool)
    af_sup = kabsch_superpose(af_xyz, exp_xyz, mask)
    e = (exp_xyz - af_sup).ravel()
    q = rigid_basis(af_sup)
    e = e - q @ (q.T @ e)
    return af_sup, e, used_core


def sym_pae(pae: np.ndarray) -> np.ndarray:
    pae = np.asarray(pae, float)
    return 0.5 * (pae + pae.T)


def _psd(c: np.ndarray) -> np.ndarray:
    w, v = np.linalg.eigh(0.5 * (c + c.T))
    return (v * np.clip(w, 0.0, None)) @ v.T


def pae_reldisp_cov(pae: np.ndarray, return_residue: bool = False):
    """(a): covariance whose Var(x_j − x_i) per coordinate is PAE_ij²/3."""
    d = sym_pae(pae) ** 2 / 3.0
    np.fill_diagonal(d, 0.0)
    n = len(d)
    j = np.eye(n) - 1.0 / n
    c = _psd(-0.5 * j @ d @ j)
    return c if return_residue else np.kron(c, np.eye(3))


def pae_domains(pae: np.ndarray, cutoff: float = 5.0, power: float = 1.0,
                resolution: float = 1.0) -> np.ndarray:
    """Residue -> domain label, ChimeraX ``alphafold pae`` clustering."""
    import networkx as nx

    p = np.asarray(pae, float)
    n = len(p)
    g = nx.Graph()
    g.add_nodes_from(range(n))
    w = 1.0 / np.maximum(p, 1e-3) ** power
    i, j = np.where(np.triu(p < cutoff, 1))
    g.add_weighted_edges_from(zip(i.tolist(), j.tolist(), w[i, j].tolist()))
    comms = nx.algorithms.community.greedy_modularity_communities(
        g, weight="weight", resolution=resolution)
    labels = np.empty(n, int)
    for k, c in enumerate(sorted(comms, key=lambda c: min(c))):
        labels[list(c)] = k
    return labels


def rigid_domain_cov(x: np.ndarray, labels: np.ndarray, sig_t: float,
                     sig_r: float) -> np.ndarray:
    """Each domain an independent rigid body: δx_i = t + ω × (x_i − c),
    t ~ N(0, σ_t² I), ω ~ N(0, σ_r² I) (radians)."""
    n = len(x)
    out = np.zeros((3 * n, 3 * n))
    for lab in np.unique(labels):
        idx = np.where(labels == lab)[0]
        r = x[idx] - x[idx].mean(0)
        jac = np.zeros((3 * len(idx), 6))
        for k in range(3):
            jac[k::3, k] = 1.0
            w = np.zeros(3)
            w[k] = 1.0
            jac[:, 3 + k] = np.cross(w, r).ravel()
        cov = jac @ np.diag([sig_t**2] * 3 + [sig_r**2] * 3) @ jac.T
        ii = (3 * idx[:, None] + np.arange(3)).ravel()
        out[np.ix_(ii, ii)] = cov
    return out


def plddt_variance(plddt: np.ndarray, table: dict) -> np.ndarray:
    """Per-coordinate variance (Å²) from pLDDT, piecewise-linear in the
    train-fitted table {'centers': [...], 'var': [...]}."""
    return np.interp(plddt, table["centers"], table["var"])


def diag_cov(var: np.ndarray) -> np.ndarray:
    return np.diag(np.repeat(np.asarray(var, float), 3))


def lowrank(c: np.ndarray, r: int) -> np.ndarray:
    w, v = np.linalg.eigh(c)
    w, v = w[::-1][:r], v[:, ::-1][:, :r]
    return (v * np.clip(w, 0, None)) @ v.T


def anm_spectrum_cov(x: np.ndarray, power: float = 1.0, cutoff: float = 15.0) -> np.ndarray:
    """ANM eigenvectors with eigenvalues λ^−power (power = 1 is plain ANM)."""
    h = anm_hessian(x, cutoff)
    w, v = np.linalg.eigh(h)
    w, v = w[6:], v[:, 6:]
    w = np.maximum(w, 1e-6 * w.max())
    return (v * w ** (-power)) @ v.T


def random_rotation_cov(c: np.ndarray, x: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Same spectrum as ``c`` (rigid modes removed), random eigenvectors in the
    non-rigid subspace: "any prior of this rank and spectrum" control."""
    q = rigid_basis(x)
    d = c.shape[0]
    p = np.eye(d) - q @ q.T
    w = np.linalg.eigvalsh(p @ c @ p)[::-1][: d - 6]
    z = rng.standard_normal((d, d - 6))
    z -= q @ (q.T @ z)
    basis, _ = np.linalg.qr(z)
    return (basis * w) @ basis.T


def nonrigid_basis(x: np.ndarray) -> np.ndarray:
    """Orthonormal (3N, 3N − 6) basis of the complement of the rigid modes at x."""
    q = rigid_basis(x)
    d = q.shape[0]
    return np.linalg.qr(np.concatenate([q, np.eye(d)], 1))[0][:, 6:d]


def whitened_quad(cov: np.ndarray, e: np.ndarray, x: np.ndarray | None = None,
                  floor: float = 1e-6, z: np.ndarray | None = None) -> tuple[float, float, int]:
    """(eᵀS⁻¹e, log det S, dim) in the non-rigid subspace at x.

    A relative floor (``floor`` × mean variance) keeps rank-deficient sources
    finite. Gaussian log-likelihood with scale s: −½(q/s + d log s + logdet + d log 2π).
    """
    z = nonrigid_basis(x) if z is None else z
    s = z.T @ cov @ z
    s += floor * np.trace(s) / len(s) * np.eye(len(s))
    ez = z.T @ e
    ch = np.linalg.cholesky(s)
    y = np.linalg.solve(ch, ez)
    return float(y @ y), float(2 * np.sum(np.log(np.diag(ch)))), len(s)


def lddt_ca(x: np.ndarray, ref: np.ndarray, r0: float = 15.0) -> float:
    """Global Cα lDDT of x against ref (no superposition), in [0, 1]."""
    dr = np.linalg.norm(ref[:, None] - ref[None], axis=-1)
    dx = np.linalg.norm(x[:, None] - x[None], axis=-1)
    m = (dr < r0) & ~np.eye(len(ref), dtype=bool)
    diff = np.abs(dr - dx)[m]
    return float(np.mean([(diff < t).mean() for t in (0.5, 1.0, 2.0, 4.0)]))
