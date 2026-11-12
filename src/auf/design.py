"""Choosing which residue-pair distances to measure (Phase B probe).

All choices use only the AlphaFold model, its pLDDT and a prior covariance;
never the experimental structure. See notes/phaseB_design_preregistration.md.
"""
from __future__ import annotations

import numpy as np

__all__ = ["candidate_pairs", "greedy_variance_reduction", "plddt_pairs", "random_pairs",
           "ca_secondary_structure", "ca_exposed", "heuristic_deer", "ensemble_mrmr",
           "linear_posterior_trace"]


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


def ca_secondary_structure(x: np.ndarray) -> np.ndarray:
    """Per-residue 'H' / 'E' / 'C' from Cα geometry only (amendment 2):
    helix if covered by a window i..i+4 with d(i,i+3) < 6.0 and d(i,i+4) < 6.8 Å;
    strand if d(i-1,i+1) > 6.2 Å and not helix; runs shorter than 3 become 'C'."""
    n = len(x)
    d = lambda i, j: np.linalg.norm(x[i] - x[j])
    ss = np.array(["C"] * n)
    for i in range(n - 4):
        if d(i, i + 3) < 6.0 and d(i, i + 4) < 6.8:
            ss[i:i + 5] = "H"
    for i in range(1, n - 1):
        if ss[i] != "H" and d(i - 1, i + 1) > 6.2:
            ss[i] = "E"
    out, i = ss.copy(), 0
    while i < n:
        j = i
        while j < n and ss[j] == ss[i]:
            j += 1
        if ss[i] != "C" and j - i < 3:
            out[i:j] = "C"
        i = j
    return out


def ca_exposed(x: np.ndarray, radius: float = 10.0) -> np.ndarray:
    """Exposed if the Cα has fewer neighbours within ``radius`` than the median."""
    d = np.linalg.norm(x[:, None] - x[None], axis=-1)
    cnt = (d < radius).sum(1) - 1
    return cnt < np.median(cnt)


def _elements(ss: np.ndarray) -> np.ndarray:
    """Element id per residue for H/E runs, -1 for coil."""
    el, k, i, n = -np.ones(len(ss), int), 0, 0, len(ss)
    while i < n:
        j = i
        while j < n and ss[j] == ss[i]:
            j += 1
        if ss[i] != "C":
            el[i:j] = k
            k += 1
        i = j
    return el


def heuristic_deer(x: np.ndarray, m: int, dmin: float = 15.0, dmax: float = 60.0,
                   min_sep: int = 3) -> np.ndarray:
    """Kazmier 2011 / DEERFold-style rule (amendment 2, item 3)."""
    ss = ca_secondary_structure(x)
    el = _elements(ss)
    ok = (el >= 0) & ca_exposed(x)
    d = np.linalg.norm(x[:, None] - x[None], axis=-1)
    i, j = np.triu_indices(len(x), min_sep)
    keep = ok[i] & ok[j] & (d[i, j] > dmin) & (d[i, j] <= dmax)
    i, j = i[keep], j[keep]
    order = np.argsort(-np.log(j - i), kind="stable")
    used, out, taken = set(), [], np.zeros(len(i), bool)
    for k in order:
        key = tuple(sorted((el[i[k]], el[j[k]])))
        if key in used:
            continue
        used.add(key)
        out.append((i[k], j[k]))
        taken[k] = True
        if len(out) == m:
            return np.array(out)
    for k in order:  # all element pairs used: fill by separation alone
        if len(out) == m:
            break
        if not taken[k]:
            out.append((i[k], j[k]))
    return np.array(out, int).reshape(-1, 2)


def ensemble_mrmr(x: np.ndarray, prior: np.ndarray, cands: np.ndarray, m: int,
                  noise_sd: float, rng: np.random.Generator, n_samples: int = 200) -> np.ndarray:
    """Hays 2018 / FPS 2.0-style mRMR on an ensemble sampled from the prior
    (amendment 2, item 4). Gaussian MI: relevance 1/2 log(1 + v / s^2),
    redundancy mean of -1/2 log(1 - rho^2) over chosen pairs."""
    w, v = np.linalg.eigh(0.5 * (prior + prior.T))
    l = v * np.sqrt(np.clip(w, 0, None))
    xs = x[None] + (rng.standard_normal((n_samples, len(w))) @ l.T).reshape(n_samples, -1, 3)
    dist = np.linalg.norm(xs[:, cands[:, 0]] - xs[:, cands[:, 1]], axis=-1)
    var = dist.var(0)
    z = (dist - dist.mean(0)) / np.sqrt(var + 1e-12)
    relevance = 0.5 * np.log1p(var / noise_sd**2)
    red_sum = np.zeros(len(cands))
    chosen = []
    for t in range(min(m, len(cands))):
        score = relevance - (red_sum / t if t else 0.0)
        score[chosen] = -np.inf
        k = int(np.argmax(score))
        chosen.append(k)
        rho = np.clip(z.T @ z[:, k] / n_samples, -0.999999, 0.999999)
        red_sum += -0.5 * np.log1p(-rho**2)
    return cands[chosen]


def linear_posterior_trace(x: np.ndarray, prior: np.ndarray, pairs: np.ndarray,
                           noise_sd) -> float:
    """tr(P_post) of the linearised distance update at x (Å², total)."""
    if len(pairs) == 0:
        return float(np.trace(prior))
    n = len(x)
    a, b = pairs[:, 0], pairs[:, 1]
    u = x[a] - x[b]
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    h = np.zeros((len(pairs), n, 3))
    h[np.arange(len(pairs)), a] = u
    h[np.arange(len(pairs)), b] = -u
    h = h.reshape(len(pairs), 3 * n)
    ph = prior @ h.T
    r = np.diag(np.broadcast_to(np.asarray(noise_sd, float) ** 2, (len(pairs),)).copy())
    s = h @ ph + r
    return float(np.trace(prior) - np.sum(ph * np.linalg.solve(s, ph.T).T))
