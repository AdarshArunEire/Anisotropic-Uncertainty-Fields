"""Phase A3: fit the PAE -> covariance recipes on the TRAIN split only, then freeze.

Everything here sees only train proteins of ``splits/afset_splits.tsv``.
Error vectors use the primary superposition (pLDDT >= 70 core).

1. pLDDT -> per-coordinate error variance: pooled mean of |e_i|²/3 in pLDDT
   bins, made monotone (non-increasing in pLDDT) by pool-adjacent-violators.
2. For every source family and parameter grid point, the Gaussian
   log-likelihood of e (non-rigid subspace) with one global scale s
   (closed-form MLE s = Σ eᵀS⁻¹e / Σ d). Best grid point per family is kept.
   Grids are fixed here, before fitting.
3. Writes ``params/phasea_a3.json`` (committed; A4 reads it and never refits)
   and ``results/phaseA/a3_fit.tsv`` (all grid points).

Usage: python3 experiments/phaseA/a3_fit_recipes.py [--workers 4]
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from multiprocessing import Pool
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np  # noqa: E402

from auf.phasea import pae as P  # noqa: E402
from auf.phasea.sources import load_protein  # noqa: E402

BIN_EDGES = [0, 50, 60, 70, 75, 80, 85, 90, 95, 101]
GRID = (
    [("pae_reldisp",), ("anm",), ("plddt_scalar",)]
    + [("pae_domains", t, r) for t in (0.5, 1.0, 2.0, 4.0, 8.0, 16.0)
       for r in (0.02, 0.05, 0.1, 0.2, 0.4)]
    + [("pae_lowrank", k, w) for k in (5, 10, 20, 40, 80)
       for w in (0.03, 0.1, 0.3, 1.0, 3.0)]
    + [("anm_recal", pw) for pw in (0.25, 0.5, 0.75, 1.5, 2.0)]
)


def _pava_decreasing(y, w):
    """Weighted pool-adjacent-violators, non-increasing fit."""
    blocks = [[yi, wi, 1] for yi, wi in zip(y, w)]
    i = 0
    while i < len(blocks) - 1:
        if blocks[i][0] < blocks[i + 1][0]:
            a, b = blocks[i], blocks[i + 1]
            wt = a[1] + b[1]
            blocks[i] = [(a[0] * a[1] + b[0] * b[1]) / wt, wt, a[2] + b[2]]
            del blocks[i + 1]
            i = max(i - 1, 0)
        else:
            i += 1
    return [b[0] for b in blocks for _ in range(b[2])]


def fit_plddt_table(prots):
    pl = np.concatenate([p["plddt"] for p in prots])
    v = np.concatenate([np.sum(p["e"].reshape(-1, 3) ** 2, 1) / 3.0 for p in prots])
    centers, means, counts = [], [], []
    for lo, hi in zip(BIN_EDGES[:-1], BIN_EDGES[1:]):
        m = (pl >= lo) & (pl < hi)
        if m.sum() == 0:
            continue
        centers.append(float(pl[m].mean()))
        means.append(float(v[m].mean()))
        counts.append(int(m.sum()))
    mono = _pava_decreasing(means, counts)
    return {"centers": centers, "var": [float(x) for x in mono],
            "raw_var": means, "count": counts}


def _worker(args):
    path, table = args
    p = load_protein(path)
    x, e = p["af"], p["e"]
    z = P.nonrigid_basis(x)
    noise = P.diag_cov(P.plddt_variance(p["plddt"], table))
    a = P.pae_reldisp_cov(p["pae"])
    labels = P.pae_domains(p["pae"])
    hess_cache = {}
    out = []
    for g in GRID:
        fam = g[0]
        if fam == "pae_reldisp":
            c = a
        elif fam == "anm":
            c = P.anm_spectrum_cov(x, 1.0)
        elif fam == "plddt_scalar":
            c = noise
        elif fam == "pae_domains":
            c = P.rigid_domain_cov(x, labels, g[1], g[2]) + noise
        elif fam == "pae_lowrank":
            key = ("lr", g[1])
            if key not in hess_cache:
                hess_cache[key] = P.lowrank(a, g[1])
            c = g[2] * hess_cache[key] + noise
        elif fam == "anm_recal":
            c = P.anm_spectrum_cov(x, g[1])
        q, logdet, d = P.whitened_quad(c, e, z=z)
        out.append((g, q, logdet, d))
    return p["name"], out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--afset", default="data/afset")
    ap.add_argument("--splits", default="splits/afset_splits.tsv")
    ap.add_argument("--out", default="params/phasea_a3.json")
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()

    with open(args.splits) as fh:
        train = sorted(r["name"] for r in csv.DictReader(fh, delimiter="\t")
                       if r["split"] == "train")
    paths = [str(Path(args.afset) / f"{n}.npz") for n in train]
    prots = [load_protein(pth) for pth in paths]
    table = fit_plddt_table(prots)
    print(f"{len(prots)} train proteins; pLDDT table: "
          + ", ".join(f"{c:.0f}:{v:.2f}" for c, v in zip(table["centers"], table["var"])))

    with Pool(args.workers) as pool:
        res = pool.map(_worker, [(pth, table) for pth in paths], chunksize=1)

    agg = {}
    for _, rows in res:
        for g, q, logdet, d in rows:
            a = agg.setdefault(g, [0.0, 0.0, 0])
            a[0] += q
            a[1] += logdet
            a[2] += d
    table_rows = []
    for g, (q, logdet, d) in agg.items():
        s = q / d
        ll = -0.5 * (d * np.log(s) + logdet + d) / d  # per-dimension, constant dropped
        table_rows.append({"family": g[0], "params": list(g[1:]), "scale": s, "ll_per_dim": ll})
    best = {}
    for r in table_rows:
        if r["family"] not in best or r["ll_per_dim"] > best[r["family"]]["ll_per_dim"]:
            best[r["family"]] = r

    params = {
        "fitted_on": f"train split of {args.splits} ({len(train)} proteins), "
                     "core (pLDDT>=70) superposition",
        "plddt_table": table,
        "dom_sig_t": best["pae_domains"]["params"][0],
        "dom_sig_r": best["pae_domains"]["params"][1],
        "lowrank_r": best["pae_lowrank"]["params"][0],
        "lowrank_w": best["pae_lowrank"]["params"][1],
        "anm_power": best["anm_recal"]["params"][0],
        "scale": {fam: r["scale"] for fam, r in best.items()},
        "ll_per_dim": {fam: r["ll_per_dim"] for fam, r in best.items()},
        "grid": [list(g) for g in GRID],
    }
    # Controls: oracle scalar is exact by construction; the PAE rigid-domain
    # control uses the calibrated pLDDT noise and a broad fixed rigid prior.
    params["scale"]["oracle_scalar"] = 1.0
    params["scale"]["pae_rigid_control"] = 1.0
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(params, indent=1) + "\n")

    rd = Path("results/phaseA")
    rd.mkdir(parents=True, exist_ok=True)
    with open(rd / "a3_fit.tsv", "w") as fh:
        fh.write("family\tparams\tscale\tll_per_dim\n")
        for r in sorted(table_rows, key=lambda r: (r["family"], -r["ll_per_dim"])):
            fh.write(f"{r['family']}\t{r['params']}\t{r['scale']:.4g}\t{r['ll_per_dim']:.4f}\n")
    for fam, r in sorted(best.items(), key=lambda kv: -kv[1]["ll_per_dim"]):
        print(f"{fam:14s} params={r['params']} scale={r['scale']:.3g} "
              f"ll/dim={r['ll_per_dim']:.3f}")


if __name__ == "__main__":
    main()
