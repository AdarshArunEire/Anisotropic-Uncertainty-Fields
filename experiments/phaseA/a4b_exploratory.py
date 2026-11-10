"""Phase A4b (POST HOC, exploratory; does not change the A5 gate).

Written after the A5 gate came out NO-GO, to locate the regime where coupled
priors could matter. Same frozen A3 parameters, same val-chosen source, no
refitting. Settings swept: distance noise 0.25 / 0.5 / 1.0 Å at 0.5 pairs per
residue; subgroups by the AF model's own error (core-superposed Cα RMSD
<= 1, 1-2, > 2 Å). Two independent draws per protein and setting.

Reported per split, setting and subgroup: mean and median gain over the AF
model (RMSD, lDDT x100) for each prior, and paired differences of the chosen
PAE source and ANM against the oracle scalar and the pLDDT scalar.

Usage: python3 experiments/phaseA/a4b_exploratory.py --workers 4
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from multiprocessing import Pool
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np  # noqa: E402

from auf.audit.fusion import fuse_distances  # noqa: E402
from auf.phasea import pae as P  # noqa: E402
from auf.phasea.sources import build_sources, load_protein  # noqa: E402

PRIORS = ["oracle_scalar", "plddt_scalar", "pae_lowrank", "pae_domains", "anm_recal"]
NOISES = (0.25, 0.5, 1.0)
DENSITY, N_DRAWS = 0.5, 2
GROUPS = (("AF RMSD <= 1 Å", 0, 1), ("1-2 Å", 1, 2), ("> 2 Å", 2, 1e9))
OUT = Path("results/phaseA")


def _worker(args):
    path, split, params = args
    p = load_protein(path)
    srcs = build_sources(p, params, PRIORS)
    exp = p["exp"]
    dist = np.linalg.norm(exp[:, None] - exp[None], axis=-1)
    i, j = np.where(np.triu((dist >= 5) & (dist <= 25), 1))
    seed = int(hashlib.sha256(f"a4b|{p['name']}".encode()).hexdigest()[:8], 16)
    rng = np.random.default_rng(seed)
    r0, l0 = P.rmsd(p["af"], exp), P.lddt_ca(p["af"], exp)
    out = {}
    for noise in NOISES:
        draws = []
        for _ in range(N_DRAWS):
            k = rng.choice(len(i), size=min(int(round(DENSITY * p["n"])), len(i)), replace=False)
            draws.append((np.stack([i[k], j[k]], 1), dist[i[k], j[k]] + rng.normal(0, noise, len(k))))
        for name, c in srcs.items():
            prior = params["scale"].get(name, 1.0) * c
            dr, dl = [], []
            for pairs, d_obs in draws:
                x, _ = fuse_distances(p["af"], prior, pairs, d_obs, noise)
                dr.append(r0 - P.rmsd(x, exp))
                dl.append(P.lddt_ca(x, exp) - l0)
            out[f"{noise}|{name}"] = (float(np.mean(dr)), float(np.mean(dl)))
    return {"name": p["name"], "split": split, "rmsd_af": r0, "lddt_af": l0, "fus": out}


def _boot(x, n=4000):
    x = np.asarray(x, float)
    rng = np.random.default_rng(0)
    b = x[rng.integers(0, len(x), (n, len(x)))].mean(1)
    return x.mean(), np.percentile(b, 2.5), np.percentile(b, 97.5)


def summarize(res, chosen):
    L = ["# Phase A4b: where could coupling matter? (POST HOC, exploratory)", "",
         "Not part of the pre-registered gate. Frozen A3 parameters; "
         f"chosen PAE source `{chosen}`. {DENSITY} pairs/residue, {N_DRAWS} draws.", ""]
    for split in ("val", "test"):
        rs = [r for r in res if r["split"] == split]
        for gname, lo, hi in GROUPS:
            g = [r for r in rs if lo < r["rmsd_af"] <= hi] if lo > 0 else \
                [r for r in rs if r["rmsd_af"] <= hi]
            if len(g) < 5:
                continue
            L += [f"## {split}, {gname} (n = {len(g)}; median AF RMSD "
                  f"{np.median([r['rmsd_af'] for r in g]):.2f} Å)", "",
                  "| noise | prior | ΔRMSD mean [CI] (Å) | ΔRMSD median | ΔlDDT mean [CI] (x100) | ΔlDDT median |",
                  "|---|---|---|---|---|---|"]
            for noise in NOISES:
                for name in PRIORS:
                    dr = [r["fus"][f"{noise}|{name}"][0] for r in g]
                    dl = [100 * r["fus"][f"{noise}|{name}"][1] for r in g]
                    a, b = _boot(dr), _boot(dl)
                    L.append(f"| {noise} | {name} | {a[0]:+.3f} [{a[1]:+.3f}, {a[2]:+.3f}] | "
                             f"{np.median(dr):+.3f} | {b[0]:+.2f} [{b[1]:+.2f}, {b[2]:+.2f}] | "
                             f"{np.median(dl):+.2f} |")
            L += ["", "Paired (source − control): mean ΔΔRMSD [CI], % better; "
                  "mean ΔΔlDDT [CI], % better", "",
                  "| noise | source | control | ΔΔRMSD (Å) | % | ΔΔlDDT (x100) | % |",
                  "|---|---|---|---|---|---|---|"]
            for noise in NOISES:
                for s in (chosen, "anm_recal"):
                    for c in ("oracle_scalar", "plddt_scalar"):
                        dr = np.array([r["fus"][f"{noise}|{s}"][0] - r["fus"][f"{noise}|{c}"][0] for r in g])
                        dl = np.array([100 * (r["fus"][f"{noise}|{s}"][1] - r["fus"][f"{noise}|{c}"][1]) for r in g])
                        a, b = _boot(dr), _boot(dl)
                        L.append(f"| {noise} | {s} | {c} | {a[0]:+.3f} [{a[1]:+.3f}, {a[2]:+.3f}] | "
                                 f"{100 * np.mean(dr > 0):.0f}% | {b[0]:+.2f} [{b[1]:+.2f}, {b[2]:+.2f}] | "
                                 f"{100 * np.mean(dl > 0):.0f}% |")
            L.append("")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--afset", default="data/afset")
    ap.add_argument("--splits", default="splits/afset_splits.tsv")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--summarize", action="store_true")
    args = ap.parse_args()
    params = json.loads(Path("params/phasea_a3.json").read_text())
    chosen = json.loads(Path("params/phasea_a4_choice.json").read_text())["chosen_pae"]
    raw = OUT / "a4b_per_protein.json"
    if not args.summarize:
        with open(args.splits) as fh:
            rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["split"] in ("val", "test")]
        jobs = [(str(Path(args.afset) / f"{r['name']}.npz"), r["split"], params) for r in rows]
        with Pool(args.workers) as pool:
            res = pool.map(_worker, jobs, chunksize=1)
        raw.write_text(json.dumps(res))
    res = json.loads(raw.read_text())
    text = summarize(res, chosen)
    (OUT / "a4b_exploratory.md").write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
