"""Phase B probe: does choosing which distances to measure from AlphaFold's
own uncertainty beat today's choices? (notes/phaseB_design_preregistration.md,
amendments 1-2)

Per protein (A2 set, frozen A3 parameters, val-chosen PAE source):
  strategies  random (3 draws), plddt_pairs, heuristic_deer, ensemble_mrmr,
              greedy_scalar, greedy_anm, greedy_pae
  budgets     m = 0 (chain only), 5, 10, 20, 50 (nested prefixes)
  fusion      prior = pae_lowrank (frozen scale) for every strategy, plus
              chain pseudo-observations (consecutive Ca-Ca at AF value,
              sd 0.05 Å); measured distances = experimental Ca-Ca + N(0, 1 Å),
              3 noise realisations. Secondary: greedy_scalar and greedy_anm
              also fused with their own priors.
  outcomes    lDDT-Ca gain (headline), Ca RMSD removed, realised squared
              error vs the linearised predicted tr(P_post) (greedy_pae).

Usage:
  python3 experiments/phaseB/b1_planning_probe.py --split val --limit 20   # debug
  python3 experiments/phaseB/b1_planning_probe.py --split test --workers 4
  python3 experiments/phaseB/b1_planning_probe.py --split test --summarize
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

from auf import design as D  # noqa: E402
from auf.audit.fusion import fuse_distances  # noqa: E402
from auf.audit.subspace import rigid_basis  # noqa: E402
from auf.phasea import pae as P  # noqa: E402
from auf.phasea.sources import build_sources, load_protein  # noqa: E402

BUDGETS = (0, 5, 10, 20, 50)
NOISE, CHAIN_SD, N_REAL, N_RANDOM, MAX_CANDS = 1.0, 0.05, 3, 3, 4000
STRATS = ["random", "plddt_pairs", "heuristic_deer", "ensemble_mrmr",
          "greedy_scalar", "greedy_anm", "greedy_pae"]
OUT = Path("results/phaseB")


def _seed(name, salt):
    return int(hashlib.sha256(f"b1|{salt}|{name}".encode()).hexdigest()[:8], 16)


def _nonrigid(cov, x):
    q = rigid_basis(x)
    pr = np.eye(len(q)) - q @ q.T
    return pr @ cov @ pr


def _sq_err(x, ref):
    return float(np.sum((P.kabsch_superpose(x, ref) - ref) ** 2))


def _worker(args):
    path, params, chosen = args
    p = load_protein(path)
    x, exp, n = p["af"], p["exp"], p["n"]
    srcs = build_sources(p, params, [chosen, "plddt_scalar", "anm_recal"], oracle=False)
    pri = {k: _nonrigid(params["scale"][k] * c, x) for k, c in srcs.items()}
    pri["greedy_pae"], pri["greedy_scalar"], pri["greedy_anm"] = (
        pri[chosen], pri["plddt_scalar"], pri["anm_recal"])
    rng = np.random.default_rng(_seed(p["name"], "cands"))
    cands = D.candidate_pairs(x)
    if len(cands) > MAX_CANDS:
        cands = cands[np.sort(rng.choice(len(cands), MAX_CANDS, replace=False))]
    mmax = max(BUDGETS)

    seqs = {
        "plddt_pairs": D.plddt_pairs(cands, p["plddt"], mmax),
        "heuristic_deer": D.heuristic_deer(x, mmax),
        "ensemble_mrmr": D.ensemble_mrmr(x, pri[chosen], cands, mmax, NOISE,
                                         np.random.default_rng(_seed(p["name"], "mrmr"))),
        "greedy_scalar": D.greedy_variance_reduction(x, pri["greedy_scalar"], cands, mmax, NOISE),
        "greedy_anm": D.greedy_variance_reduction(x, pri["greedy_anm"], cands, mmax, NOISE),
        "greedy_pae": D.greedy_variance_reduction(x, pri["greedy_pae"], cands, mmax, NOISE),
    }
    rr = np.random.default_rng(_seed(p["name"], "random"))
    for r in range(N_RANDOM):
        seqs[f"random#{r}"] = D.random_pairs(cands, mmax, rr)

    chain = np.stack([np.arange(n - 1), np.arange(1, n)], 1)
    chain_d = np.linalg.norm(x[1:] - x[:-1], axis=1)
    dexp = lambda pr: np.linalg.norm(exp[pr[:, 0]] - exp[pr[:, 1]], axis=1)
    l0, r0, e0 = P.lddt_ca(x, exp), P.rmsd(x, exp), _sq_err(x, exp)
    nrng = np.random.default_rng(_seed(p["name"], "noise"))
    noise_draws = nrng.normal(0, NOISE, (N_REAL, mmax))

    def fuse(prior, pairs, real):
        allp = np.concatenate([chain, pairs]) if len(pairs) else chain
        dobs = np.concatenate([chain_d, dexp(pairs) + noise_draws[real, :len(pairs)]]) \
            if len(pairs) else chain_d
        sd = np.concatenate([np.full(len(chain), CHAIN_SD), np.full(len(pairs), NOISE)])
        xf, _ = fuse_distances(x, prior, allp, dobs, sd)
        return xf

    res = {}
    jobs = [(s, s, chosen) for s in seqs] + [
        ("greedy_scalar@own", "greedy_scalar", "greedy_scalar"),
        ("greedy_anm@own", "greedy_anm", "greedy_anm")]
    for label, sname, prior_key in jobs:
        prior = pri[prior_key] if prior_key in pri else pri[chosen]
        for m in BUDGETS:
            if m == 0 and label != "greedy_pae":
                continue  # chain-only reference computed once
            pairs = seqs[sname][:m]
            out = []
            for real in range(N_REAL):
                xf = fuse(prior, pairs, real)
                out.append((P.lddt_ca(xf, exp) - l0, r0 - P.rmsd(xf, exp), _sq_err(xf, exp)))
            res[f"{label}|{m}"] = np.mean(out, axis=0).tolist()
    # Predicted curve (greedy_pae): linearised tr(P_post) incl. chain observations.
    pred = {}
    for m in BUDGETS:
        pairs = seqs["greedy_pae"][:m]
        allp = np.concatenate([chain, pairs]) if m else chain
        sd = np.concatenate([np.full(len(chain), CHAIN_SD), np.full(m, NOISE)])
        pred[m] = D.linear_posterior_trace(x, pri[chosen], allp, sd)
    return {"name": p["name"], "n": n, "rmsd_af": r0, "lddt_af": l0, "sqerr_af": e0,
            "n_cands": int(len(cands)), "n_heur": int(len(seqs["heuristic_deer"])),
            "res": res, "pred": pred}


def _boot(x, n=4000):
    x = np.asarray(x, float)
    rng = np.random.default_rng(0)
    b = x[rng.integers(0, len(x), (n, len(x)))].mean(1)
    return float(x.mean()), float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))


def _get(r, label, m, k):
    if label == "random":
        return float(np.mean([r["res"][f"random#{i}|{m}"][k] for i in range(N_RANDOM)]))
    key = f"{label}|{m}" if m else "greedy_pae|0"
    return r["res"][key][k]


def summarize(res, split, chosen):
    L = [f"# Phase B probe: measurement planning ({split}, {len(res)} proteins)", "",
         "Generated by `b1_planning_probe.py --summarize`. Pre-registration: "
         "[phaseB_design_preregistration.md](../../notes/phaseB_design_preregistration.md). "
         f"Fusion prior `{chosen}` + chain restraints; noise {NOISE} Å.", "",
         f"Median AF Cα RMSD {np.median([r['rmsd_af'] for r in res]):.2f} Å; median AF lDDT-Cα "
         f"{100 * np.median([r['lddt_af'] for r in res]):.1f}. Proteins with < 50 heuristic "
         f"pairs available: {sum(r['n_heur'] < 50 for r in res)}.", "",
         "## Gain over the AlphaFold model (mean [95% CI])", "",
         "| strategy | m | lDDT-Cα gain (x100) | Cα RMSD removed (Å) |", "|---|---|---|---|"]
    labels = STRATS + ["greedy_scalar@own", "greedy_anm@own"]
    a = _boot([100 * _get(r, "greedy_pae", 0, 0) for r in res])
    b = _boot([_get(r, "greedy_pae", 0, 1) for r in res])
    L.append(f"| chain only | 0 | {a[0]:+.2f} [{a[1]:+.2f}, {a[2]:+.2f}] | "
             f"{b[0]:+.3f} [{b[1]:+.3f}, {b[2]:+.3f}] |")
    for m in BUDGETS[1:]:
        for s in labels:
            a = _boot([100 * _get(r, s, m, 0) for r in res])
            b = _boot([_get(r, s, m, 1) for r in res])
            L.append(f"| {s} | {m} | {a[0]:+.2f} [{a[1]:+.2f}, {a[2]:+.2f}] | "
                     f"{b[0]:+.3f} [{b[1]:+.3f}, {b[2]:+.3f}] |")
    L += ["", "## Paired: greedy_pae − other (m = 20 primary; mean [95% CI], % proteins better)", "",
          "| m | vs | ΔlDDT (x100) | % better | ΔRMSD removed (Å) | % better |",
          "|---|---|---|---|---|---|"]
    gate = {}
    for m in BUDGETS[1:]:
        for s in [x for x in STRATS if x != "greedy_pae"]:
            dl = np.array([100 * (_get(r, "greedy_pae", m, 0) - _get(r, s, m, 0)) for r in res])
            dr = np.array([_get(r, "greedy_pae", m, 1) - _get(r, s, m, 1) for r in res])
            a, b = _boot(dl), _boot(dr)
            pl, pr = 100 * np.mean(dl > 0), 100 * np.mean(dr > 0)
            if m == 20:
                gate[s] = (a, pl)
            L.append(f"| {m} | {s} | {a[0]:+.2f} [{a[1]:+.2f}, {a[2]:+.2f}] | {pl:.0f}% | "
                     f"{b[0]:+.3f} [{b[1]:+.3f}, {b[2]:+.3f}] | {pr:.0f}% |")
    L += ["", "## Pre-registered criterion (m = 20, lDDT-Cα)", ""]
    ok_all = True
    for s in ("random", "greedy_scalar", "heuristic_deer"):
        a, pl = gate[s]
        ok = a[1] > 0 and pl >= 60
        ok_all &= ok
        L.append(f"- greedy_pae vs {s}: {a[0]:+.2f} [{a[1]:+.2f}, {a[2]:+.2f}] lDDT points, "
                 f"better in {pl:.0f}% -> {'pass' if ok else 'fail'}")
    a, pl = gate["ensemble_mrmr"]
    L.append(f"- (reported, not required) vs ensemble_mrmr: {a[0]:+.2f} [{a[1]:+.2f}, {a[2]:+.2f}], "
             f"better in {pl:.0f}%")
    L += ["", f"**Verdict: {'PASS' if ok_all else 'FAIL'}**", "",
          "## Predicted vs realised error (greedy_pae)", "",
          "Predicted = linearised tr(P_post) with the frozen scale; realised = squared Cα "
          "error after fusion (superposed), mean over noise realisations.", "",
          "| m | Spearman ρ across proteins | median realised / predicted |", "|---|---|---|"]
    from scipy.stats import spearmanr
    for m in BUDGETS:
        pred = np.array([r["pred"][str(m)] for r in res])
        real = np.array([_get(r, "greedy_pae", m, 2) for r in res])
        L.append(f"| {m} | {spearmanr(pred, real)[0]:.2f} | {np.median(real / pred):.2f} |")
    L.append("")
    return "\n".join(L), ok_all


def figure(res, split):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    surf, grid, ink2 = "#fcfcfb", "#e4e3df", "#52514e"
    plt.rcParams.update({"figure.facecolor": surf, "axes.facecolor": surf, "axes.edgecolor": grid,
                         "axes.labelcolor": ink2, "xtick.color": ink2, "ytick.color": ink2,
                         "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.color": grid, "grid.linewidth": 0.6})
    style = {"greedy_pae": ("#eb6834", "o"), "random": ("#8a8984", "s"),
             "heuristic_deer": ("#2a78d6", "^"), "ensemble_mrmr": ("#1baf7a", "D"),
             "greedy_scalar": ("#52514e", "v"), "plddt_pairs": ("#b3b2ad", "x")}
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    for s, (c, mk) in style.items():
        ys = [100 * np.mean([_get(r, s if m else "greedy_pae", m, 0) for r in res]) for m in BUDGETS]
        ax.plot(BUDGETS, ys, color=c, marker=mk, ms=6, lw=2, label=s)
    ax.set_xlabel("measured pairs (m)")
    ax.set_ylabel("lDDT-Cα gain over AlphaFold (points)")
    ax.legend(frameon=False, fontsize=8)
    fig.suptitle(f"Measurement planning: lDDT gained vs number of pairs ({split}, n={len(res)})",
                 fontsize=9, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(OUT / f"b1_curve_{split}.png", dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--split", default="test")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--summarize", action="store_true")
    args = ap.parse_args()
    params = json.loads(Path("params/phasea_a3.json").read_text())
    chosen = json.loads(Path("params/phasea_a4_choice.json").read_text())["chosen_pae"]
    OUT.mkdir(parents=True, exist_ok=True)
    tag = f"{args.split}" + (f"_dev{args.limit}" if args.limit else "")
    raw = OUT / f"b1_per_protein_{tag}.json"
    if not args.summarize:
        with open("splits/afset_splits.tsv") as fh:
            names = sorted(r["name"] for r in csv.DictReader(fh, delimiter="\t")
                           if r["split"] == args.split)
        if args.limit:
            names = names[: args.limit]
        jobs = [(f"data/afset/{n}.npz", params, chosen) for n in names]
        with Pool(args.workers) as pool:
            res = pool.map(_worker, jobs, chunksize=1)
        raw.write_text(json.dumps(res))
    res = json.loads(raw.read_text())
    text, _ = summarize(res, tag, chosen)
    (OUT / f"b1_report_{tag}.md").write_text(text)
    figure(res, tag)
    print(text)


if __name__ == "__main__":
    main()
