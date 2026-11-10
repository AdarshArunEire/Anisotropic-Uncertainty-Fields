"""Phase A4: do AlphaFold's errors align with predicted directions? (+ A5 gate)

Primary set only (A2 AlphaFold error set, ``splits/afset_splits.tsv``).
Parameters are read from ``params/phasea_a3.json`` (frozen on train); nothing
is fitted here. The val-chosen PAE source goes to ``params/phasea_a4_choice.json``.

Per protein of the validation and test splits:
- capture ratio C_k (k = 5, 10, 20) of e under every source and control,
  core superposition (primary) and all-residue superposition (secondary);
- Gaussian log-likelihood per dimension with the frozen global scale;
- fusion (test only for the gate; val also run): prior mean = AF model,
  prior covariance = scale x source; data = Cα-Cα distances of the
  experimental structure, pairs 5-25 Å apart, 0.5 per residue, noise sd
  1.0 Å; 3 independent draws per protein, averaged. Outcomes: Cα RMSD
  (superposed) and lDDT-Cα gain over the AF model.

Then (``--summarize``) the PAE source is chosen on validation by geometric
mean C_10, and the pre-registered gate is evaluated on test:
(i)  geometric mean of C_10(source)/C_10(pLDDT scalar) with 95% bootstrap
     CI lower bound > 1.5 (amendment 1);
(ii) fusion gain over oracle scalar and PAE rigid-domain control, both RMSD
     and lDDT, mean paired difference CI excluding zero and better in >= 60%.

Usage:
  python3 experiments/phaseA/a4_alignment.py --workers 4      # compute
  python3 experiments/phaseA/a4_alignment.py --summarize      # gate + report
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
from auf.audit.subspace import capture_ratio  # noqa: E402
from auf.phasea import pae as P  # noqa: E402
from auf.phasea.sources import CONTROLS, PAE_SOURCES, SOURCES, build_sources, load_protein  # noqa: E402

KS = (5, 10, 20)
PAIRS_PER_RES, NOISE_SD, DMIN, DMAX, N_DRAWS = 0.5, 1.0, 5.0, 25.0, 3
OUT = Path(os.environ.get("AUF_A4_OUT", "results/phaseA"))


def _seed(name: str, salt: str) -> int:
    return int(hashlib.sha256(f"{salt}|{name}".encode()).hexdigest()[:8], 16)


def _random_cov(n3: int, rng) -> np.ndarray:
    g = rng.standard_normal((n3, n3))
    return g @ g.T


def fusion_draws(p: dict, rng) -> list[tuple[np.ndarray, np.ndarray]]:
    exp = p["exp"]
    dist = np.linalg.norm(exp[:, None] - exp[None], axis=-1)
    i, j = np.where(np.triu((dist >= DMIN) & (dist <= DMAX), 1))
    m = max(1, int(round(PAIRS_PER_RES * p["n"])))
    out = []
    for _ in range(N_DRAWS):
        k = rng.choice(len(i), size=min(m, len(i)), replace=False)
        pairs = np.stack([i[k], j[k]], 1)
        d_obs = dist[i[k], j[k]] + rng.normal(0, NOISE_SD, len(k))
        out.append((pairs, d_obs))
    return out


def _worker(args):
    path, split, params, do_fusion = args
    p = load_protein(path, "core")
    p_all = load_protein(path, "all")
    srcs = build_sources(p, params)
    rng = np.random.default_rng(_seed(p["name"], "random"))
    caps = {}
    n3 = 3 * p["n"]
    for mode, pp in (("core", p), ("all", p_all)):
        s2 = srcs if mode == "core" else build_sources(pp, params)
        s2 = dict(s2, random=_random_cov(n3, rng))
        for name, c in s2.items():
            cr = capture_ratio(c, pp["e"], KS, rigid_ref=pp["af"])
            caps[(mode, name)] = {k: float(cr[k]) for k in KS}
    z = P.nonrigid_basis(p["af"])
    ll = {}
    for name, c in srcs.items():
        s = params["scale"].get(name, 1.0)
        q, logdet, d = P.whitened_quad(s * c, p["e"], z=z)
        ll[name] = -0.5 * (q + logdet + d * np.log(2 * np.pi)) / d

    fus = {}
    if do_fusion:
        chosen = params.get("chosen_pae")
        priors = {n: params["scale"].get(n, 1.0) * c for n, c in srcs.items()}
        if chosen is not None:
            priors["random_rotation"] = P.random_rotation_cov(
                priors[chosen], p["af"], np.random.default_rng(_seed(p["name"], "rot")))
        draws = fusion_draws(p, np.random.default_rng(_seed(p["name"], "fusion")))
        r0 = P.rmsd(p["af"], p["exp"])
        l0 = P.lddt_ca(p["af"], p["exp"])
        for name, prior in priors.items():
            dr, dl = [], []
            for pairs, d_obs in draws:
                x, _ = fuse_distances(p["af"], prior, pairs, d_obs, NOISE_SD)
                dr.append(r0 - P.rmsd(x, p["exp"]))
                dl.append(P.lddt_ca(x, p["exp"]) - l0)
            fus[name] = (float(np.mean(dr)), float(np.mean(dl)))
        fus["_base"] = (r0, l0)
    return {"name": p["name"], "split": split, "n": p["n"], "used_core": p["used_core"],
            "mean_plddt": float(p["plddt"].mean()),
            "rmsd_af": P.rmsd(p["af"], p["exp"]), "caps": caps, "ll": ll, "fus": fus}


def compute(args):
    params = _load_params(args)
    with open(args.splits) as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["split"] in ("val", "test")]
    # Stage 1: capture on val chooses the PAE source; stage 2 runs everything
    # (fusion included) once the choice is fixed.
    if "chosen_pae" not in params:
        jobs = [(str(Path(args.afset) / f"{r['name']}.npz"), "val", params, False)
                for r in rows if r["split"] == "val"]
        with Pool(args.workers) as pool:
            val = pool.map(_worker, jobs, chunksize=1)
        gm = {s: float(np.exp(np.mean([np.log(v["caps"][("core", s)][10]) for v in val])))
              for s in PAE_SOURCES}
        chosen = max(gm, key=gm.get)
        print("val geometric-mean C_10:", {k: round(v, 2) for k, v in gm.items()},
              "-> chosen", chosen)
        choice = {"chosen_pae": chosen, "chosen_on_val_gm_c10": gm}
        Path(args.choice).write_text(json.dumps(choice, indent=1) + "\n")
        params.update(choice)

    jobs = [(str(Path(args.afset) / f"{r['name']}.npz"), r["split"], params, True)
            for r in rows]
    with Pool(args.workers) as pool:
        res = pool.map(_worker, jobs, chunksize=1)
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "a4_per_protein.json", "w") as fh:
        json.dump([{**r, "caps": {f"{m}|{s}": v for (m, s), v in r["caps"].items()}}
                   for r in res], fh)
    print(f"wrote {len(res)} proteins to {OUT / 'a4_per_protein.json'}")


def _boot(x, stat=np.mean, n=5000, seed=0):
    x = np.asarray(x, float)
    rng = np.random.default_rng(seed)
    b = [stat(x[rng.integers(0, len(x), len(x))]) for _ in range(n)]
    return float(stat(x)), float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))


def _gmean(x):
    return float(np.exp(np.mean(np.log(np.asarray(x, float)))))


def _load_params(args) -> dict:
    params = json.loads(Path(args.params).read_text())
    if Path(args.choice).exists():
        params.update(json.loads(Path(args.choice).read_text()))
    return params


def figures(res, chosen, split="test"):
    """Dot + 95% CI charts: capture ratio vs pLDDT scalar, and fusion gains."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    blue, orange, gray, ink2, grid, surf = ("#2a78d6", "#eb6834", "#8a8984", "#52514e",
                                            "#e4e3df", "#fcfcfb")
    plt.rcParams.update({
        "figure.facecolor": surf, "axes.facecolor": surf, "axes.edgecolor": grid,
        "axes.labelcolor": ink2, "xtick.color": ink2, "ytick.color": ink2,
        "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": grid, "grid.linewidth": 0.6})
    rs = [r for r in res if r["split"] == split]

    def color(name):
        return orange if name == chosen else gray if name in CONTROLS + [
            "random", "random_rotation", "plddt_scalar"] else blue

    def dotci(ax, names, stats, ref=None):
        for y, (n, (m, lo, hi)) in enumerate(zip(names, stats)):
            c = color(n)
            ax.plot([lo, hi], [y, y], color=c, lw=2, solid_capstyle="round")
            ax.plot(m, y, "o", ms=8, color=c, mec=surf, mew=2)
        ax.set_yticks(range(len(names)))
        ax.set_yticklabels([n + (" (val-chosen)" if n == chosen else "") for n in names])
        ax.invert_yaxis()
        ax.grid(axis="y", visible=False)
        for r in ref or []:
            ax.axvline(r[0], color=ink2, lw=1, ls=r[1], zorder=0)

    names = SOURCES + CONTROLS + ["random"]
    stats = [_boot([r["caps"][f"core|{n}"]["10"] / r["caps"]["core|plddt_scalar"]["10"]
                    for r in rs], _gmean) for n in names]
    order = np.argsort([-s[0] for s in stats])
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    dotci(ax, [names[i] for i in order], [stats[i] for i in order],
          ref=[(1.0, "-"), (1.5, "--")])
    ax.set_xscale("log")
    ax.set_xlabel("C_10 / C_10(pLDDT scalar), geometric mean over proteins (95% CI)")
    fig.suptitle(f"Share of AlphaFold error in each source's top 10 directions, relative to "
                 f"the pLDDT scalar\n({split}, n={len(rs)}; solid = pLDDT, dashed = gate 1.5)",
                 fontsize=9, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(OUT / f"a4_capture_{split}.png", dpi=200)
    plt.close(fig)

    fnames = [n for n in rs[0]["fus"] if n != "_base"]
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.6), sharey=True)
    for ax, k, lab, mult in ((axes[0], 0, "Cα RMSD removed (Å)", 1.0),
                             (axes[1], 1, "lDDT-Cα gained (x100)", 100.0)):
        st = [_boot([mult * r["fus"][n][k] for r in rs]) for n in fnames]
        if k == 0:
            order = np.argsort([-s[0] for s in st])
        dotci(ax, [fnames[i] for i in order], [st[i] for i in order], ref=[(0.0, "-")])
        ax.set_xlabel(lab)
    axes[1].tick_params(labelleft=False)
    fig.suptitle(f"Fusion with sparse distances: gain over the AlphaFold model "
                 f"({split}, n={len(rs)}, mean, 95% CI)", fontsize=9, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(OUT / f"a4_fusion_{split}.png", dpi=200)
    plt.close(fig)


def summarize(args):
    params = _load_params(args)
    res = json.loads((OUT / "a4_per_protein.json").read_text())
    chosen = params["chosen_pae"]
    lines = ["# Phase A4 results (generated by a4_alignment.py --summarize)", ""]
    gate = {}
    for split in ("val", "test"):
        rs = [r for r in res if r["split"] == split]
        lines += [f"## {split} split ({len(rs)} proteins)", "",
                  f"Median AF Cα RMSD (core-superposed, all residues): "
                  f"{np.median([r['rmsd_af'] for r in rs]):.2f} Å; median mean pLDDT "
                  f"{np.median([r['mean_plddt'] for r in rs]):.1f}.", "",
                  "### Capture ratio (geometric mean over proteins [95% CI])", "",
                  "| source | C_5 | C_10 | C_20 | C_10 / C_10(pLDDT scalar) | C_10, all-residue superposition | log-lik/dim |",
                  "|---|---|---|---|---|---|---|"]
        names = SOURCES + CONTROLS + ["random"]
        for s in names:
            cells = []
            for k in KS:
                m, lo, hi = _boot([r["caps"][f"core|{s}"][str(k)] for r in rs], _gmean)
                cells.append(f"{m:.1f} [{lo:.1f}, {hi:.1f}]")
            ratio = [r["caps"][f"core|{s}"]["10"] / r["caps"]["core|plddt_scalar"]["10"]
                     for r in rs]
            m, lo, hi = _boot(ratio, _gmean)
            if split == "test":
                gate[("i", s)] = (m, lo, hi)
            alls = _gmean([r["caps"][f"all|{s}"]["10"] for r in rs])
            llv = [r["ll"][s] for r in rs if s in r["ll"]]
            llc = f"{np.mean(llv):.3f}" if llv else "—"
            tag = " (chosen)" if s == chosen else ""
            lines.append(f"| {s}{tag} | " + " | ".join(cells)
                         + f" | {m:.2f} [{lo:.2f}, {hi:.2f}] | {alls:.1f} | {llc} |")
        lines += ["", "### Fusion: gain over the AlphaFold model (mean [95% CI])", "",
                  f"Pairs 5–25 Å, {PAIRS_PER_RES}/residue, noise {NOISE_SD} Å, "
                  f"{N_DRAWS} draws/protein. Positive = better.", "",
                  "| prior | ΔRMSD (Å) | ΔlDDT-Cα (x100) |", "|---|---|---|"]
        fnames = [n for n in rs[0]["fus"] if n != "_base"]
        for s in fnames:
            dr = [r["fus"][s][0] for r in rs]
            dl = [100 * r["fus"][s][1] for r in rs]
            a, b = _boot(dr), _boot(dl)
            lines.append(f"| {s} | {a[0]:+.3f} [{a[1]:+.3f}, {a[2]:+.3f}] | "
                         f"{b[0]:+.2f} [{b[1]:+.2f}, {b[2]:+.2f}] |")
        lines += ["", "### Paired fusion differences vs controls "
                  "(source − control; mean [95% CI], % proteins better)", "",
                  "| source | control | ΔΔRMSD (Å) | % better | ΔΔlDDT (x100) | % better |",
                  "|---|---|---|---|---|---|"]
        for s in [chosen, "anm", "anm_recal"] + [x for x in PAE_SOURCES if x != chosen]:
            for c in ("oracle_scalar", "pae_rigid_control", "pae_domains", "plddt_scalar",
                      "random_rotation"):
                if (c == "random_rotation" and s != chosen) or c == s:
                    continue
                dr = np.array([r["fus"][s][0] - r["fus"][c][0] for r in rs])
                dl = np.array([100 * (r["fus"][s][1] - r["fus"][c][1]) for r in rs])
                a, b = _boot(dr), _boot(dl)
                pr, pl = 100 * np.mean(dr > 0), 100 * np.mean(dl > 0)
                if split == "test":
                    gate[("ii", s, c)] = (a, pr, b, pl)
                lines.append(f"| {s} | {c} | {a[0]:+.3f} [{a[1]:+.3f}, {a[2]:+.3f}] | "
                             f"{pr:.0f}% | {b[0]:+.2f} [{b[1]:+.2f}, {b[2]:+.2f}] | {pl:.0f}% |")
        lines.append("")

    lines += ["## Gate (A5), test split, pre-registered criteria", ""]
    verdict = {}
    for s in (chosen, "anm"):
        m, lo, hi = gate[("i", s)]
        g1 = lo > 1.5
        g2 = True
        for c in ("oracle_scalar", "pae_rigid_control"):
            a, pr, b, pl = gate[("ii", s, c)]
            ok = a[1] > 0 and b[1] > 0 and pr >= 60 and pl >= 60
            g2 &= ok
            lines.append(f"- {s} vs {c}: RMSD {a[0]:+.3f} [{a[1]:+.3f}, {a[2]:+.3f}] Å "
                         f"({pr:.0f}% better), lDDT {b[0]:+.2f} [{b[1]:+.2f}, {b[2]:+.2f}] "
                         f"({pl:.0f}% better) -> {'pass' if ok else 'fail'}")
        verdict[s] = g1 and g2
        lines.append(f"- {s}: (i) ratio {m:.2f} [{lo:.2f}, {hi:.2f}] -> "
                     f"{'pass' if g1 else 'fail'}; (ii) -> {'pass' if g2 else 'fail'}")
        # Stricter control decided before any real data (not pre-registered):
        # the train-fitted rigid-domain recipe instead of the fixed broad one.
        if s != "pae_domains":
            a, pr, b, pl = gate[("ii", s, "pae_domains")]
            ok = a[1] > 0 and b[1] > 0 and pr >= 60 and pl >= 60
            lines.append(f"- (stricter, not in the gate) {s} vs fitted pae_domains: RMSD "
                         f"{a[0]:+.3f} [{a[1]:+.3f}, {a[2]:+.3f}] Å ({pr:.0f}%), lDDT "
                         f"{b[0]:+.2f} [{b[1]:+.2f}, {b[2]:+.2f}] ({pl:.0f}%) -> "
                         f"{'pass' if ok else 'fail'}")
    go = any(verdict.values())
    lines += ["", f"**Primary-set verdict: {'GO' if go else 'NO-GO on the primary set'}** "
              "(PARTIAL requires the ATLAS/MD secondary set, reported separately).", ""]
    for split in ("val", "test"):
        figures(res, chosen, split)
    lines += ["![capture](a4_capture_test.png)", "", "![fusion](a4_fusion_test.png)", ""]
    (OUT / "a4_report.md").write_text("\n".join(lines))
    print("\n".join(lines))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--afset", default="data/afset")
    ap.add_argument("--splits", default="splits/afset_splits.tsv")
    ap.add_argument("--params", default="params/phasea_a3.json")
    ap.add_argument("--choice", default="params/phasea_a4_choice.json",
                    help="PAE source chosen on validation (written by stage 1)")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--summarize", action="store_true")
    args = ap.parse_args()
    summarize(args) if args.summarize else compute(args)


if __name__ == "__main__":
    main()
