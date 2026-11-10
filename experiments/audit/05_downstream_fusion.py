"""Stage 1, downstream test: does anisotropy matter when fusing experimental data?

Scenario. A predictor supplies a structure (prior mean) and an uncertainty
(prior covariance). An experiment then measures a few noisy inter-residue
distances (crosslinks, FRET, NMR restraints). The structure is updated by
conditioning the prior on the data (``auf.audit.fusion``). Better
uncertainty should give a better updated structure and honest error bars.

Truth is real MD: frames of replica j. The prior comes from replica i (its
mean structure and covariance), so the prior is a *prediction*, not an
oracle. Priors compared, all with the same per-residue magnitude:

  scalar          per-residue sphere                       (what pLDDT-style gives)
  ellipsoid       per-residue 3x3 from replica i           (what DynaProt/BackFlip-2 predict)
  ellipsoid_shrunk  same, shrunk to alpha = 0.74 (Stage-1 audit value)
  anm_ellipsoid   per-residue ANM shape, no learning
  joint_anm       full 3N x 3N ANM (residues move together), one global scale
  joint_md        full 3N x 3N from replica i (upper reference)

Reported per setting (measurement density x noise), mean over proteins with
bootstrap CIs: RMSE of the updated structure to the true frame, the share of
the prior error removed, and calibration of the posterior (mean Mahalanobis
/ 3; 1 is honest, > 1 overconfident) and 90% coverage.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
from multiprocessing import Pool
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

PRIORS = ["scalar", "ellipsoid", "ellipsoid_shrunk", "anm_ellipsoid", "joint_anm", "joint_md"]
ALPHA = 0.74
SETTINGS = [(0.1, 0.5), (0.1, 2.0), (0.5, 0.5), (0.5, 2.0)]  # (pairs per residue, noise Å)
REPLICA_PAIRS = [(0, 1), (2, 0)]


def build_priors(train: np.ndarray, start_ca: np.ndarray) -> tuple[np.ndarray, dict]:
    from auf.audit import isotropize
    from auf.audit.direction import shape_power, unit_trace
    from auf.audit.fusion import block_diag_3
    from auf.data.anm import anm_covariance_full
    from auf.data.atlas import superpose

    t, n, _ = train.shape
    mu = train.mean(0)
    dev = (train - mu).reshape(t, 3 * n)
    full = dev.T @ dev / (t - 1) + 1e-3 * np.eye(3 * n)
    idx = np.arange(n)
    blocks = full.reshape(n, 3, n, 3)[idx, :, idx, :]
    tr = np.trace(blocks, axis1=-2, axis2=-1)[:, None, None]

    start_in, _ = superpose(start_ca[None], reference=mu)
    anm = anm_covariance_full(start_in[0])
    anm_blocks = anm.reshape(n, 3, n, 3)[idx, :, idx, :]
    anm_scale = np.trace(full) / np.trace(anm)
    return mu, {
        "scalar": block_diag_3(isotropize(blocks)),
        "ellipsoid": block_diag_3(blocks),
        "ellipsoid_shrunk": block_diag_3(shape_power(blocks, ALPHA) * tr),
        "anm_ellipsoid": block_diag_3(unit_trace(anm_blocks) * tr),
        "joint_anm": anm * anm_scale + 1e-3 * np.eye(3 * n),
        "joint_md": full,
    }


def run_one(args) -> dict | None:
    pdb_chain, raw_dir, n_frames, seed, keep_raw = args
    import mdtraj as md

    from auf.audit.fusion import fuse_distances, posterior_blocks
    from auf.data.atlas import fetch_atlas, read_ca_trajectories, superpose

    rng = np.random.default_rng(seed)
    try:
        entry = fetch_atlas(pdb_chain, raw_dir)
        ca = read_ca_trajectories(entry, pdb_chain)
        start = md.load(str(entry / f"{pdb_chain}.pdb"))
        start_ca = start.xyz[0, start.top.select("name CA")].astype(float) * 10.0
    except Exception as exc:
        return {"pdb_chain": pdb_chain, "error": repr(exc)}
    n_rep, t, n, _ = ca.shape
    aligned, _ = superpose(ca.reshape(-1, n, 3))
    aligned = aligned.reshape(n_rep, t, n, 3)
    c90 = stats.chi2.ppf(0.90, 3)

    out = {"pdb_chain": pdb_chain, "n_res": n, "results": []}
    for i, j in REPLICA_PAIRS:
        mu, priors = build_priors(aligned[i], start_ca)
        frames = aligned[j][np.linspace(0, t - 1, n_frames).astype(int)]
        base_rmse = np.sqrt(np.mean(np.sum((frames - mu) ** 2, -1), -1))
        dmat = np.linalg.norm(mu[:, None] - mu[None], axis=-1)
        cand = np.argwhere(np.triu((dmat > 5) & (dmat < 25), 1))
        for density, noise in SETTINGS:
            m = max(2, int(round(density * n)))
            pairs = cand[rng.choice(len(cand), size=min(m, len(cand)), replace=False)]
            true_d = np.linalg.norm(frames[:, pairs[:, 0]] - frames[:, pairs[:, 1]], axis=-1)
            d_obs = true_d + rng.normal(scale=noise, size=true_d.shape)
            measured = np.zeros(n, bool)
            measured[pairs.ravel()] = True
            for name in PRIORS:
                x, cov = fuse_distances(mu, priors[name], pairs, d_obs, noise)
                err = x - frames
                rmse = np.sqrt(np.mean(np.sum(err**2, -1), -1))
                rmse_meas = np.sqrt(np.mean(np.sum(err[:, measured] ** 2, -1), -1))
                blocks = posterior_blocks(cov)
                m2 = np.einsum("fri,frij,frj->fr", err, np.linalg.inv(blocks), err)
                out["results"].append({
                    "pair": f"{i}->{j}", "density": density, "noise": noise, "prior": name,
                    "rmse_prior_mean": float(base_rmse.mean()),
                    "rmse": float(rmse.mean()), "rmse_measured": float(rmse_meas.mean()),
                    "error_removed": float(np.mean(1 - rmse / base_rmse)),
                    "calib_scale": float(m2.mean() / 3), "coverage90": float((m2 <= c90).mean()),
                })
    if not keep_raw:
        shutil.rmtree(entry, ignore_errors=True)
    return out


def boot(x, rng, n_boot=4000):
    x = np.asarray(x)
    b = x[rng.integers(0, len(x), (n_boot, len(x)))].mean(1)
    return x.mean(), np.percentile(b, 2.5), np.percentile(b, 97.5)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--entries", default="data/audit/entries.txt")
    p.add_argument("--n-proteins", type=int, default=40)
    p.add_argument("--max-len", type=int, default=250)
    p.add_argument("--frames", type=int, default=40)
    p.add_argument("--raw-dir", default="data/atlas")
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--keep-raw", action="store_true")
    p.add_argument("--out", default="results/05_downstream_fusion")
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    entries = Path(args.entries).read_text().split()
    if args.max_len:
        lengths = {}
        for e in entries:
            f = Path("data/audit") / f"{e}.npz"
            if f.exists():
                with np.load(f) as z:
                    lengths[e] = len(z["reference"])
        entries = [e for e in entries if lengths.get(e, 0) <= args.max_len]
    entries = entries[: args.n_proteins]
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    jobs = [(e, args.raw_dir, args.frames, args.seed + k, args.keep_raw) for k, e in enumerate(entries)]
    per = []
    with Pool(args.workers) as pool:
        for k, r in enumerate(pool.imap_unordered(run_one, jobs), 1):
            print(f"[{k}/{len(jobs)}] {r['pdb_chain']}: {'ERROR ' + r['error'] if 'error' in r else 'ok'}",
                  flush=True)
            if "error" not in r:
                per.append(r)
    (out / "per_protein.json").write_text(json.dumps(per))

    rng = np.random.default_rng(args.seed)
    lines = [f"# Downstream fusion: {len(per)} proteins "
             f"({sum(r['n_res'] for r in per)} residues), {args.frames} true frames x "
             f"{len(REPLICA_PAIRS)} replica pairs each",
             "Mean over proteins [95% bootstrap CI].", ""]
    for density, noise in SETTINGS:
        lines.append(f"## {density:g} distances per residue, noise {noise:g} Å")
        lines.append(f"{'prior':<18}{'RMSE (Å)':>22}{'error removed':>24}{'calib scale':>14}{'cover90':>9}")
        base = None
        for name in PRIORS:
            rows = {}
            for r in per:
                v = [x for x in r["results"] if x["prior"] == name and x["density"] == density
                     and x["noise"] == noise]
                rows[r["pdb_chain"]] = {k: np.mean([x[k] for x in v]) for k in
                                        ("rmse", "error_removed", "calib_scale", "coverage90", "rmse_prior_mean")}
            vals = list(rows.values())
            if base is None:
                base = np.mean([v["rmse_prior_mean"] for v in vals])
                lines.append(f"{'(no data)':<18}{base:>10.3f}")
            rm = boot([v["rmse"] for v in vals], rng)
            er = boot([v["error_removed"] for v in vals], rng)
            lines.append(
                f"{name:<18}{rm[0]:>8.3f} [{rm[1]:.3f},{rm[2]:.3f}]"
                f"{100 * er[0]:>9.1f}% [{100 * er[1]:.1f},{100 * er[2]:.1f}]"
                f"{np.mean([v['calib_scale'] for v in vals]):>12.2f}"
                f"{100 * np.mean([v['coverage90'] for v in vals]):>8.0f}%")
        # paired comparisons that answer the question directly
        def paired(a, b):
            d = []
            for r in per:
                def get(name):
                    return np.mean([x["rmse"] for x in r["results"] if x["prior"] == name
                                    and x["density"] == density and x["noise"] == noise])
                d.append(get(b) - get(a))
            m, lo, hi = boot(d, rng)
            wins = np.mean(np.array(d) > 0)
            return f"{m:+.3f} Å [{lo:+.3f}, {hi:+.3f}], better in {100 * wins:.0f}% of proteins"
        lines.append(f"  ellipsoid vs scalar (RMSE gain):        {paired('ellipsoid', 'scalar')}")
        lines.append(f"  shrunk ellipsoid vs scalar:             {paired('ellipsoid_shrunk', 'scalar')}")
        lines.append(f"  joint ANM vs best per-residue ellipsoid: {paired('joint_anm', 'ellipsoid_shrunk')}")
        lines.append(f"  joint MD vs joint ANM:                  {paired('joint_md', 'joint_anm')}")
        lines.append("")
    report = "\n".join(lines)
    print(report)
    (out / "report.md").write_text(report + "\n")


if __name__ == "__main__":
    main()
