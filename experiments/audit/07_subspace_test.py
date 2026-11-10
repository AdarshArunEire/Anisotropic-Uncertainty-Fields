"""Phase A1: do simulation errors concentrate in predicted coupled directions?

For each ATLAS protein (Stage-1 random draw) and replica pair (i -> j),
errors are measured from replica i's mean structure:
  - mean shift:  mean_j - mean_i                (a conformational error)
  - frames:      x_t(j) - mean_i, 50 frames      (what the fusion test fits)
and scored with the capture ratio C_k (``auf.audit.subspace``) against
direction sources built from replica i or from structure alone:

  md_joint    full 3N x 3N covariance of replica i      (coupling)
  md_blocks   its per-residue 3x3 blocks only            (no coupling)
  md_scalar   per-residue variance only                  (no direction, no coupling)
  anm         full ANM of the starting structure         (cheap coupling)
  random      Wishart random covariance                  (C = 1 by construction)

md_blocks and md_scalar are the controls: if coupling carries the error,
md_joint must clearly beat them. Rigid-body directions are removed.
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

SOURCES = ["md_joint", "md_blocks", "md_scalar", "anm", "random"]
KS = (5, 10, 20)
PAIRS = [(0, 1), (1, 2), (2, 0)]


def run_one(args):
    pdb_chain, raw_dir, n_frames, seed, keep_raw = args
    import mdtraj as md

    from auf.audit.fusion import block_diag_3
    from auf.audit.subspace import capture_ratio
    from auf.data.anm import anm_covariance_full
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
    idx = np.arange(n)

    res = {s: {"mean": {k: [] for k in KS}, "frames": {k: [] for k in KS}} for s in SOURCES}
    for i, j in PAIRS:
        mu = aligned[i].mean(0)
        dev = (aligned[i] - mu).reshape(t, 3 * n)
        full = dev.T @ dev / (t - 1)
        blocks = full.reshape(n, 3, n, 3)[idx, :, idx, :]
        var = np.trace(blocks, axis1=-2, axis2=-1) / 3
        start_in, _ = superpose(start_ca[None], reference=mu)
        g = rng.normal(size=(3 * n, 3 * n + 10))
        srcs = {
            "md_joint": full,
            "md_blocks": block_diag_3(blocks),
            "md_scalar": np.diag(np.repeat(var, 3)),
            "anm": anm_covariance_full(start_in[0]),
            "random": g @ g.T,
        }
        e_mean = (aligned[j].mean(0) - mu).ravel()
        sel = np.linspace(0, t - 1, n_frames).astype(int)
        e_frames = (aligned[j][sel] - mu).reshape(n_frames, 3 * n)
        for name, cov in srcs.items():
            cm = capture_ratio(cov, e_mean, KS, rigid_ref=mu)
            cf = capture_ratio(cov, e_frames, KS, rigid_ref=mu)
            for k in KS:
                res[name]["mean"][k].append(float(cm[k]))
                res[name]["frames"][k].append(float(np.mean(cf[k])))
    if not keep_raw:
        shutil.rmtree(entry, ignore_errors=True)
    return {"pdb_chain": pdb_chain, "n_res": n,
            "res": {s: {e: {k: float(np.mean(v)) for k, v in d.items()}
                        for e, d in r.items()} for s, r in res.items()}}


def boot(x, rng, n_boot=4000):
    x = np.asarray(x)
    b = x[rng.integers(0, len(x), (n_boot, len(x)))].mean(1)
    return float(x.mean()), float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--entries", default="data/audit/entries.txt")
    p.add_argument("--max-len", type=int, default=300)
    p.add_argument("--frames", type=int, default=50)
    p.add_argument("--raw-dir", default="data/atlas")
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--keep-raw", action="store_true")
    p.add_argument("--out", default="results/07_subspace_test")
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    from auf.splits import load_splits

    splits = load_splits()
    entries = Path(args.entries).read_text().split()
    jobs = [(e, args.raw_dir, args.frames, args.seed + k, args.keep_raw)
            for k, e in enumerate(entries)]
    per = []
    with Pool(args.workers) as pool:
        for k, r in enumerate(pool.imap_unordered(run_one, jobs), 1):
            if "error" in r:
                print(f"[{k}/{len(jobs)}] {r['pdb_chain']}: ERROR {r['error']}", flush=True)
                continue
            if r["n_res"] > args.max_len:
                continue
            r["split"] = splits.get(r["pdb_chain"], "?")
            per.append(r)
            print(f"[{k}/{len(jobs)}] {r['pdb_chain']}: ok", flush=True)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "per_protein.json").write_text(json.dumps(per))

    rng = np.random.default_rng(args.seed)
    lines = [f"# A1 subspace test: {len(per)} ATLAS proteins, replica pairs {PAIRS}",
             "Capture ratio C_k = share of error in top-k directions / random share "
             "(1 = random). Mean over proteins [95% bootstrap CI].", ""]
    for subset in ("all", "test"):
        rows = per if subset == "all" else [r for r in per if r["split"] == "test"]
        for err in ("mean", "frames"):
            label = "replica mean shift" if err == "mean" else "individual frames"
            lines.append(f"## {label}, {subset} proteins (n={len(rows)})")
            lines.append(f"{'source':<11}" + "".join(f"{'C_' + str(k):>24}" for k in KS))
            for s in SOURCES:
                cells = []
                for k in KS:
                    m, lo, hi = boot([r["res"][s][err][str(k)] if str(k) in r["res"][s][err]
                                      else r["res"][s][err][k] for r in rows], rng)
                    cells.append(f"{m:>8.2f} [{lo:5.2f},{hi:5.2f}]")
                lines.append(f"{s:<11}" + "".join(f"{c:>24}" for c in cells))
            lines.append("")
    report = "\n".join(lines)
    print(report)
    (out / "report.md").write_text(report + "\n")


if __name__ == "__main__":
    main()
