"""Build the Stage-1 audit dataset from ATLAS.

For each entry in a seeded random draw (not selected on flexibility or on
any outcome), downloads the 1000-frame archive, computes everything the
audit needs, writes ``<out>/<pdb_chain>.npz`` (a few hundred kB), and deletes
the raw trajectories. Resumable: finished entries are skipped.

Per entry:
  sigma_rep_{mean,start,core} (3, n, 3, 3)  replica covariances under three
      superposition choices: iterative mean, ATLAS starting structure, and
      mean-fit on the rigid half of the protein only (pooled RMSF <= median).
  sigma (n, 3, 3)        pooled covariance, mean fit.
  n_eff_rep (3, n)       effective frames per replica.
  reference (n, 3)       mean-fit reference (the frame everything is in).
  anm (n, 3, 3)          ANM covariance of the starting structure, mean frame.
  ss (n,)                DSSP (H/E/C) of the starting structure.
  atlas_rmsf (3, n)      ATLAS published RMSF.
  crystal_u (n, 3, 3)    deposited ANISOU of Cα, mean frame, NaN if absent.
  crystal_rmsd           Cα RMSD between crystal and MD reference.
"""
from __future__ import annotations

import argparse
import csv
import shutil
import sys
import traceback
from multiprocessing import Pool
from pathlib import Path

import numpy as np


def sample_entries(info_tsv: Path, n: int, lo: int, hi: int, seed: int) -> list[str]:
    with open(info_tsv) as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if lo <= int(r["length"]) <= hi]
    rng = np.random.default_rng(seed)
    return [rows[i]["PDB"] for i in rng.choice(len(rows), size=min(n, len(rows)), replace=False)]


def build_one(args) -> str:
    pdb_chain, raw_dir, out_dir, pdb_dir, keep_raw = args
    out = Path(out_dir) / f"{pdb_chain}.npz"
    if out.exists():
        return f"{pdb_chain}: cached"
    try:
        import mdtraj as md

        from auf.data.anm import anm_covariance
        from auf.data.atlas import (
            covariances_from_trajectories, fetch_atlas, kabsch_rotation,
            read_atlas_rmsf, read_ca_trajectories, superpose,
        )
        from auf.data.crystal import crystal_for_atlas

        entry = fetch_atlas(pdb_chain, raw_dir)
        ca = read_ca_trajectories(entry, pdb_chain)
        start = md.load(str(entry / f"{pdb_chain}.pdb"))
        start_ca = start.xyz[0, start.top.select("name CA")].astype(float) * 10.0
        ss = md.compute_dssp(start, simplified=True)[0]
        ss_ca = np.array([ss[a.residue.index] for a in start.top.atoms if a.name == "CA"])

        mean = covariances_from_trajectories(pdb_chain, ca)
        first = covariances_from_trajectories(pdb_chain, ca, reference=start_ca)
        rmsf = np.sqrt(np.trace(mean.sigma, axis1=-2, axis2=-1))
        core = covariances_from_trajectories(pdb_chain, ca, fit_mask=rmsf <= np.median(rmsf))

        ref = mean.reference
        start_in_ref, _ = superpose(start_ca[None], reference=ref)
        anm = anm_covariance(start_in_ref[0])

        crystal_u = np.full_like(mean.sigma, np.nan)
        crystal_rmsd = np.nan
        cr = crystal_for_atlas(pdb_chain, entry, pdb_dir)
        if cr is not None and len(cr[0]) == len(ref):
            xyz, u = cr
            m = np.isfinite(xyz[:, 0])
            x = xyz[m] - xyz[m].mean(0)
            y = ref[m] - ref[m].mean(0)
            rot = kabsch_rotation(x, y)
            crystal_rmsd = float(np.sqrt(np.mean(np.sum((x @ rot - y) ** 2, axis=1))))
            crystal_u[m] = rot.T @ u[m] @ rot

        np.savez_compressed(
            out, pdb_chain=pdb_chain,
            sigma_rep_mean=mean.sigma_rep, sigma_rep_start=first.sigma_rep,
            sigma_rep_core=core.sigma_rep, sigma=mean.sigma,
            n_eff_rep=mean.n_eff_rep, reference=ref, anm=anm, ss=ss_ca,
            atlas_rmsf=read_atlas_rmsf(entry, pdb_chain),
            crystal_u=crystal_u, crystal_rmsd=crystal_rmsd,
        )
        if not keep_raw:
            shutil.rmtree(entry, ignore_errors=True)
        return f"{pdb_chain}: ok ({len(ref)} res, crystal {'yes' if np.isfinite(crystal_rmsd) else 'no'})"
    except Exception as exc:  # keep going; failures are listed at the end
        return f"{pdb_chain}: FAILED {type(exc).__name__}: {exc}\n{traceback.format_exc(limit=2)}"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--info", default="data/2024_11_18_ATLAS_info.tsv",
                   help="ATLAS release table (downloaded if missing)")
    p.add_argument("--n", type=int, default=150)
    p.add_argument("--min-len", type=int, default=60)
    p.add_argument("--max-len", type=int, default=300)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--entries", nargs="*", help="explicit entries instead of a draw")
    p.add_argument("--raw-dir", default="data/atlas")
    p.add_argument("--out", default="data/audit")
    p.add_argument("--pdb-dir", default="data/pdb")
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--keep-raw", action="store_true")
    args = p.parse_args()

    info = Path(args.info)
    if not info.exists():
        import urllib.request

        info.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(
            "https://www.dsimb.inserm.fr/ATLAS/data/download/distributions/"
            "2024_11_18_ATLAS_info.tsv", info)
    entries = args.entries or sample_entries(info, args.n, args.min_len, args.max_len, args.seed)
    Path(args.out).mkdir(parents=True, exist_ok=True)
    (Path(args.out) / "entries.txt").write_text("\n".join(entries) + "\n")

    jobs = [(e, args.raw_dir, args.out, args.pdb_dir, args.keep_raw) for e in entries]
    failed = []
    with Pool(args.workers) as pool:
        for k, msg in enumerate(pool.imap_unordered(build_one, jobs), 1):
            print(f"[{k}/{len(jobs)}] {msg}", flush=True)
            if "FAILED" in msg:
                failed.append(msg.split(":")[0])
    print(f"\ndone; {len(failed)} failed: {' '.join(failed)}")
    sys.exit(0)


if __name__ == "__main__":
    main()
