"""Stage 1, question 1: how much information is in the anisotropy?

For each ATLAS entry, per-residue Cα covariances are computed after joint
superposition of all three replicas (see ``auf.data.atlas``). Three tests:

A. Magnitude of anisotropy vs sampling noise.
   signal = AIRM(Sigma, isotropize(Sigma)), residue by residue.
   null_W = AIRM(Sigma, Wishart resample of Sigma at n_eff): how far a
            finite-trajectory estimate sits from its own truth. (The
            original scaffold null.)

B. Is the anisotropy real?  (isotropic-truth null)
   AIRM(S, isotropize(S)) for S ~ W(I/n_eff, n_eff): the anisotropy a
   *sphere* shows when sampled at n_eff independent frames. Gives a
   per-residue p-value; reported with Benjamini-Hochberg at q = 0.05.

C. Is there directional structure beyond a per-residue scale?  (replica split)
   For every ordered pair of replicas (i, j): does Sigma_i predict Sigma_j
   better than isotropize(Sigma_i) does? Both keep residue i's trace, so a
   win can only come from direction.

   Chance is *not* 50%. Superposition bias is a deterministic function of the
   structure, so it reproduces across replicas and can make isotropic noise
   look directionally consistent. The control replays C on isotropic
   Gaussian fluctuations about the entry's own reference structure, at each
   residue's real RMSF, with many independent frames so the artifact
   dominates sampling noise. That is the ceiling the artifact alone reaches.

   Kill criterion: if C does not clear its control, the anisotropy is a
   per-residue scale plus noise and artifact, and Stage 2 is not justified.

Caveat: replicas start from the same structure, so C measures MD
reproducibility of direction, not agreement with experiment.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np

from auf.audit import (
    airm_distance,
    fractional_anisotropy,
    isotropize,
    wishart_isotropic_null_airm,
    wishart_null_airm,
)

# Fixed, length-stratified random draw from the 2024-11-18 ATLAS release
# (seed 0, 60-200 residues). Not chosen by flexibility or by outcome.
DEFAULT_ENTRIES = [
    "3czz_B", "1jyo_F", "1ng6_A", "2qev_A",
    "1dvo_A", "2ygs_A", "2qsk_A", "2fp1_A",
]


def bh_reject(p: np.ndarray, q: float = 0.05) -> np.ndarray:
    """Benjamini-Hochberg rejections at FDR q."""
    p = np.asarray(p)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    passed = np.nonzero(ranked <= q)[0]
    out = np.zeros(len(p), bool)
    if passed.size:
        out[order[: passed.max() + 1]] = True
    return out


def iso_null_pvalues(obs: np.ndarray, n_eff: np.ndarray, rng, n_draws: int) -> np.ndarray:
    """Per-residue p-value of the observed anisotropy under isotropic truth."""
    p = np.empty(len(obs))
    cache: dict[int, np.ndarray] = {}
    for k, (o, n) in enumerate(zip(obs, n_eff)):
        n = int(max(n, 4))
        if n not in cache:
            cache[n] = np.sort(wishart_isotropic_null_airm(n, n_draws=n_draws, rng=rng))
        null = cache[n]
        p[k] = (n_draws - np.searchsorted(null, o, side="left") + 1) / (n_draws + 1)
    return p


def replica_direction_wins(sigma_rep: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Test C: per (ordered replica pair, residue), does Sigma_i beat iso(Sigma_i)?"""
    wins, ratios = [], []
    for i, j in itertools.permutations(range(len(sigma_rep)), 2):
        si, sj = sigma_rep[i], sigma_rep[j]
        d_full = airm_distance(si, sj)
        d_iso = airm_distance(isotropize(si), sj)
        wins.append(d_full < d_iso)
        ratios.append(d_full / d_iso)
    return np.stack(wins), np.stack(ratios)


def isotropic_control(cov, rng, n_frames: int) -> np.ndarray:
    """Test C on isotropic truth: superposition artifact ceiling, per residue."""
    from auf.data.atlas import covariances_from_trajectories

    n_rep = len(cov.sigma_rep)
    scale = np.sqrt(np.trace(cov.sigma, axis1=-2, axis2=-1) / 3.0)
    noise = rng.normal(size=(n_rep, n_frames, cov.n_res, 3))
    frames = cov.reference + scale[:, None] * noise
    ctl = covariances_from_trajectories("control", frames)
    wins, _ = replica_direction_wins(ctl.sigma_rep)
    return wins.mean(axis=0)


def analyse(cov, rng, n_draws: int, control_frames: int) -> dict:
    sigma = cov.sigma
    n_eff = cov.n_eff
    signal = airm_distance(sigma, isotropize(sigma))

    # A: Wishart null about Sigma, at each residue's pooled n_eff.
    null_w = np.array([
        np.median(wishart_null_airm(s, n_eff=int(max(n, 4)), n_draws=64, rng=rng))
        for s, n in zip(sigma, n_eff)
    ])

    # B: isotropic-truth null, per replica (each replica is an independent
    # estimate at its own n_eff), then pooled.
    p_rep = []
    for s_r, n_r in zip(cov.sigma_rep, cov.n_eff_rep):
        obs_r = airm_distance(s_r, isotropize(s_r))
        p_rep.append(iso_null_pvalues(obs_r, n_r, rng, n_draws))
    p_rep = np.stack(p_rep)
    p_pool = iso_null_pvalues(signal, n_eff, rng, n_draws)

    # C: replica split, direction only, against the superposition control.
    wins, ratios = replica_direction_wins(cov.sigma_rep)
    control = isotropic_control(cov, rng, control_frames)

    return {
        "pdb_chain": cov.pdb_chain,
        "n_res": cov.n_res,
        "n_frames_per_replica": cov.n_frames,
        "n_eff_median": float(np.median(n_eff)),
        "fa_median": float(np.median(fractional_anisotropy(sigma))),
        "A_signal_median": float(np.median(signal)),
        "A_wishart_null_median": float(np.median(null_w)),
        "A_signal_over_null": float(np.median(signal) / np.median(null_w)),
        "A_frac_residues_signal_gt_null": float(np.mean(signal > null_w)),
        "B_frac_bh_pooled": float(bh_reject(p_pool).mean()),
        "B_frac_bh_per_replica": [float(bh_reject(p).mean()) for p in p_rep],
        "C_frac_aniso_beats_iso": float(wins.mean()),
        "C_median_dist_ratio_full_over_iso": float(np.median(ratios)),
        "C_control_isotropic": float(control.mean()),
        "rmsf_corr_vs_atlas": (
            None if cov.atlas_rmsf is None else
            float(np.corrcoef(cov.rmsf_rep.ravel(), cov.atlas_rmsf.ravel())[0, 1])
        ),
        "_per_residue": {
            "signal": signal, "null_w": null_w, "p_pool": p_pool,
            "c_win_frac": wins.mean(axis=0), "c_ratio": np.median(ratios, axis=0),
            "c_control": control,
            "fa": fractional_anisotropy(sigma), "n_eff": n_eff,
            "rmsf": np.sqrt(np.trace(sigma, axis1=-2, axis2=-1)),
        },
    }


def report(r: dict) -> None:
    print(f"== {r['pdb_chain']}  ({r['n_res']} residues, "
          f"{r['n_frames_per_replica']} frames x 3 replicas, "
          f"median n_eff {r['n_eff_median']:.0f} pooled)")
    if r["rmsf_corr_vs_atlas"] is not None:
        print(f"   RMSF vs ATLAS published:   r = {r['rmsf_corr_vs_atlas']:.4f}")
    print(f"   fractional anisotropy:     {r['fa_median']:.3f} (median)")
    print(f"   A  AIRM(Sigma, iso):       {r['A_signal_median']:.3f}  <- signal")
    print(f"      AIRM Wishart null:      {r['A_wishart_null_median']:.3f}  <- noise floor")
    print(f"      signal / noise:         {r['A_signal_over_null']:.2f}   "
          f"(residues above floor: {100 * r['A_frac_residues_signal_gt_null']:.0f}%)")
    per = ", ".join(f"{100 * f:.0f}%" for f in r["B_frac_bh_per_replica"])
    print(f"   B  anisotropic at FDR 5%:  {100 * r['B_frac_bh_pooled']:.0f}% pooled; "
          f"per replica {per}")
    print(f"   C  Sigma_i beats iso(Sigma_i) at predicting Sigma_j: "
          f"{100 * r['C_frac_aniso_beats_iso']:.0f}% of residue-pairs; "
          f"median distance ratio {r['C_median_dist_ratio_full_over_iso']:.2f}")
    verdict = "clears" if r["C_frac_aniso_beats_iso"] > r["C_control_isotropic"] + 0.10 else "DOES NOT clear"
    print(f"      isotropic-truth control: {100 * r['C_control_isotropic']:.0f}%  "
          f"-> {verdict} the superposition artifact (margin 10 pts)")


def synthetic(n: int, seed: int):
    """Scaffold-era synthetic input, kept for running without data."""
    from types import SimpleNamespace

    rng = np.random.default_rng(seed)
    a = rng.normal(size=(3, n, 3, 3))
    reps = np.swapaxes(a, -1, -2) @ a + 3.0 * np.eye(3)
    return SimpleNamespace(
        reference=rng.normal(scale=8.0, size=(n, 3)),
        pdb_chain="synthetic", sigma=reps.mean(0), sigma_rep=reps,
        n_eff=np.full(n, 300), n_eff_rep=np.full((3, n), 100), n_res=n,
        n_frames=1001, atlas_rmsf=None,
    )


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("entries", nargs="*", default=None,
                   help=f"ATLAS pdb_chain ids (default: {' '.join(DEFAULT_ENTRIES)})")
    p.add_argument("--data-dir", default="data/atlas")
    p.add_argument("--reference", choices=["mean", "start"], default="mean")
    p.add_argument("--synthetic", action="store_true",
                   help="run on synthetic covariances instead of ATLAS")
    p.add_argument("--null-draws", type=int, default=2000)
    p.add_argument("--control-frames", type=int, default=5000,
                   help="independent frames per replica in the isotropic control")
    p.add_argument("--out", default="results/01_information_content")
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()
    rng = np.random.default_rng(args.seed)

    if args.synthetic:
        covs = [synthetic(256, args.seed)]
    else:
        from auf.data.atlas import load_atlas_covariances

        covs = [load_atlas_covariances(e, args.data_dir, reference=args.reference)
                for e in (args.entries or DEFAULT_ENTRIES)]

    results = [analyse(c, rng, args.null_draws, args.control_frames) for c in covs]
    for r in results:
        report(r)

    if len(results) > 1:
        per = [r["_per_residue"] for r in results]
        cat = {k: np.concatenate([q[k] for q in per]) for k in per[0]}
        print(f"\n== all {len(results)} entries, {len(cat['signal'])} residues")
        print(f"   A  signal / noise:         "
              f"{np.median(cat['signal']) / np.median(cat['null_w']):.2f}")
        print(f"   B  anisotropic at FDR 5%:  {100 * bh_reject(cat['p_pool']).mean():.0f}%")
        print(f"   C  Sigma_i beats iso:      {100 * cat['c_win_frac'].mean():.0f}% "
              f"(median ratio {np.median(cat['c_ratio']):.2f}); "
              f"isotropic control {100 * cat['c_control'].mean():.0f}%")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    summary = [{k: v for k, v in r.items() if not k.startswith("_")} for r in results]
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    for r in results:
        np.savez_compressed(out / f"{r['pdb_chain']}_per_residue.npz", **r["_per_residue"])
    print(f"\nwrote {out}/")
    print("\nKill criterion: C not clearing its isotropic control means no")
    print("reproducible direction beyond per-residue scale; stop before Stage 2.")


if __name__ == "__main__":
    main()
