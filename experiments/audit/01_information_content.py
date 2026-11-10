"""Stage 1, question 1: how much information is in the anisotropy?

Replaces each ground-truth covariance with its trace-matched isotropic
surrogate and measures the distance lost, against the Wishart sampling-noise
null for the trajectory length the covariance was estimated from.

Reads nothing real yet: the loader is a stub. Run it to see the shape of the
comparison on synthetic covariances.
"""
from __future__ import annotations

import argparse

import numpy as np

from auf.audit import airm_distance, fractional_anisotropy, isotropize, wishart_null_airm


def load_covariances(path: str | None, n: int, seed: int) -> np.ndarray:
    """Placeholder. Replace with ATLAS-derived per-residue covariances.

    Target shape is (n_residues, 3, 3), computed after superposing trajectory
    frames onto a reference conformation.
    """
    if path is not None:
        raise NotImplementedError("ATLAS loader not written yet")
    rng = np.random.default_rng(seed)
    a = rng.normal(size=(n, 3, 3))
    return np.swapaxes(a, -1, -2) @ a + 3.0 * np.eye(3)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--path", default=None, help="ATLAS covariances (not yet wired)")
    p.add_argument("--n", type=int, default=512, help="synthetic residue count")
    p.add_argument("--n-eff", type=int, default=100,
                   help="effectively independent trajectory frames")
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    sigma = load_covariances(args.path, args.n, args.seed)
    signal = airm_distance(sigma, isotropize(sigma))

    rng = np.random.default_rng(args.seed + 1)
    null = np.concatenate([
        wishart_null_airm(s, n_eff=args.n_eff, n_draws=32, rng=rng) for s in sigma[:64]
    ])

    print(f"residues:                  {len(sigma)}")
    print(f"fractional anisotropy:     {np.median(fractional_anisotropy(sigma)):.4f}")
    print(f"AIRM(Sigma, isotropic):    {np.median(signal):.4f}  <- signal")
    print(f"AIRM Wishart null (n={args.n_eff}): {np.median(null):.4f}  <- noise floor")
    ratio = np.median(signal) / np.median(null)
    print(f"signal / noise:            {ratio:.2f}")
    print()
    print("Ratio near 1 means the apparent anisotropy is consistent with")
    print("finite-sampling noise, and Stage 2 is not justified.")


if __name__ == "__main__":
    main()
