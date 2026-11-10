"""Stage 1 audit in physical units: Å, degrees, Å^3 and coverage.

Answers "is the scalar a bad guess?" without abstract distances:
- RMS motion along each principal axis vs the scalar's single value;
- angle between predicted and actual principal axis (random = 60 deg);
- held-out coverage of a nominal 95% region (calibration, as-is);
- volume needed for a region that truly covers 95% of held-out positions,
  ellipsoid vs sphere (efficiency once both are calibrated).

"Held-out" is replica j, predicted from replica i. Positions are Gaussian
draws from Sigma_j, so this isolates covariance shape and size from
between-replica mean shifts.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
from scipy import stats

from auf.audit import isotropize


def axis_angle(a, b):
    va = np.linalg.eigh(a)[1][..., -1]
    vb = np.linalg.eigh(b)[1][..., -1]
    return np.degrees(np.arccos(np.clip(np.abs(np.sum(va * vb, -1)), 0, 1)))


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--data", default="data/audit")
    p.add_argument("--draws", type=int, default=2000)
    p.add_argument("--out", default="results/04_physical_units.json")
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()
    rng = np.random.default_rng(args.seed)
    c95 = stats.chi2.ppf(0.95, 3)

    sig, ang_rep, ang_anm = [], [], []
    cov = {"ellipsoid": [], "sphere": []}
    vol = {"ellipsoid": [], "sphere": []}
    for f in sorted(Path(args.data).glob("*.npz")):
        with np.load(f) as z:
            s, reps, anm = z["sigma"], z["sigma_rep_mean"], z["anm"]
        sig.append(s)
        ang_anm.append(axis_angle(anm, s))
        for i, j in itertools.permutations(range(len(reps)), 2):
            ang_rep.append(axis_angle(reps[i], reps[j]))
        for i, j in ((0, 1), (1, 2), (2, 0)):
            chol = np.linalg.cholesky(reps[j])
            x = np.einsum("rab,nrb->nra", chol, rng.normal(size=(args.draws, s.shape[0], 3)))
            for name, pred in (("ellipsoid", reps[i]), ("sphere", isotropize(reps[i]))):
                m = np.einsum("nra,rab,nrb->nr", x, np.linalg.inv(pred), x)
                cov[name].append((m <= c95).mean(0))
                k = np.quantile(m, 0.95, axis=0)
                vol[name].append(4 / 3 * np.pi * k**1.5 * np.sqrt(np.linalg.det(pred)))

    s = np.concatenate(sig)
    sd = np.sqrt(np.linalg.eigvalsh(s))[:, ::-1]  # long, mid, short
    iso_sd = np.sqrt(np.trace(s, axis1=-2, axis2=-1) / 3)
    a_rep, a_anm = np.concatenate(ang_rep), np.concatenate(ang_anm)
    ce, cs = np.concatenate(cov["ellipsoid"]), np.concatenate(cov["sphere"])
    ve, vs = np.concatenate(vol["ellipsoid"]), np.concatenate(vol["sphere"])
    ratio = vs / ve

    res = {
        "n_proteins": len(sig), "n_residues": int(len(s)),
        "rms_long_mid_short_A": np.median(sd, 0).round(3).tolist(),
        "long_over_short_quartiles": np.percentile(sd[:, 0] / sd[:, 2], [25, 50, 75]).round(2).tolist(),
        "scalar_rms_A": float(np.median(iso_sd).round(3)),
        "scalar_understates_long_axis_pct": float(np.median(100 * (1 - iso_sd / sd[:, 0])).round(1)),
        "scalar_overstates_short_axis_x": float(np.median(iso_sd / sd[:, 2]).round(2)),
        "axis_angle_deg_median": {"random": 60.0, "anm": float(np.median(a_anm).round(1)),
                                  "other_replica": float(np.median(a_rep).round(1))},
        "coverage_of_nominal_95": {"ellipsoid": float(np.median(ce).round(3)),
                                   "sphere": float(np.median(cs).round(3))},
        "calibrated_95_volume_A3_median": {"ellipsoid": float(np.median(ve).round(1)),
                                           "sphere": float(np.median(vs).round(1))},
        "sphere_over_ellipsoid_volume_quartiles": np.percentile(ratio, [25, 50, 75]).round(2).tolist(),
        "share_ellipsoid_smaller": float((ratio > 1).mean().round(3)),
        "share_ellipsoid_1p5x_smaller": float((ratio > 1.5).mean().round(3)),
        "share_ellipsoid_2x_smaller": float((ratio > 2).mean().round(3)),
    }
    print(json.dumps(res, indent=2))
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
