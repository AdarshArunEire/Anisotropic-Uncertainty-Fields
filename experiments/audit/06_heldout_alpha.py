"""Phase 0.3: shrinkage alpha fitted on train proteins, scored on test.

Replaces the in-sample alpha of the Stage-1 audit (fitted on the same
replica pairs it was scored on). One alpha for all proteins, chosen to
maximise mean held-out Gaussian log-likelihood across train proteins'
replica pairs; validation picks nothing here (single parameter), and the
test number is reported once.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np

from auf.audit import isotropize
from auf.audit.direction import gaussian_loglik, shape_power
from auf.splits import load_splits

ALPHAS = np.linspace(0.0, 1.5, 61)


def gains(sig_rep: np.ndarray, alpha: float) -> float:
    out = []
    for i, j in itertools.permutations(range(len(sig_rep)), 2):
        si, sj = sig_rep[i], sig_rep[j]
        tr = np.trace(si, axis1=-2, axis2=-1)[:, None, None]
        out.append(np.mean(gaussian_loglik(shape_power(si, alpha) * tr, sj)
                           - gaussian_loglik(isotropize(si), sj)))
    return float(np.mean(out))


def main() -> None:
    splits = load_splits()
    data = {"train": [], "val": [], "test": []}
    for f in sorted(Path("data/audit").glob("*.npz")):
        with np.load(f) as z:
            data[splits[f.stem]].append(z["sigma_rep_mean"])

    curve = np.array([np.mean([gains(s, a) for s in data["train"]]) for a in ALPHAS])
    a_train = float(ALPHAS[np.argmax(curve)])
    rng = np.random.default_rng(0)
    res = {"n": {k: len(v) for k, v in data.items()}, "alpha_train": a_train}
    for split in ("val", "test"):
        per = {
            "raw": [gains(s, 1.0) for s in data[split]],
            "alpha_train": [gains(s, a_train) for s in data[split]],
            "alpha_0.74_in_sample": [gains(s, 0.74) for s in data[split]],
        }
        res[split] = {}
        for k, v in per.items():
            v = np.asarray(v)
            b = v[rng.integers(0, len(v), (4000, len(v)))].mean(1)
            res[split][k] = [round(float(v.mean()), 4),
                             round(float(np.percentile(b, 2.5)), 4),
                             round(float(np.percentile(b, 97.5)), 4)]
    print(json.dumps(res, indent=1))
    Path("results").mkdir(exist_ok=True)
    Path("results/06_heldout_alpha.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
