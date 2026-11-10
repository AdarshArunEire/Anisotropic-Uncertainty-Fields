# Downstream test: does uncertainty structure matter when fusing data?

Run 2026-11-10: `python experiments/audit/05_downstream_fusion.py`
(40 proteins from the Stage-1 random draw, <= 250 residues, 5,838 residues;
40 true frames x 2 replica pairs per protein). Literature context:
[literature_coupling.md](literature_coupling.md).

## Setup

- **Prediction**: mean structure of MD replica i. Its distance to the truth
  is 2.9 Å RMSE on average (between-replica conformational differences).
- **Truth**: individual frames of a different replica j.
- **Data**: noisy distances between random residue pairs 5–25 Å apart,
  0.1 or 0.5 per residue, noise 0.5 or 2 Å.
- **Update**: Gaussian conditioning on the distances (Gauss-Newton /
  iterated Kalman, `auf.audit.fusion`), with each prior below. All priors
  have the same per-residue magnitude except `joint_anm`.

| prior | what it knows |
|---|---|
| scalar | per-residue size only (pLDDT-style) |
| ellipsoid | + per-residue direction (what DynaProt / BackFlip-2 output per residue) |
| ellipsoid_shrunk | ellipsoid shrunk to alpha = 0.74 (Stage-1 audit) |
| anm_ellipsoid | per-residue shape from ANM, no learning |
| joint_anm | full 3N x 3N ANM: residues move together, one global scale |
| joint_anm_matched | ANM coupling, rescaled to MD per-residue sizes |
| joint_md | full 3N x 3N from replica i (near-oracle: needs a simulation of the same protein) |

## Results: share of the prediction's error removed by the data

Mean over proteins, [95% bootstrap CI].

| prior | 0.5/res, 0.5 Å | 0.5/res, 2 Å | 0.1/res, 0.5 Å | 0.1/res, 2 Å |
|---|---|---|---|---|
| scalar | 3.9% [3.0, 4.7] | 2.3% | 0.9% | 0.5% |
| ellipsoid | 4.3% [3.1, 5.4] | 2.8% | 0.8% | 0.8% |
| ellipsoid_shrunk | 4.7% [3.7, 5.8] | 3.0% | 1.0% | 0.8% |
| anm_ellipsoid | 1.6% | 0.6% | 0.5% | 0.3% |
| joint_anm | 5.1% [2.4, 7.6] | 1.9% | −0.6% | −0.2% |
| joint_anm_matched | 6.0% [3.5, 8.5] | 2.6% | 0.2% | −0.1% |
| **joint_md** | **25.0% [21.8, 28.3]** | **16.2% [13.5, 19.3]** | **6.6%** | **3.7%** |

Paired, densest/most precise setting:
- ellipsoid (shrunk) vs scalar: +0.03 Å, better in 80% of proteins. **Ignorable.**
- joint MD vs joint ANM at MD magnitudes: **+0.68 Å [0.47, 0.92], better in 98% of proteins.**

Calibration (90% posterior regions, all settings): every prior is
overconfident, covering the truth 57–70% of the time (calib scale 2.4–6.3,
1 = honest). The priors describe fluctuation about a structure, not the error
of the structure itself.

## Reading

1. **Per-residue anisotropy is ignorable downstream** (~0.02–0.03 Å).
2. **Inter-residue coupling is the effect**: correct joint covariance removes
   ~6x more error than per-residue priors, in 98% of proteins.
3. **ANM coupling does not capture it**, even with correct per-residue sizes.
   This conflicts with older ENM successes (Zheng & Brooks 2006, DEN), which
   used a few low modes toward a known target state; worth reconciling.
4. **All priors are badly miscalibrated** for this use.

## Caveats

- joint_md is near-oracle (another simulation of the same protein). The open
  question is how much of its benefit a *predicted* coupling recovers.
- Simulated data on simulated truth; real experimental restraints next.
- DynaProt and BackFlip-2 already predict pairwise coupling; their coupling
  outputs have not been tested here yet.

## Next

1. Put DynaProt / BackFlip-2 predicted coupling into this benchmark.
2. Reconcile with ENM literature: low-mode subspace priors, other cutoffs.
3. PAE -> joint covariance, calibrated; test in the same benchmark.
4. Real experimental restraints (crosslinks, NMR, DEER) on proteins with
   known structures.
