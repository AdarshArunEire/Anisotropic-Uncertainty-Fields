# Stage 1, Q1: is there directional signal in ATLAS? (first pass)

Run: `python experiments/audit/01_information_content.py` (2026-11-10).
Data: 8 ATLAS entries, a seeded random draw (seed 0) from the 996 entries of
60-200 residues in the 2024-11-18 release. Not selected on flexibility or
outcome. 1001 frames x 3 replicas each (the `analysis` archive), Cα only,
1004 residues total.

## Loader check

`sqrt(tr Σ)` per replica reproduces ATLAS's published RMSF at r >= 0.998 on
7/8 entries (max abs diff ~0.02 Å on 2erl_A). The exception is 1jyo_F
(r = 0.875), a largely disordered entry (ATLAS avg RMSF 6.35 Å), where the
joint three-replica superposition used here differs materially from ATLAS's
per-replica fit.

## Results

| entry  | n_res | FA (med) | A: signal/noise | B: aniso @FDR5% pooled | C: Σ beats iso(Σ) | C control (iso truth) |
|--------|------:|---------:|----------------:|-----------------------:|------------------:|----------------------:|
| 3czz_B |   101 | 0.57 | 2.18 | 78% | 89% | 19% |
| 1jyo_F |   105 | 0.73 | 2.55 | 84% | 49% | 24% |
| 1ng6_A |   148 | 0.81 | 4.60 | 99% | 91% | 22% |
| 2qev_A |   145 | 0.65 | 2.35 | 92% | 79% | 26% |
| 1dvo_A |   152 | 0.72 | 5.11 | 100% | 89% | 22% |
| 2ygs_A |    92 | 0.44 | 2.62 | 88% | 95% | 22% |
| 2qsk_A |    95 | 0.61 | 6.27 | 100% | 97% | 27% |
| 2fp1_A |   166 | 0.56 | 1.87 | 67% | 65% | 46% |
| **all**| 1004 |      | 2.87 | 88% | 81% | 27% |

Control numbers vary by a few points between runs (Monte Carlo).

**Kill criterion not triggered on this sample.** Each residue's covariance
from one replica predicts another replica's covariance better than its
trace-matched isotropic surrogate does, for 81% of residue-pairs, against
27% for isotropic truth passed through the same superposition. The
anisotropy is not a global or per-residue scale factor: its direction
reproduces across independent 100 ns trajectories.

## Things learned along the way

1. **Chance for test C is not 50%.** When the truth is isotropic, a noisy
   Σ loses to its own isotropic surrogate, and C sits around 15-25%.
2. **Superposition bias is directional.** Fitting removes rigid motion but
   also absorbs part of each residue's fluctuation into the fitted rotation
   and leaks other residues' into it. It distorts *eigenvalues*, i.e. shape
   (see `test_superposition_bias_is_directional_on_small_proteins`: up to
   2x on the smallest eigenvalue at 30 residues). Because it is a
   deterministic function of the structure, it reproduces across replicas
   and can fake directional agreement. On 2fp1_A it alone reaches 46%, so
   that entry's 65% is a thin margin. Any published directional metric
   computed on globally superposed MD inherits this. Per-residue local
   frames (BackFlip-2's `R_i Σ_loc R_iᵀ`, and the Stage 2 se(3) design)
   sidestep global superposition, which is an argument for them
   independent of expressiveness.
3. **Per-replica test B is underpowered**, not negative. With the
   conservative n_eff (max autocorrelation time over x/y/z, Sokal window),
   per-replica n_eff is often 10-40, where even a sphere shows AIRM ~0.5-1.
   Test C does not depend on n_eff and is the more reliable read.
4. **Pooled Σ includes between-replica conformational shifts.** Test A and
   pooled B therefore credit some anisotropy to replicas sitting in
   different basins. Test C uses within-replica Σ only.

## Caveats

- Eight proteins. Scale to a few hundred before quoting numbers.
- Replicas share a starting structure. C measures MD reproducibility of
  direction, not agreement with experiment or with a predictor.
- 1jyo_F (disordered): C 49% vs control 24% clears formally, but its median
  distance ratio is 1.01, so direction there is barely better than scale.
- This establishes that the *ground truth* carries direction. It says
  nothing yet about whether DynaProt / BackFlip-2 *predictions* capture it
  (Q2) or are calibrated (Q3).

## Next

1. Scale Q1 to the ATLAS test split used by DynaProt / BackFlip-2.
2. Obtain their predicted covariances on those entries and run Q2 (the
   isotropic baseline) and Q3 (χ²₃ calibration of held-out frame
   deviations, with `scale_factor` separating magnitude from shape).
3. Repeat C with a local-frame (per-residue N-Cα-C) covariance to see how
   much of the directional signal survives without global superposition.
