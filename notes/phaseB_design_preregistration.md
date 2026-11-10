# Phase B probe: choosing which distances to measure (pre-registration)

Committed 2026-11-10, before any selection result was computed. Changes after
this point go under "Amendments", dated, with the reason.

## Question

Given an AlphaFold model and its covariance, does choosing *which* residue
pairs to measure by expected error reduction remove more of AlphaFold's real
error per measurement than the choices people make today?

## Data

The A2 AlphaFold error set (`splits/afset_splits.tsv`). Development and
debugging on **val**; the result is reported on **test** (266 proteins).
Ground truth = the experimental structure. Measurements are synthetic:
Cα–Cα distances of the experimental structure plus Gaussian noise.

## Candidate pairs (chosen from the AlphaFold model only)

- **Generic**: residue pairs with AF-model Cα–Cα distance 5–25 Å, |i − j| ≥ 3.
- **Crosslinkable** (secondary): both residues Lys (or the N-terminal residue),
  AF-model Cα–Cα distance ≤ 25 Å, |i − j| ≥ 3.

No strategy ever sees the experimental structure.

## Selection strategies (m measurements per protein)

1. **random**: uniform over candidates (3 draws, averaged). Today's default.
2. **plddt_pairs**: greedy on (100 − pLDDT_i) + (100 − pLDDT_j), each residue
   used at most once. The "measure the low-confidence parts" heuristic.
3. **greedy_scalar**: greedy expected-variance reduction (below) under the
   pLDDT scalar covariance. Same algorithm, no coupling: isolates the value
   of coupling from the value of the algorithm.
4. **greedy_pae** (primary): greedy expected-variance reduction under the
   val-chosen PAE covariance (`pae_lowrank`, frozen A3 parameters).
5. **greedy_anm** (secondary): same under ANM with recalibrated spectrum.

Greedy expected-variance reduction (A-optimal, linearised at the AF model):
pick the candidate with the largest |P h|² / (hᵀ P h + σ²), where h is the
distance Jacobian row, then update P ← P − P h hᵀ P / (hᵀ P h + σ²); repeat.

## Fusion and outcome

All strategies are fused with the **same** prior (`pae_lowrank`, frozen scale)
so only the choice of pairs differs (secondary: each greedy strategy with its
own prior). Noise sd 1.0 Å (primary); 3 noise realisations per selection.
Budgets m = 5, 10, 20, 50; **primary m = 20**. Outcome: Cα RMSD removed
relative to the AF model (Å); secondary lDDT-Cα, and measurements needed to
match random at m = 50.

## Success criterion (decided now)

On test, generic candidates, m = 20: **greedy_pae** beats **random** and
**greedy_scalar** in RMSD removed, paired mean difference with 95% bootstrap
CI excluding zero, better in ≥ 60% of proteins. Reported for all proteins
and for AF error > 2 Å (subgroup, not part of the criterion).

If greedy_pae beats random but not greedy_scalar, the benefit comes from the
selection algorithm and pLDDT, not from coupling, and is reported that way.
