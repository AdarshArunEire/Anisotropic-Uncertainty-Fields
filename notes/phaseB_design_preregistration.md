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

## Amendments

### 2026-11-10, amendment 1 (before any selection result was computed)

**Added strategies** after reading the prior art
([literature_experiment_design.md](literature_experiment_design.md)):
- **heuristic_deer**: the Kazmier 2011 / DEERFold rule. Candidates restricted
  to residues in helices or strands (DSSP on the AF model) and exposed; pairs
  > 15 Å apart; chosen greedily to maximise sequence separation while
  covering each pair of secondary-structure elements at most once.
- **ensemble_fps**: FPS 2.0-style. 200 conformers sampled from the same prior
  covariance as greedy_pae; greedy forward choice of the pair that most
  reduces the ensemble-averaged expected RMSD.

**Success criterion extended:** greedy_pae must also beat heuristic_deer
(same paired test). Beating ensemble_fps is reported but not required, since
it shares greedy_pae's prior; the comparison measures closed form vs sampling
at their respective costs.

**Secondary:** all strategies with the feasibility filter (exposed, helix or
strand, > 15 Å). **Second benchmark (separate pre-registration before it
runs):** DEERFold's two-state targets, outcome = reaching the other state.

### 2026-11-12, amendment 2 (before any selection result was computed)

Decided in discussion on 2026-11-10/11 and after reading Hays et al. 2018,
Optimal Probes and Labelizer; no selection or fusion result of this probe
had been computed.

1. **Headline outcome: lDDT-Cα gain** over the AF model (what users read,
   same scale as pLDDT). Cα RMSD removed becomes secondary. The success
   criterion is otherwise unchanged: at m = 20, generic candidates, test
   split, greedy_pae beats random, greedy_scalar and heuristic_deer in lDDT-Cα
   gain (paired mean, 95% bootstrap CI excluding zero, better in >= 60% of
   proteins).
2. **Chain restraints in fusion** for every strategy: consecutive Cα–Cα
   distances fixed to their AF-model values (sd 0.05 Å) as pseudo-
   observations, so fusion stops damaging local geometry (Phase A).
   A **chain-only** run (m = 0) is reported as the reference.
3. **heuristic_deer, fixed definition (Cα only):** helix residue if
   d(i, i+3) < 6.0 Å and d(i, i+4) < 6.8 Å; strand residue if
   d(i, i−2, i+2 span) > 6.2 Å for both neighbours (d(i−1, i+1) > 6.2 Å);
   exposed if its Cα has fewer neighbours within 10 Å than the protein's
   median. Candidates: both residues helix/strand and exposed, AF distance
   15–60 Å (> 15 Å as in DEERFold; ≤ 60 Å DEER window), |i − j| ≥ 3. Greedy:
   each step takes the candidate maximising log sequence separation, skipping
   pairs whose two secondary-structure elements are already connected, until
   all element pairs are used, then by separation alone.
4. **ensemble_mrmr (Hays 2018 / FPS 2.0 MI mode):** 200 conformers drawn from
   N(AF model, scale × pae_lowrank); per candidate, distance variance v and
   correlation ρ with chosen pairs. Relevance ½ log(1 + v / σ²); redundancy =
   mean over chosen of −½ log(1 − ρ²); greedy on relevance − redundancy.
5. **Predicted-curve check (secondary):** for greedy_pae, the predicted
   expected squared error after m pairs, tr(P_post) (linearised, frozen
   scale), against the realised squared Cα error after fusion. Reported:
   Spearman correlation across proteins and the median ratio
   realised / predicted at m = 5, 10, 20, 50. A ratio far from 1 means the
   "how many pairs" curve would need calibration before use.
6. Fixed details: noise sd 1.0 Å; 3 noise realisations per selection; random:
   3 draws × 3 realisations; budgets nested (the m = 5 set is the first 5 of
   the m = 50 sequence for greedy methods); candidate pool for generic
   strategies capped at 4,000 pairs per protein by a seeded subsample if
   larger (compute); development and debugging on 20 validation proteins only.
