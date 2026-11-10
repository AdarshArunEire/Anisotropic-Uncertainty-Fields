# Phase A: Are AlphaFold's errors low-dimensional and predictable?

Run 2026-11-10 on `phaseA/af-error-set`. Criteria: [phaseA_preregistration.md](phaseA_preregistration.md).
Full tables: [results/phaseA/a4_report.md](../results/phaseA/a4_report.md) (pre-registered)
and [results/phaseA/a4b_exploratory.md](../results/phaseA/a4b_exploratory.md) (post hoc).

## Verdict

**Pre-registered gate: NO-GO on the primary set.** Criterion (i) passes by a
wide margin; criterion (ii) fails: no prior beats the oracle scalar in the
primary fusion setting. The ATLAS/MD secondary set, which decides between
NO-GO and PARTIAL, has not been run yet.

In plain terms: AlphaFold's errors *do* concentrate in directions that PAE
and a simple elastic network predict, but for the typical AlphaFold model
(median Cα error 1.0 Å, median pLDDT 95) sparse distance data cannot improve
it, coupled prior or not. The effect is ignorable there (≤ 0.07 Å).
Where AlphaFold is substantially wrong (> 2 Å, about one protein in five),
coupled priors remove roughly twice as much error as any scalar prior. But
the train-fitted PAE rigid-domain prior already gets that benefit, and every
prior damages local geometry (lDDT).

## Data (A2)

1,282 proteins with an experimental structure released after their AlphaFold
DB model was created (1,248 X-ray, 34 cryo-EM; median 261 residues; median
mean pLDDT 95). MMseqs2 splits (30% id): train 743 / val 273 / test 266. The
test split is below the pre-registered ≥ 300: the RCSB candidate list
(10,661 entities) ran out first. Biggest losses: AF model newer than the
entry (5,129), duplicate UniProt (2,439), no SIFTS mapping or no AFDB model
(1,029). Manifest: [results/phaseA/afset_manifest.tsv](../results/phaseA/afset_manifest.tsv).

## Recipes (A3, train only, frozen before A4)

Gaussian log-likelihood per dimension of the AF error on train: PAE low-rank
(rank 80, weight 0.3) −0.71; PAE rigid domains (σ_t = 2 Å, σ_r = 0.1 rad)
−0.78; ANM −0.85; pLDDT scalar −0.94. The raw relative-displacement
embedding (a) is near-singular on its own (−3.8). Rank 80 is at the edge of
the pre-fixed grid, which was not widened. Validation chose PAE low-rank by
capture (C_10 16.7 vs 16.2 domains vs 16.0 relative-displacement), committed
before any test metric.

## Results on test (266 proteins)

**Capture (gate i).** Share of the error in the top 10 directions, relative
to the pLDDT scalar (geometric mean [95% CI]):

| source | C_10 | ratio to pLDDT scalar |
|---|---|---|
| ANM (structure only) | 19.8 | **3.87 [3.24, 4.62]** |
| PAE low-rank (val-chosen) | 16.0 | **3.12 [2.70, 3.63]** |
| PAE rigid domains (fitted) | 15.3 | 2.99 [2.52, 3.58] |
| oracle scalar (knows each residue's true error) | 15.2 | 2.96 [2.62, 3.36] |
| pLDDT scalar | 5.1 | 1 |
| random | 0.9 | 0.18 |

Gate (i) needs a CI lower bound > 1.5: both pass. Notably, ANM's top
directions hold more of the error than an oracle that knows which residues
are wrong. *Direction* is informative, not only *which residues*.

**Fusion (gate ii).** Primary setting: 0.5 distances per residue, 1 Å noise.
Gain over the AF model, mean [95% CI]:

| prior | Cα RMSD removed (Å) | lDDT-Cα change (points) |
|---|---|---|
| oracle scalar | +0.087 [+0.069, +0.107] | −0.85 |
| PAE low-rank | +0.059 [+0.013, +0.121] | −4.15 |
| PAE rigid domains (fitted) | +0.064 [+0.004, +0.137] | −4.44 |
| ANM | +0.054 [−0.016, +0.135] | −3.71 |
| pLDDT scalar | −0.010 [−0.027, +0.009] | −4.49 |
| PAE rigid-domain control (broad, pre-registered) | −0.042 [−0.133, +0.066] | −6.84 |

PAE low-rank vs oracle scalar: RMSD −0.027 [−0.060, +0.020] Å (better in
21%), lDDT −3.3 points (better in 5%). **Fail.** Vs the pre-registered
broad rigid-domain control: +0.10 Å, +2.7 points, better in 83% / 94%
(pass), but against the *fitted* rigid-domain recipe the difference is nil
(−0.005 [−0.047, +0.030] Å).

A fusion sanity check (near-noiseless, 3 distances per residue) moves every
prior strongly toward the experiment (e.g. 1.04 → 0.35 Å), so the negative
result is not a bug. It comes from the regime: AF error ≈ data noise.

## Post hoc (A4b, exploratory, does not change the gate)

Same frozen parameters, noise 0.25 / 0.5 / 1 Å, split by the AF model's own
error. Test, noise 0.5 Å, Cα RMSD removed (mean, median):

| AF error group | n | oracle scalar | pLDDT scalar | PAE low-rank | PAE domains | ANM recal |
|---|---|---|---|---|---|---|
| ≤ 1 Å (median 0.63) | 130 | +0.03 (0.02) | −0.02 (−0.02) | −0.01 (−0.01) | −0.01 (−0.03) | +0.01 (0.01) |
| 1–2 Å (median 1.36) | 78 | +0.10 (0.09) | +0.07 (0.05) | +0.12 (0.10) | +0.14 (0.11) | +0.13 (0.12) |
| > 2 Å (median 3.81) | 58 | +0.26 (0.15) | +0.19 (0.12) | **+0.46 (0.30)** | **+0.59 (0.22)** | **+0.47 (0.31)** |

Validation shows the same pattern (> 2 Å, 59 proteins: PAE low-rank +0.86,
domains +0.73, oracle +0.48, pLDDT +0.40 Å). In the large-error group,
coupled priors remove about 12–15% of the error vs 5–7% for scalar priors,
and beat even the oracle scalar in 71–83% of proteins. lDDT still drops for
every prior in every setting: Cα-only Gaussian fusion moves domains the
right way but distorts local distances.

## What this means

1. **Real, citable finding:** AlphaFold's errors are low-dimensional, and
   PAE and ANM predict the directions (3–4x the pLDDT scalar, at or above an
   oracle scalar). This is the Phase A write-up's headline.
2. **Ignorable for most models:** for well-predicted proteins (≈ 80% of this
   set), no prior lets sparse 1 Å-noise distances improve the model.
3. **Adjacent positive (post hoc, to be confirmed prospectively):** when
   AlphaFold is wrong by > 2 Å, coupled priors roughly double the error
   removed. The PAE rigid-domain prior, *fitted* on train, already captures
   this. So the bar for anything new is "beat fitted PAE rigid domains on
   large-error cases", not "beat a scalar".
4. **Method gap:** every prior lowers lDDT. Adding local geometry (Cα–Cα
   3.8 Å bonds, a stiffer short-range prior) is needed before any fusion
   result is usable.

## Next steps (proposed, need sign-off)

- **A5 secondary set:** AF models of ATLAS test proteins with MD coupling,
  to settle NO-GO vs PARTIAL as pre-registered.
- **Pre-register a large-error follow-up** (roadmap B3 territory): proteins
  with AF error > 2 Å or known conformational changes; control = fitted PAE
  rigid domains; outcome = RMSD and lDDT with local geometry restraints.
- Otherwise Phase C as planned (C2 audit paper, C1 real-motion data).
