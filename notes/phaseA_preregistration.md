# Phase A pre-registration

Committed 2026-11-10, **before** any AlphaFold error data (A2) was collected
or any Phase A metric computed. Changes after this point are added below
under "Amendments", dated, with the reason; the original text is not edited.

## Question

Do AlphaFold's real errors concentrate in a small number of *predictable*
coupled directions, enough that a coupled prior beats the strongest simple
alternatives when fitting predicted structures to sparse data?

## Data (A2)

- **Proteins**: experimental structures released on or after 2023-01-01
  (after AlphaFold DB v4 was built, so no model used its own structure as a
  template), X-ray <= 2.5 Å or cryo-EM <= 3.0 Å, one chain, 60–500
  residues, with an AlphaFold DB model covering >= 90% of the chain.
- **Splits**: MMseqs2 clustering at 30% identity / 80% coverage; clusters
  assigned 60/20/20 train/val/test by hash (as `auf.splits`). Target
  >= 300 test proteins.
- **Secondary set**: ATLAS proteins (MD coupling available), split by
  `splits/atlas_splits.tsv`. AlphaFold may have trained on these; reported
  separately, never pooled with the primary set.
- **Exclusions** (fixed now): chains with > 10% missing residues in the
  modelled range; NMR structures; proteins whose AlphaFold model has
  mean pLDDT < 50 (no structure to correct).

## Error vectors

Cα only. AlphaFold model superposed onto the experimental structure on
residues with pLDDT >= 70 (**primary**) and on all residues (secondary).
e = x_exp − x_AF (3N), with the six rigid-body directions of the AF model
projected out.

## Direction sources (3N x 3N covariances, rigid modes removed)

1. PAE-derived, three constructions (roadmap A3: relative-displacement
   embedding, PAE rigid domains + pLDDT noise, low-rank truncation).
   **One** is chosen on the validation split by the A4 capture metric; only
   that one enters the gate. All three are reported on test as secondary.
2. ANM (Cα, 15 Å cutoff) of the AlphaFold model.
3. MD coupling: ATLAS secondary set only.
4. Random: isotropic Gaussian random subspace, same rank.

## Metrics

**Capture ratio** (primary k = 10; also 5, 20):
C_k(S) = (|P_k e|² / |e|²) / (k / (3N − 6)), where P_k projects onto the top-k
eigenvectors of S. C = 1 is random; C = 2 means the top-k directions carry
twice their share of the error. Per protein; mean over test proteins with a
95% bootstrap CI over proteins.

**Fusion**: prior mean = AlphaFold model; prior covariance = source, scaled
by one constant fitted on train (held-out likelihood of e). Data = Cα–Cα
distances measured on the *experimental* structure, pairs 5–25 Å apart,
0.5 per residue, noise sd 1.0 Å (**primary setting**; others secondary).
Outcome = improvement over the prior mean in Cα RMSD and in lDDT-Cα, versus
the experimental structure. Paired per protein.

## Controls that a coupled prior must beat

- **Oracle scalar**: per-residue isotropic prior with variance |e_i|²/3 set
  from the *true* error (an upper bound for any scalar prior).
- **PAE rigid domains + rigid fitting**: PAE graph clustering (Phenix
  convention, PAE cutoff 5 Å), each domain free to move rigidly (6 DOF,
  broad prior), per-residue noise from pLDDT; the current practical
  standard.
- **Random subspace** of the same rank as the chosen source.

## Gate (A5)

On the test split, primary settings:

- **GO** if the val-chosen PAE source (or ANM) has
  (i) mean C_10 with 95% CI lower bound > 2.0, **and**
  (ii) in fusion, a paired improvement over **both** the oracle scalar and
  PAE-rigid-domains, on **both** RMSD and lDDT, with 95% CI excluding zero,
  better in >= 60% of test proteins.
- **PARTIAL** if only MD coupling (ATLAS secondary set) meets (i) and (ii):
  coupling is real but not yet predictable from cheap inputs → roadmap B4 /
  B4b first.
- **NO-GO** otherwise → Phase C.

Whatever the outcome, all sources, controls and settings are reported, and
the Phase A write-up is produced.

## Amendments

(none)
