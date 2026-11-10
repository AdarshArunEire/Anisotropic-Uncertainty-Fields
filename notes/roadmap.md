# Roadmap

Agreed 2026-11-10. Goal (see CLAUDE.md): a result that changes practice,
found honestly. Every phase ends in a written result, so nothing is wasted
if a later gate fails.

Notation: N residues; a "direction source" is any 3N x 3N covariance whose
top eigenvectors claim to say *which way* a structure is wrong or moves.

---

## Phase 0: Ground rules and housekeeping (before any new result)

0.1 **Merge to main.** Review and merge `stage1/atlas-loader` and
    `stage1/downstream-fusion`. New work branches from `main` with
    `<phase>/<what>` names.

0.2 **Protein splits.** Fix train / validation / test splits by protein
    (and by sequence cluster, 30% identity, so homologues do not leak).
    Anything fitted (α, calibration constants, PAE-to-covariance mapping)
    is fitted on train, chosen on validation, reported on test only.

0.3 **Held-out α.** Refit the Stage-1 shrinkage α on train proteins and
    report the downstream numbers on test proteins. Replaces the current
    in-sample α.

0.4 **Pre-register the Phase A gate.** Before running A3, commit
    `notes/phaseA_preregistration.md` with the metric, the controls, and
    the go/no-go thresholds (draft in A5). No changing them after seeing
    results; if they turn out to be wrong, say so in the write-up.

0.5 **Licences.** ATLAS is CC BY-NC (non-commercial): fine for research,
    but anything *sold* must not ship ATLAS-derived weights or data without
    permission. AlphaFold DB is CC BY 4.0, PDB is CC0. Record what each
    deliverable depends on.

0.6 **Reproducibility.** One script per result, seeds fixed, data manifests
    (entry lists, versions) committed, CI running the unit tests.

---

## Phase A: Does the big version survive? (cheap, decisive)

A1. **Subspace test on simulations (~1 hour).**
    For each ATLAS protein and replica pair (i, j):
    - error vector e = mean_j − mean_i (3N, after joint superposition),
      and per-frame errors x_t(j) − mean_i;
    - direction sources from replica i: MD covariance, ANM (MD-matched
      magnitudes), random;
    - metric: fraction of |e|² in the top-k eigenvectors, k = 5, 10, 20,
      divided by the random expectation k / (3N − 6). Also the whitened
      log-likelihood of e under each source.
    Confirms (or not) the mechanism behind the 25% fusion result without
    any fusion step.

A2. **AlphaFold error set.**
    - Proteins with an experimental structure *released after AlphaFold's
      training cutoff* (so AlphaFold has not seen the answer), plus the
      ATLAS proteins as a secondary set with that caveat stated.
    - AlphaFold DB model + PAE + pLDDT per protein (no GPU needed).
    - Map residues (SIFTS / sequence alignment), Cα only to start,
      superpose (on high-pLDDT core and on all residues; report both),
      error vector per residue.
    - Where several experimental structures exist, keep all: the spread
      between them is the experimental "noise floor" for error.
    - Target size: >= 300 proteins test, separate train/validation.

A3. **Define PAE → covariance (on train only, before looking at test).**
    PAE[i, j] is the expected error of residue j *when aligned on residue i*
    — a frame-relative quantity, which is why it connects to se(3) later.
    Candidate constructions, all fixed on train:
    - (a) relative-displacement Gaussian: treat PAE²/3 as the variance of
      x_j − x_i, embed into a 3N x 3N covariance (nearest PSD);
    - (b) PAE domains (same clustering as Phenix) as rigid bodies with
      a 6-DOF covariance each, plus pLDDT-sized per-residue noise;
    - (c) (a) with low-rank truncation.

A4. **Do AlphaFold's errors align with predicted directions?**
    Same metric as A1 on the A2 set. Direction sources:
    - PAE covariance (a, b, c), ANM (low modes and full),
      MD coupling (ATLAS subset only), random.
    Controls that coupling must beat ("not duh"):
    - **oracle scalar**: each residue's true error size, no direction;
    - **PAE rigid domains + rigid fitting** (the current practical
      standard: Phenix process/dock_predicted_model style);
    - **random subspace of the same rank** (guards against "any low-rank
      prior helps").
    Then the same comparison inside the fusion benchmark: distances
    measured on the *experimental* structure, applied to the *AlphaFold*
    model. Real errors, synthetic data; no GPU.

A5. **Go / no-go (thresholds pre-registered in 0.4; draft):**
    Go to Phase B if, on test proteins,
    - the best PAE- or ANM-derived source captures >= 2x random in the
      top 10 directions, with the CI over proteins excluding 2x; **and**
    - in the fusion test it beats both controls (oracle scalar; PAE rigid
      domains) on RMSD and lDDT, CI excluding zero, in >= 60% of proteins.
    If only MD coupling passes: coupling is real but not yet predictable
    → Phase B with step B4 (learned coupling) moved first.
    If nothing passes: Phase C.

**Phase A write-up** regardless of outcome: "Are AlphaFold's errors
low-dimensional and predictable?" This is a result either way.

---

## Phase B: Fix predicted structures (if A5 = go)

B1. **Calibrate PAE-covariance** on train: a scale and possibly a rank or
    shrinkage parameter, chosen by held-out likelihood and coverage.
    Report coverage of 50/90/95% regions and NLL on test. Same for the
    scalar baseline (calibrated pLDDT), so the comparison is fair.

B2. **Data-efficiency curves.** In the fusion benchmark, error removed vs
    number of distances (0.02–1 per residue) and noise (0.5–4 Å). The
    adoption argument is "same accuracy with N× fewer restraints".

B3. **Large-error cases.** Apo/holo and open/closed pairs from the PDB
    (curated set, e.g. from published conformational-change benchmarks);
    start from the wrong state or from the AlphaFold model, fuse distances
    from the right state. Score RMSD, lDDT, GDT, TM.

B4. **Learned coupling: DynaProt and BackFlip-2.** Their pairwise coupling
    outputs as priors in A4/B2/B3. Needs their code and weights; GPU
    machine (Adarsh's). One-command scripts prepared here first (H3).
    Outcome decides whether *our* model is needed: if they recover most of
    the MD-coupling benefit, the contribution is benchmark + calibration +
    PAE version; if not, train a coupling predictor (B4b).

B4b. **(Conditional) Train a coupling predictor** targeting the joint
    covariance (low-rank + diagonal), loss = held-out Gaussian NLL, so it
    is calibrated by construction. Inputs: structure + AlphaFold
    embeddings/PAE. Training data: ATLAS (non-commercial; see 0.5) and/or
    mdCATH.

B5. **Reconcile with the older ENM successes.** Low-mode subspace priors
    (k = 5–20 modes), other cutoffs, mass-weighted / pfANM, and a
    Zheng & Brooks 2006 style setup (few distances, known target state).
    Report where ENM does and does not work, not only where it fails.

B6. **Real experimental restraints.** Proteins with a known structure and
    deposited crosslinks (PRIDE / published XL-MS sets), NMR restraints
    (BMRB), DEER distances. Start from AlphaFold, fuse, score against the
    deposited structure. This is the result that drives adoption.

B7. **se(3) Stage 2.** Rigid-body (6-DOF) pose covariance between residue
    frames or PAE domains, transported with Ad_g. Compare against the
    3N x 3N Cα version at equal parameter count on A4, B2, B3, B6.
    Motivated by: crystal anisotropy is rigid-body (Stage 1, R4) and PAE is
    frame-relative.

B8. **Integration.** Prototype, in order of effort:
    - a restraint/prior file format + converter (AlphaFold JSON → prior);
    - ISOLDE plug-in (it already weights restraints by PAE);
    - Phenix / IMP (PMI) restraint;
    Contact the maintainers with the B6 numbers before building deeply.

B9. **Paper(s).** Main paper: PAE-derived calibrated coupling for fitting
    predicted structures to data (A + B). Companion: Stage-1 audit (C2).

---

## Phase C: Fallback results (and partly done regardless)

C1. **Fitting to experiments that measure real motion**: DEER
    distributions, NMR order parameters / RDCs, SAXS, cryo-EM
    heterogeneity (e.g. cryoDRGN outputs), where the MD-style
    25% vs 4% result applies directly.

C2. **Audit reference paper** (write in any branch of the plan):
    per-residue anisotropy real but ignorable downstream; coupling matters;
    free physics falls short; all uncertainty overconfident; MD directions
    do not match crystal data; superposition artifact.

C3. **Calibration method** for scalar and coupled uncertainty as its own
    deliverable (conformal or held-out-likelihood recalibration).

---

## Housekeeping in parallel

H1. Merge + CI (0.1, 0.6).
H2. Held-out α (0.3).
H3. GPU scripts: one command each for DynaProt and BackFlip-2 inference
    on a list of PDB/AlphaFold files, writing a standard `.npz`
    (per-residue 3x3 + pairwise coupling) that our benchmark reads.
    Include environment files and a smoke test on 2 proteins.
H4. Results registry: one table (`notes/results_index.md`) listing every
    result, the script, the commit, and the headline number.
H5. Disk hygiene: raw trajectories deleted after processing; derived
    `.npz` only.

## Order of execution

```
0.1 0.2 0.3 0.4 0.5 0.6        (housekeeping, ~1 day)
A1                              (1 hour)
A2  A3                          (A3 frozen on train before A4)
A4  → A5 gate                   (+ Phase A write-up)
  go:      B1 B2 → B3 → B4 (H3 first) → B5 → B6 → B7 → B8 → B9
  partial: B4/B4b first, then as "go"
  no-go:   C1 C2 C3
C2 drafted alongside whichever branch runs.
```
