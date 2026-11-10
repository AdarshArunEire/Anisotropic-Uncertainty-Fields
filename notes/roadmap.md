# Roadmap

Agreed 2026-11-10. Goal (see CLAUDE.md): a result that changes practice,
found honestly. Every phase ends in a written result, so nothing is wasted
if a later gate fails.

Notation: N residues; a "direction source" is any 3N x 3N covariance whose
top eigenvectors claim to say *which way* a structure is wrong or moves.


## Log (done, in order)

Times are UTC commit times. The first four commits were re-authored on
2026-11-10, which reset their timestamps to 17:16–17:17; their original
order is kept.

| # | Date / time (UTC) | Step | What happened | Headline result | Branch | Commit |
|---|---|---|---|---|---|---|
| 1 | 2026-11-10 16:58 | setup | Repository created on GitHub | — | `main` | `da67d1b` Initial commit |
| 2 | 2026-11-10 17:16 | Stage 1 scaffold | SPD primitives, calibration tests, Wishart null imported from handoff tarball | 15 tests pass | `stage1/atlas-loader` | `cd95af5` Import Stage 1 scaffold: SPD primitives, calibration tests, Wishart null |
| 3 | 2026-11-10 17:16 | Stage 1 Q1 | ATLAS loader (download, joint superposition, per-residue covariance, n_eff); Q1 first pass on 8 proteins | RMSF matches ATLAS r >= 0.998; direction reproducible 81% vs 27% control | `stage1/atlas-loader` | `d3bec63` Add ATLAS loader and run Q1 on real data |
| 4 | 2026-11-10 17:16 | housekeeping | CLAUDE.md: branch naming, commit identity, test command | — | `stage1/atlas-loader` | `5c30aed` Add CLAUDE.md: branch naming, commit identity, test command |
| 5 | 2026-11-10 17:17 | housekeeping | Co-Authored-By trailer turned off | — | `stage1/atlas-loader` | `3908d36` Turn off Co-Authored-By trailer in project settings |
| 6 | 2026-11-10 17:27 | Stage 1 audit | Dataset builder, ANM baseline, crystal ANISOU + TLS, direction-vs-amount metrics | — | `stage1/atlas-loader` | `97c486f` Add strong-audit pipeline: dataset builder, ANM baseline, crystal ANISOU + TLS, directional metrics |
| 7 | 2026-11-10 17:40 | Stage 1 audit | Results in physical units | Ellipsoid region 17% smaller than sphere (median); axis error 26° replica / 38° ANM / 60° random | `stage1/atlas-loader` | `8bb03c5` Add physical-units view of the audit (Å, degrees, coverage, volume) |
| 8 | 2026-11-10 17:45 | Stage 1 audit | Full audit, 150 random proteins + 99 with crystal ANISOU ([stage1_audit.md](stage1_audit.md)) | Direction real (72% vs 32% control); raw anisotropy overconfident; ANM closes 47%; MD ≠ crystal (rigid body) | `stage1/atlas-loader` | `962b34f` Stage 1 audit results: 150 random ATLAS proteins + 99 with crystal ANISOU |
| 9 | 2026-11-10 17:59 | downstream | Fusion benchmark (sparse distances + prior covariance) | — | `stage1/downstream-fusion` | `9c5ce43` Add downstream fusion test: distance data + prior covariance |
| 10 | 2026-11-10 18:00 | housekeeping | Project goal and honesty bar in CLAUDE.md | — | `stage1/downstream-fusion` | `db258c2` CLAUDE.md: record project goal and honesty bar |
| 11 | 2026-11-10 18:05 | downstream | Fairness control: ANM coupling at MD magnitudes | — | `stage1/downstream-fusion` | `28a910a` Fusion test: add ANM couplings at MD per-residue magnitudes |
| 12 | 2026-11-10 18:10 | literature | Review of coupling priors for data fusion ([literature_coupling.md](literature_coupling.md)) | Benchmark, calibration, PAE→covariance appear open; DynaProt/BackFlip-2 already predict coupling | `stage1/downstream-fusion` | `0cb9c13` Add literature review on coupling priors for data fusion |
| 13 | 2026-11-10 18:11 | downstream | Fusion results, 40 proteins ([downstream_fusion.md](downstream_fusion.md)) | Coupling removes 25% of error vs 4–5% per-residue; ANM coupling 6%; all priors overconfident | `stage1/downstream-fusion` | `9b871a3` Downstream fusion results: coupling matters, per-residue anisotropy does not |
| 14 | 2026-11-10 18:18 | planning | Roadmap agreed (this file) | — | `stage1/downstream-fusion` | `26b4130` Add roadmap: phases 0/A/B/C with gates, controls and order |
| 15 | 2026-11-10 18:19 | housekeeping | Chronological log + step status; maintenance rules in CLAUDE.md | — | `stage1/downstream-fusion` | `e9dc503` Roadmap: add chronological log and step status; CLAUDE.md: how to maintain it |
| 16 | 2026-11-10 18:21 | 0.2 | Homology-aware splits for all ATLAS entries ([splits/](../splits/README.md)) | train 1166 / val 377 / test 395; random-150 → 83/33/34 | `phase0/splits-and-preregistration` | `075db44` Phase 0.2: homology-aware ATLAS splits (ECOD H-group components) |
| 17 | 2026-11-10 18:21 | 0.3 | Shrinkage α fitted on train, scored on test | α = 0.55; test: raw −0.118, shrunk +0.155 [0.092, 0.220] nats/res/frame | `phase0/splits-and-preregistration` | `fe38372` Phase 0.3: held-out shrinkage alpha (train-fitted 0.55; test gain +0.155 nats) |
| 18 | 2026-11-10 18:22 | 0.4 | Phase A pre-registration ([phaseA_preregistration.md](phaseA_preregistration.md)) | — | `phase0/splits-and-preregistration` | `5388611` Phase 0.4: pre-register Phase A metrics, controls and go/no-go gate |
| 19 | 2026-11-10 18:22 | 0.6 | CI: pytest on every push | — | `phase0/splits-and-preregistration` | `886f2a8` Phase 0.6: run the test suite on GitHub Actions |

## Status of planned steps

| Step | Status | Branch | Notes |
|---|---|---|---|
| 0.1 merge to main | waiting for Adarsh's review | — | `stage1/atlas-loader`, `stage1/downstream-fusion` |
| 0.2 splits | done | `phase0/splits-and-preregistration` | `splits/atlas_splits.tsv` |
| 0.3 held-out α | done | `phase0/splits-and-preregistration` | α = 0.55 |
| 0.4 pre-registration | done | `phase0/splits-and-preregistration` | |
| 0.5 licences | done (noted in roadmap) | — | |
| 0.6 CI | done | `phase0/splits-and-preregistration` | GitHub Actions |
| A1 subspace test | in progress | `phaseA/subspace-test` | |
| A2–A5 | not started | — | |
| B, C | not started | — | |

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
