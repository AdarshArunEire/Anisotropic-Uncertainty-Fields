# Literature: who chooses which measurements (label sites / pairs) to make?

Searched 2026-11-10 (web search, abstracts and summaries; full texts of the
key papers still need reading before any claim of novelty). Question: who
already tells experimentalists *which* residue pairs to label or measure in
order to determine a structure, and how?

## Direct competitors: pair / site selection for structure

| Work | Experiment | How pairs are chosen | Scale of evidence | Relation to us |
|---|---|---|---|---|
| Kazmier, Alexander, Meiler, Mchaourab, *J Struct Biol* 2011 ([PMC3073550](https://pmc.ncbi.nlm.nih.gov/articles/PMC3073550)) | DEER/EPR | Pick spin-label pairs giving the most informative restraint set for Rosetta de novo folding; ~1/3 of the restraints needed | T4 lysozyme and a few proteins | Same goal (fewer, better pairs); prior is generic folding, not a predicted model's own uncertainty |
| DEERFold (Mchaourab / Meiler groups), *Nat Commun* 2025 ([article](https://www.nature.com/articles/s41467-025-62582-4), [preprint](https://www.biorxiv.org/content/10.1101/2024.10.30.621127v1.full)) | DEER into fine-tuned AlphaFold2 | Reuses the information-content pair selection above; it improved models, but **random pairs already succeeded 58.6%** of the time | Several conformational-switch proteins | **Closest competitor and natural partner.** Their selection is not driven by AlphaFold's PAE; their random baseline is strong, which is the bar |
| SpinCenter-DEER-Diffusion (Kao, Meiler, Mchaourab), preprint ([secondary summary](https://github.com/KarelBerka/Alphafoldology/issues/1161)) | DEER + diffusion model | Not checked | Not checked | To read |
| Dimura, Peulen, …, Seidel, Gohlke, *Nat Commun* 2020, FPS 2.0 ([arXiv](https://arxiv.org/pdf/1909.02148), [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7589535)) | smFRET | Automated FRET-pair selection to minimise structural uncertainty over an ensemble of sampled candidate models (greedy forward/backward, mutual-information mRMR) | T4 lysozyme (33 pairs chosen from ~13k) | **Closest in concept**: uncertainty-driven pair choice already exists for FRET. Their uncertainty comes from sampling, not from AlphaFold's PAE |
| Krishnamohan, Hamilton et al., *J Phys Chem B* 2023 ([abstract](https://www.aminer.org/pub/63db33fb90e50fcafd1c442b/coevolution-and-smfret-enhances-conformation-sampling-and-fret-experimental-design-in-tandem)) | smFRET | Coevolution + structure-based sampling of conformations, then FRET network designed on them | PDZ1–2 tandem | Same pattern: sample candidate states, then pick discriminating pairs |
| MtsslWizard ([Springer](https://link.springer.com/content/pdf/10.1007/s00723-012-0314-0.pdf)); chiLife (Tessmer & Stoll); MMM (Jeschke) | DEER | Feasibility (steric, distance window) and distance *contrast between known states* | Tools, widely used | Complementary: feasibility filter we need; contrast needs both states known |
| PCS tag-site optimisation, *Magn Reson* 2022 ([article](https://mr.copernicus.org/articles/3/65/2022/)) | NMR (pseudocontact shifts) | Optimal tag placement for one local site, rest of structure known | Method paper | Goal-directed version of the idea (a site of interest) |

## Consumers of measurements (not selectors)

These take restraints and produce better structures; they do not decide what
to measure. They are where our chosen measurements would be used.

- AlphaLink / AlphaLink2 (Rappsilber lab), crosslinks into AlphaFold2 /
  Multimer, *Nat Biotechnol* 2023 ([doi](https://www.doi.org/10.1038/s41587-023-01704-z)),
  *Nat Commun* 2024 ([AlphaLink2](https://github.com/Rappsilber-Laboratory/AlphaLink2)).
- Distance-AF, user distance constraints into AlphaFold2, *Commun Biol* 2025
  ([article](https://www.nature.com/articles/s42003-025-08783-5)).
- A general framework to inject distances into AlphaFold-like predictors,
  bioRxiv 2026 ([preprint](https://www.biorxiv.org/content/10.64898/2026.07.02.736010v1.full)).
- Experiment-guided AlphaFold as an inverse problem, arXiv 2025 ([arXiv 2502.09372](https://arxiv.org/html/2502.09372v1)).
- AlphaFold-Metainference filters AlphaFold distances by PAE ([article](https://www.nature.com/articles/s41467-025-56572-9)):
  PAE used to trust predictions, not to plan experiments.
- IMP Bayesian integrative modelling (Sali lab; review by Habeck,
  [De Gruyter](https://www.degruyterbrill.com/document/doi/10.1515/hsz-2023-0145/html?lang=en)):
  restraints as likelihood; no design module found.

## A different route: mutations that reveal structure directly

- Rollins, …, Marks, *Nat Genet* 2019 ([PMC7610650](https://pmc.ncbi.nlm.nih.gov/articles/PMC7610650))
  and Schmiedel & Lehner, *Nat Genet* 2019 ([Springer](https://link.springer.com/article/10.1038/s41588-019-0431-x)):
  epistasis in deep mutational scans identifies 3D contacts. Massive
  libraries, not a few chosen labels; competes only for "use mutations to
  get structure", and AlphaFold has largely displaced it for folds.

## General method

Bayesian optimal experimental design (expected information gain) is
standard in statistics (e.g. [arXiv 2604.21849](https://arxiv.org/pdf/2604.21849));
our greedy variance reduction is its Gaussian special case. Nothing found
that applies it with AlphaFold's PAE as the prior.

## Where the gap is (provisional)

1. **Uncertainty-driven pair selection already exists** (Kazmier 2011 for
   DEER; FPS 2.0 for FRET). The idea alone is not new; a pitch must not
   claim it is.
2. **What we did not find:** selection that uses *the predicted model's own
   error estimate* (PAE / a calibrated AlphaFold covariance), with no
   sampling step and no second state known in advance. Existing tools get
   uncertainty from sampling ensembles, a folding prior, or two known
   states.
3. **What we did not find either:** a benchmark of selection strategies on
   *real AlphaFold errors* across hundreds of proteins. Published tests are
   on one to a few proteins (often T4 lysozyme).
4. **Goal-weighted selection** (region, one coordinate, which-state) exists
   only piecemeal (PCS local-site paper; state-contrast tools).
5. **The bar is high:** DEERFold's random pairs already succeed 58.6% of the
   time. Our probe must beat random *and* an information-content baseline
   in the Kazmier/FPS style, not only pLDDT heuristics.

## Consequences for the plan

- Add two baselines to the probe: (a) **ensemble-based selection** in the FPS
  style (pairs with the largest distance spread across sampled models, e.g.
  ANM normal-mode samples), and (b) feasibility filtering (surface
  exposure) before any selection.
- Read in full before claiming novelty: Kazmier 2011, Dimura 2020,
  DEERFold 2025, SpinCenter-DEER-Diffusion.
- Natural first contacts if the probe works: the Mchaourab / Meiler groups
  (DEERFold, pair selection) and the Seidel / Gohlke groups (FPS).

---

## Full-text reading (2026-11-10)

Read via PMC / Europe PMC. DEERFold: first ~40k characters verbatim, the
rest from an extracted summary, so numbers marked (s) need checking before
citing. SpinCenter: abstract only (bioRxiv rate-limited).

### Kazmier et al. 2011 (J Struct Biol 173:549), the selection DEERFold still uses
- **Inputs: sequence, secondary structure, predicted burial only.** No 3D
  structure, no ensemble, no uncertainty.
- Score (Monte Carlo over pair sets): log sequence separation + "one restraint
  per pair of secondary-structure elements" (1:1 weight); buried residues
  excluded. A coverage heuristic, not an information measure.
- Evidence: **one protein** (T4 lysozyme C-terminal domain, 107 residues),
  Rosetta de novo folding. Correct fold (< 7.5 Å): 54.4% of models with
  optimised pairs vs 21.0% random vs 8.0% none; 21 restraints; 19 pairs
  measured by DEER for validation.
- Stated limits: cone model too soft (12 Å-wrong topologies satisfy the
  restraints); predicted secondary structure hurts; membrane proteins need
  more restraints.

### Dimura et al. 2020 (Nat Commun 11:5394), FPS 2.0, the strongest prior art
- **Genuinely uncertainty-driven.** Objective = expected Cα RMSD of the
  final model, ⟨⟨RMSD⟩⟩ = average over reference conformers of the
  RMSD of conformers weighted by the χ² p-value of their FRET distances
  (i.e. how well the data could tell them apart). Greedy forward (best below
  ~6 pairs), backward elimination, or mutual-information (mRMR) selection.
- Uncertainty comes from a **sampled ensemble** (NMSim normal-mode sampling,
  10 × 10,000 conformers; or CASP / other-state models). No AlphaFold, no
  pLDDT/PAE. Selection takes **about a day** on a desktop.
- Six proteins (148–409 residues): 10–19 guiding pairs reach Cα RMSD
  2.3–2.8 Å from seeds 4–19 Å off; T4L validated by real smFRET.
  **No random or expert baseline reported.**
- Number of pairs needed tied to model complexity:
  N_fit = 6·(N_bodies − 1) − N_bonds for rigid bodies (the se(3) count).
- Stated limits: result depends on the prior ensemble; overfitting with few
  pairs; buried and loop regions stay poor (lDDT 0.7); best for few
  subdomains.

### DEERFold (Wu, Stein, Kao, Brown, Mchaourab, Nat Commun 16:7107, 2025)
- DEER distance distributions fed into a fine-tuned AlphaFold2 (OpenFold).
  A ~2 Å-wide Gaussian per pair works about as well as the full distribution (s).
- **Pair selection = Kazmier-style heuristic**: ends of long helices
  (≥ 20/30 res) or strands, exposed, > 15 Å apart (s). **pLDDT used only as a
  post-hoc model filter; PAE not used.**
- Benchmark: 29 targets / 59 conformations, two-state transporters and
  PDBe-KB monomers with alternate states (s). Constraints needed vary hugely:
  14/29 targets ≤ 5, some ~120 (s).
- The 58.6%: Pgp, 8 random pairs out of 54 switched the model to the wide
  inward-facing state (TM ≥ 0.9) in 293/500 runs (s). Not a comparison with
  selected pairs at equal budget.
- Authors call constraint information content "a critical parameter" and
  suggest learning pair selection (s). An open invitation.

### SpinCenter-DEER-Diffusion (Kao, Meiler, Mchaourab, bioRxiv 2026-11-01)
- Abstract: differentiable spin-label model guiding a frozen Boltz-2 to
  generate weighted ensembles from full DEER distributions; synthetic tests
  plus four transporters. No pair selection mentioned in the abstract.

## What the reading changes

1. **The selection used in practice is a 2011 sequence heuristic** tested on
   one protein, still used in 2025. Nobody uses the predicted model's own
   error estimate to choose sites. The gap is real for DEER.
2. **FPS 2.0 is the real competitor for the idea** (expected-RMSD greedy
   selection), but needs a sampled ensemble and a day of compute, has no
   random baseline, and six proteins. A calibrated AlphaFold covariance gives
   the same kind of objective in closed form, in seconds: our
   tr(W·P) greedy is the Gaussian version of their ⟨⟨RMSD⟩⟩.
3. **Ready-made hard-case benchmark:** DEERFold's 29 targets / 59
   conformations are two-state proteins with deposited structures, exactly
   the "AlphaFold gives one state, the data comes from the other" setting.
4. **New deliverable idea: predict how many measurements a target needs.**
   DEERFold needed 3 to ~120; FPS ties the count to rigid-body complexity.
   The spectrum of our goal-weighted factors (or PAE domains × 6) could
   predict it in advance, which is something labs would plan budgets with.
5. **Feasibility rules to adopt:** labels on exposed residues, in helices or
   strands, > 15 Å apart; distances inside the DEER window.

## Plan changes for the probe (before it runs)

- Baselines added: (a) Kazmier/DEERFold heuristic (helix/strand ends,
  exposed, > 15 Å, maximise sequence separation and element coverage);
  (b) FPS-style ensemble selection with an ensemble sampled from the same
  prior (ANM or PAE covariance samples), greedy on expected RMSD.
- Feasibility filter (exposed, secondary structure, > 15 Å) applied to all
  strategies in a secondary analysis.
- Second benchmark: DEERFold two-state targets, scored on reaching the
  other state.
