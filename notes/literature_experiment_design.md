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
