# Literature, coupling priors for data fusion (2026-11-10)

Question: has anyone used a (predicted) joint inter-residue covariance as a
prior for fusing structures with sparse experimental data, and benchmarked
it downstream against scalar / per-residue / ENM priors?

Compiled by a literature-search agent; "verified" means the abstract, full
text or docs were fetched. Not exhaustive.

## Verdict

- **Done:** ENM / normal modes as a deformation prior for sparse distances or
  density; PAE-weighted restraints and PAE rigid-domain splitting;
  independent per-pair Gaussian restraints from AF2 distograms; learned
  coupling prediction (DynaProt, BackFlip-2).
- **Partly done:** ENM + sparse distances on a handful of proteins
  (Zheng & Brooks 2006); Kalman structure estimation for NMR (DIKF, 1990s,
  lost to restrained MD).
- **Apparently open (medium confidence):**
  1. A many-protein, ground-truth benchmark of scalar vs per-residue
     anisotropic vs ENM vs learned / MD joint covariance as Gaussian priors
     for sparse-data fusion.
  2. Calibration (coverage, likelihood) of predicted joint covariance or PAE.
  3. PAE converted into a proper joint or se(3) covariance.
  4. Residue-frame / domain se(3) pose covariance in proteins (mature in
     robotics: Barfoot & Furgale 2014).

## Things that change our plan

- **DynaProt and BackFlip-2 already predict pairwise coupling** (DynaProt:
  pairwise scalar covariances; BackFlip-2: pairwise dynamic correlations).
  The novelty is not "predict coupling"; it is whether predicted coupling is
  calibrated and useful downstream. Their coupling heads must be tested as
  priors in our fusion benchmark.
- **Our ANM-null result conflicts with older successes** (Zheng & Brooks 2006;
  DEN). Likely differences: they use a few low modes as a subspace from the
  starting structure toward a target state, not a full calibrated Gaussian.
  Needs checking before claiming "cheap physics fails".
- **Kalman / Gaussian conditioning lost to MD sampling once** (DIKF, 1992).
  Reviewers will raise it; the answer has to be accuracy plus calibration at
  negligible cost, shown on a benchmark.

## Key works

| Work | Link | What | Scalar/isotropic control? | Downstream vs truth, many proteins? |
|---|---|---|---|---|
| Schröder, Levitt & Brunger 2010 Nature (DEN) | https://pmc.ncbi.nlm.nih.gov/articles/PMC2859093 | ENM-like spring prior, low-res crystallography | No (noDEN only) | Rfree, 19 entries; 1 synthetic vs truth |
| Zheng & Brooks 2006 Biophys J | https://pmc.ncbi.nlm.nih.gov/articles/PMC1471861 | ANM modes + <= 10 sparse distances | No | Several proteins, RMSD to end state |
| Altman / Jardetzky DIKF; Liu et al. 1992 | not re-fetched | Kalman NMR structure estimation | — | Lost to restrained MD |
| López-Blanco & Chacón 2013 (iMODFIT) | PubMed 23999189 | NMA cryo-EM flexible fitting | No | Limited |
| Alshammari, He & Wriggers 2024 Bioinf Adv | https://pmc.ncbi.nlm.nih.gov/articles/PMC11783307 | ENM fitting of AF2 models to EM maps | No | 4 proteins |
| Oeffner et al. 2022 Acta D (Phenix + ISOLDE) | https://journals.iucr.org/d/issues/2022/11/00/ai5009/ | pLDDT -> B, PAE domains, PAE-weighted restraints | No | No |
| bAIes (Bonomi lab, bioRxiv 2025) | https://doi.org/10.1101/2025.06.25.661454 | Distogram -> per-pair Gaussian restraints, Bayesian | No | Task-specific |
| Jussupow & Kaila 2023 JCTC (AF-ENM) | https://pmc.ncbi.nlm.nih.gov/articles/PMC11181330 | ENM from AF scores | n/a | Dynamics only |
| Sweeney et al. 2024 MLSB | MLSB 2024 proceedings | GP covariance head on AF2 | Full vs diagonal mentioned | No fusion |
| DynaProt (arXiv:2509.01038) | https://arxiv.org/abs/2509.01038 | 3x3 + pairwise scalar coupling | — | No fusion |
| BackFlip-2 (arXiv:2609.08474) | https://arxiv.org/abs/2609.08474 | Directional flexibility + pairwise correlations | — | No fusion |
| Stahl et al. 2023 Nat Biotech (AlphaLink) | https://doi.org/10.1038/s41587-023-01704-z | Crosslinks into AF2 | n/a | Yes, no covariance prior |
| Maddipatla et al. (arXiv:2502.09372) | https://arxiv.org/abs/2502.09372 | AF3 as prior guided by experiment | n/a | Yes, implicit prior |
| Barfoot & Furgale 2014 IEEE T-RO | https://doi.org/10.1109/TRO.2014.2298059 | se(3) pose uncertainty and fusion | — | Robotics |

Not re-checked by the agent: Tama 2004 NMFF, HEMNMA, TLS papers, AF2 PAE
calibration plots, CalPro (arXiv:2601.07201) numbers, 2026 journal versions
of ROCKET and experiment-guided AF3, resTrain (bioRxiv 2026). Possible
missed area: FRET-positioning (Seidel / Kalinin FPS) and IMP Gaussian priors.

## Likely adopters

ISOLDE and Phenix (already consume PAE); integrative modelling (IMP/PMI,
bAIes/PLUMED, HADDOCK); crosslink-MS and DEER/FRET modelling groups.
