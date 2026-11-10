# Literature, Stage 1

## Prior work that predicts anisotropic covariance

**DynaProt** (Bafna, Jing, Berger, Sep 2026, arXiv:2509.01038). Per-residue
3x3 Ca covariance via Cholesky factor, plus scalar pairwise couplings. Trained
on ATLAS MD (1390 proteins). SE(3)-*invariant* IPA backbone. The mean is not
learned; it is fixed to the input coordinate. Joint 3Nx3N covariance is
reconstructed by a heuristic the paper calls "not exact". RMSF Pearson r
median 0.865 vs FlexPert-3D 0.830 and NMA 0.697.

**BackFlip-2** (Viliuga et al., Sep 2026, arXiv:2609.08474). Same output
structure but SE(3)-*equivariant*: Sigma = R_i Sigma_loc R_i^T. Trained on
ATLAS, evaluated on ATLAS test, mdCATH, NMR ensembles, and de novo designs.
RMSF per-target r 0.87.

## What neither paper does

1. **No calibration of any kind.** Metrics are Pearson r, MAE, RMWD,
   symmetric KL, ellipsoid IoU/Dice. No coverage, chi-squared, NEES or
   likelihood. The covariances are scored as point predictions of a
   covariance, never as uncertainty.
2. **No isotropic baseline.** BackFlip-2's nearest equivalent is an
   equivariant-vs-invariant ablation (IoU 0.49 vs 0.42), where the invariant
   model is noted to produce "nearly isotropic" covariances. That is an
   architecture ablation, not a scalar baseline.
3. **No directional metric separated from magnitude**, and no
   sampling-noise ceiling on either.
4. **3x3 positional only.** No rotational block, no rotation-translation
   coupling, so nothing in se(3).
5. **Target is MD fluctuation about equilibrium, not prediction error**
   against experiment. pLDDT/PAE calibration is not discussed in either.

## Evidence that anisotropy may not pay

This is the risk, and it is real.

- BackFlip-2's gain over its own scalar predecessor BackFlip on RMSF is
  r 0.87 vs 0.84, MAE 0.58 vs 0.62 A. Small.
- Its couplings, the relational/directional part, are its weakest result:
  DCCM r 0.80 against AFMD-T's 0.89.
- arXiv:q-bio/0611039 compared computed anisotropic displacement tensors to
  crystallographic data: magnitude agreement ~0.66-0.68 correlation across
  methods, but *directional* agreement varied considerably.
- Crystallography has litigated this for thirty years. Merritt
  (PMC3322606) argues the isotropic / anisotropic / TLS choice must be settled
  by statistical analysis of refinement, not a resolution heuristic, and that
  combining treatments "is not always statistically justifiable". Parameter
  counting binds: isotropic B gives ~2 observations per parameter, a full
  ellipsoid adds six per atom. Afonine et al. (IUCr D 2018, rr5149) warn
  unrestrained TLS refinement does not guarantee physically realistic
  parameters. TLS works because it is few parameters per rigid *group*.

## Evidence that it does pay, elsewhere

- MomentBA (arXiv:2609.13691) reports lower rotational error for anisotropic
  correspondence uncertainty than isotropic models in bundle adjustment.
- GARPS (arXiv:2509.13652) reports a significant drop when its anisotropic
  covariance is replaced by identity.
- Both are self-ablations, not independent benchmarks. Note also that in these
  settings the covariance is consumed by a *fusion* step, so the decision is
  directional by construction. That is the condition to look for: anisotropy
  pays when the downstream decision is directional.

## Design consequence

Make the contribution the audit, not the win. The published gap is an
evaluation gap, and an audit is informative whichever way the measurement
falls. If anisotropy turns out to carry little signal beyond its trace, that
is a finding about 21 parameters per residue buying one parameter's worth of
information.
