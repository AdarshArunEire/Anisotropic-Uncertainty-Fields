# Anisotropic Uncertainty Fields

Per-residue uncertainty in protein structure is a direction-dependent object: a
covariance on the tangent space, not a scalar. Recent models predict such
covariances ([DynaProt](https://arxiv.org/abs/2509.01038),
[BackFlip-2](https://arxiv.org/abs/2609.08474)) but evaluate them with
correlation and distribution distances only. Neither reports calibration, and
neither reports an isotropic baseline.

This repository works through that in stages.

## Stage 1 — Audit

Does the anisotropy carry signal, and are predicted covariances calibrated?

Three questions, answerable on existing data with no new model:

1. **Information content.** How much is lost by replacing each ground-truth
   covariance with `(tr Σ / 3) I`, against a null that accounts for
   finite-trajectory sampling noise?
2. **Isotropic baseline.** The comparison the published papers omit.
3. **Calibration.** Do standardised residuals follow their reference
   distribution (χ²₃ for 3×3 positional covariance), or are the ellipsoids
   systematically the wrong shape?

A negative result here is informative: it would mean the field is fitting 21
parameters per residue for one parameter's worth of signal.

Later stages (6D se(3) covariance including rotation–translation coupling;
transport between residue frames) are gated on Stage 1 finding directional
signal.

## Layout

```
src/auf/audit/     Stage-1 library code
src/auf/data/      ATLAS loader (download, superposition, covariances, n_eff)
experiments/audit/ Runnable Stage-1 scripts
notes/             Literature review, design notes
results/           Generated figures and reports (not committed)
tests/             Unit tests
```

## Status

Stage 1 audit run on 150 random ATLAS proteins plus 99 with crystal ANISOU.
Full write-up: [notes/stage1_audit.md](notes/stage1_audit.md).

- Anisotropy is real and reproducible across MD replicas (kill criterion
  not triggered), but worth something only after shrinkage toward isotropy.
- A one-parameter elastic network gets ~47% of the way from scalar to the
  replica ceiling, a baseline neither published model reports.
- MD direction does not match experimental crystal anisotropy, which is
  mostly whole-molecule rigid-body motion.

Q2/Q3 need the published models' predictions and are not started.

## Running

```
pip install -e ".[test,atlas]"
pytest
python experiments/audit/00_build_dataset.py --n 150   # ~30 min, raw data deleted as it goes
python experiments/audit/02_audit.py
python experiments/audit/03_figures.py
python experiments/audit/04_physical_units.py
```

ATLAS data are CC BY-NC 4.0. They are downloaded into `data/` and never
committed; derived numbers in `notes/` cite the source.

## License

MIT.
