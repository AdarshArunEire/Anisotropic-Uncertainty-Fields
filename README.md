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

ATLAS loader written and validated (reproduces ATLAS's published RMSF at
r >= 0.998). First pass of Q1 on 8 randomly drawn ATLAS entries: the
anisotropy is directional and reproducible across replicas, well above a
superposition-artifact control, so the Stage 1 kill criterion is not
triggered. See [notes/stage1_q1_results.md](notes/stage1_q1_results.md).
Q2 and Q3 need the published models' predictions and are not started.

## Running

```
pip install -e ".[test,atlas]"
pytest
python experiments/audit/01_information_content.py            # downloads ~250 MB of ATLAS
python experiments/audit/01_information_content.py 1ng6_A    # any ATLAS pdb_chain
```

ATLAS data are CC BY-NC 4.0. They are downloaded into `data/` and never
committed; derived numbers in `notes/` cite the source.

## License

MIT.
