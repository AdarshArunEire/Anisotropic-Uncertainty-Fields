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
experiments/audit/ Runnable Stage-1 scripts
notes/             Literature review, design notes
results/           Generated figures and reports (not committed)
tests/             Unit tests
```

## Status

Scaffold. Nothing measured yet.

## License

MIT.
