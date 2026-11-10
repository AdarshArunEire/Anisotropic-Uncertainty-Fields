# Phase A1: do simulation errors concentrate in coupled directions?

Run 2026-11-10: `python experiments/audit/07_subspace_test.py`
(150 ATLAS proteins, Stage-1 random draw; replica pairs 0→1, 1→2, 2→0;
50 frames each). Errors are measured from replica i's mean structure, rigid
directions removed. C_k = share of error in the top-k directions of a source,
divided by the share a random k-dimensional subspace would get (1 = random).

## Results (individual frames; mean over proteins [95% CI])

| source | what it knows | C_10, all 150 | C_10, test split (34) |
|---|---|---|---|
| md_joint | coupled directions from replica i | **30.3 [28.4, 32.5]** | **29.6 [26.3, 33.3]** |
| anm | coupled directions from structure alone | **20.7 [18.7, 22.7]** | **21.3 [17.8, 25.2]** |
| md_blocks | per-residue ellipsoids, no coupling | 12.1 [11.1, 13.2] | 11.7 [9.4, 14.2] |
| md_scalar | per-residue size only | 10.3 [9.3, 11.3] | 9.9 [7.9, 12.2] |
| random | nothing | 0.99 | 0.96 |

Replica mean-shift errors give the same picture (md_joint 32.4, anm 20.9,
md_blocks 13.0, md_scalar 10.5). Full tables for k = 5, 10, 20:
`results/07_subspace_test/report.md` (regenerated, not committed).

## Reading

1. **Mechanism confirmed.** Coupled directions carry ~3x the error that a
   no-coupling scalar does (30 vs 10), and per-residue direction adds little
   on top of the scalar (12 vs 10). This is the geometry behind the 25% vs
   4% fusion result.
2. **Errors are extremely concentrated.** The top 10 of ~300–900
   directions hold about 30x their share; even a scalar that knows which
   residues are floppy gets ~10x. This is why the original "2x random"
   gate was meaningless (amendment 1).
3. **Correction to the fusion conclusion: ANM directions are good.** ANM
   reaches C_10 ≈ 21, twice the scalar, from structure alone. Yet in the
   fusion test ANM coupling removed only ~6% of error vs 25% for MD. So
   ANM points the right way but weights its directions wrongly (eigenvalue
   spectrum and scale), and the fusion step is sensitive to that weighting.
   That reconciles with the older ENM successes (Zheng & Brooks; DEN),
   which use ENM low modes as a *subspace* rather than as a calibrated
   Gaussian. It also points to a cheap fix to test: ANM eigenvectors with a
   recalibrated spectrum.

## Consequences for the plan

- A4 should include "ANM subspace, recalibrated spectrum" (fitted on train)
  as an additional source, reported as secondary (the pre-registered gate
  source is still chosen on validation from the registered list).
- Roadmap B5 (reconcile with ENM) moves earlier in Phase B.
