# Phase B probe: choosing which distances to measure

Run 2026-11-12 on `phaseB/planning-probe`, test split (266 proteins), criteria in
[phaseB_design_preregistration.md](phaseB_design_preregistration.md)
(amendments 1–2, committed before any result). Full tables:
[results/phaseB/b1_report_test.md](../results/phaseB/b1_report_test.md);
curve: [b1_curve_test.png](../results/phaseB/b1_curve_test.png).

## Verdict

**Pre-registered criterion: FAIL.** At 20 measured pairs, choosing pairs from
AlphaFold's own uncertainty (greedy_pae) does not raise lDDT-Cα more than
random, the pLDDT-scalar greedy or the DEER heuristic. But **no strategy raises
lDDT at all**: every one lowers it, and more measurements lower it more.
The headline metric failed because of the fusion step, not because of the
choice of pairs.

On the secondary metric (Cα RMSD removed), the choice of pairs matters clearly.

## Results (test, m = 20, mean [95% CI])

Mean AF error 1.97 Å (median 1.04 Å); median AF lDDT-Cα 95.7.

| strategy | Cα RMSD removed (Å) | lDDT-Cα change (points) |
|---|---|---|
| random | +0.004 [−0.011, +0.019] | −1.21 |
| pLDDT pairs | +0.033 | −1.22 |
| DEER heuristic (Kazmier / DEERFold rule) | +0.114 | −1.51 |
| ensemble mRMR (Hays 2018 / FPS 2.0 style) | +0.132 | −1.35 |
| greedy, pLDDT scalar covariance | +0.120 | −1.12 |
| **greedy, PAE covariance (ours)** | **+0.183 [+0.143, +0.225]** | −1.71 |
| greedy, ANM, fused with its own prior | +0.215 | −0.81 |

Paired, greedy_pae minus the alternative, RMSD removed: vs random +0.18 Å
(better in 82% of proteins), vs pLDDT pairs +0.15 (71%), vs DEER heuristic
+0.07 [+0.003, +0.12] (70%), vs mRMR +0.05 (60%), vs pLDDT-scalar greedy +0.06
(58%).

**Post hoc, by AlphaFold's own error** (not part of the criterion):

| AF error | n | random | DEER heuristic | greedy_pae | greedy_pae vs heuristic |
|---|---|---|---|---|---|
| ≤ 1 Å | 130 | −0.03 | −0.02 | +0.01 | +0.03 Å, 61% better |
| 1–2 Å | 78 | +0.005 | +0.04 | **+0.22** | **+0.18 Å, 91% better** |
| > 2 Å | 58 | +0.08 | +0.51 | +0.52 | +0.005 Å (tie) |

Where the model is moderately wrong (1–2 Å), uncertainty-guided pairs remove
about five times more error than the DEER rule, while random pairs remove
essentially none.
Where it is badly wrong (> 2 Å), the DEER rule (long-range pairs, 15–60 Å)
does as well; ANM with its own prior is best there (+0.68 Å).

**Predicted vs realised error (the "how many pairs" curve):** across proteins
the predicted error ranks the realised error with Spearman ρ ≈ 0.62–0.65, and
the realised error is a **constant ~25% of the predicted** at every budget
(0, 5, 10, 20, 50 pairs). The curve's shape is right and its scale is off by
one factor (the frozen covariance scale was fitted by likelihood on
heavy-tailed errors, which inflates it). One calibration constant fitted on
train would fix the scale; not done here.

## Reading

1. **Pair choice is real.** Random measurements fix essentially nothing;
   20 well-chosen pairs remove ~0.2 Å on average and ~0.5 Å where AlphaFold
   is > 2 Å off. AlphaFold-uncertainty selection beats today's practice
   (random, the DEER heuristic) on RMSD, decisively in the 1–2 Å group.
2. **The fusion step is the bottleneck.** Cα-only Gaussian fusion of 1 Å
   distances lowers local accuracy (lDDT) for every strategy, even with chain
   restraints, and the better a strategy moves the structure, the more it
   tends to cost in lDDT. Real use would fuse with a method that keeps local
   geometry by construction (DEERFold / AlphaLink style, or all-atom
   refinement). Whether the RMSD advantage survives that is untested.
3. **Coupling helps modestly over a scalar greedy** (+0.06 Å overall, +0.24 Å
   for > 2 Å errors), and the closed-form PAE greedy beats sampling-based mRMR
   at equal prior (+0.05 Å).
4. **The size is small in absolute terms** on this mostly well-predicted set:
   about 9% of the mean error at 20 pairs. The use case is proteins
   AlphaFold gets moderately or badly wrong.

## What it would take to continue (not agreed)

- Fusion that preserves local geometry (more local restraints, or a
  DEERFold/AlphaLink-based fusion on a GPU), then re-test with lDDT.
- Calibrate the predicted curve on train (one constant), check on test.
- Two-state hard cases (DEERFold targets), where the 1–2 Å and > 2 Å pattern
  would be tested properly.
