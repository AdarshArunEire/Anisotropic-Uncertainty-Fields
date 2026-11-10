# Stage 1 audit: 150 ATLAS proteins (23480 residues), random draw
Mean over proteins, [95% bootstrap CI over proteins].

## R1  Reproducibility of direction (kill criterion)
Sigma_i beats iso(Sigma_i) at predicting Sigma_j:  72% [69, 75] (n=150)
  same, isotropic truth through same superposition: 32% [30, 35] (n=150)
  margin over control:                              40% [36, 43] (n=150)
  proteins clearing control by >10 pts:             132/150
  fit on starting structure:                        72% [70, 75] (n=150)
  fit on rigid half only:                           79% [76, 81] (n=150)

## R2  How much anisotropy is real, and what it is worth
fractional anisotropy (median per protein):        0.65 [0.63, 0.66] (n=150)
optimal shape power alpha, replica -> replica:     0.74 [0.71, 0.76] (n=150)
principal-axis |cos|, replica -> replica (rand .5): 0.77 [0.76, 0.79] (n=150)
held-out log-lik gain over scalar, nats/res/frame: -0.055 [-0.134, 0.018] (n=150)
  same with alpha-shrunk Sigma:                     0.167 [0.144, 0.193] (n=150)

## R3  Predictable from structure alone? (ANM, no learning)
ANM direction beats scalar vs MD:                  66% [62, 69] (n=150)
ANM optimal shape power:                           0.95 [0.88, 1.01] (n=150)
ANM principal-axis |cos| vs MD:                    0.69 [0.67, 0.70] (n=150)
median unit-trace AIRM to replica j (lower is better):
  scalar (floor)         1.32 [1.26, 1.38] (n=150)
  ANM raw                1.26 [1.17, 1.36] (n=150)
  ANM, alpha tuned       1.10 [1.05, 1.16] (n=150)
  replica i (ceiling)    0.98 [0.91, 1.05] (n=150)
fraction of floor->ceiling gap closed by ANM:      0.47 [0.40, 0.56] (n=150) (median)

## R4  Against experiment: crystallographic ANISOU (104 proteins)
FA crystal vs FA MD:                               0.32 [0.29, 0.34] (n=104) vs 0.65 [0.63, 0.67] (n=104)
MD   beats scalar vs crystal:                    7% [5, 9] (n=104)
     optimal shape power:                         0.11 [0.09, 0.13] (n=104)
     principal-axis |cos| (random 0.5):           0.58 [0.56, 0.60] (n=104)
ANM  beats scalar vs crystal:                    20% [17, 23] (n=104)
     optimal shape power:                         0.20 [0.17, 0.24] (n=104)
     principal-axis |cos| (random 0.5):           0.59 [0.57, 0.60] (n=104)
TLS (1 rigid group, 20 params, in-sample) beats scalar: 67% [63, 71] (n=104)
     principal-axis |cos|:                        0.80 [0.78, 0.82] (n=104)
TLS + MD shape fits better than TLS + scalar: 24/104 proteins; median RSS ratio 1.04 [1.02, 1.05] (n=104)
TLS + ANM shape fits better than TLS + scalar: 31/104 proteins; median RSS ratio 1.02 [1.01, 1.03] (n=104)

## R5  Where the signal lives (residue level, pooled; descriptive)
stratum               n  rep win  control  axis cos   loglik  ANM win
helix              9927      73%      30%      0.78    0.011      68%
strand             4452      73%      33%      0.79    0.055      66%
coil               9101      73%      28%      0.77   -0.084      67%
termini (5)        1500      62%      28%      0.73   -0.476      52%
interior          21980      74%      30%      0.78    0.014      68%
RMSF low third     7827      76%      32%      0.80    0.095      67%
RMSF mid third     7826      78%      33%      0.80    0.189      77%
RMSF high third    7827      65%      24%      0.73   -0.336      57%
