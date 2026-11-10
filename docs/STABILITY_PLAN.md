# Stability analyses: plan (draft for review, 2026-10-24)

## Why

Each universe in the main multiverse is one fit. The spec curve shows how results
change *across* analytic choices, but not whether a given solution would come back
if we refit it, resampled the children, used other sites, nudged the method's
settings, or added realistic measurement noise. This section of the paper asks
three questions:

1. **Are the classes real?** Do the same subgroups reappear across seeds,
   resamples, site halves, settings and noise, or do they dissolve?
2. **Optimization reliability.** Is 20 starts enough? Do repeated starts land on
   the same solution?
3. **Null benchmark.** Do the methods "find" tidy classes in data that has no true
   subgroups?

The poverty RRR is tracked as a secondary outcome in every replicate.

Placement: **main paper**, as a stability section next to the spec curve.

## Anchors

The stability work runs on two anchor configurations × all 7 methods × k = 2-6:

| Anchor | Scale | Family | Sample | Time (GBTM/GMM) |
|---|---|---|---|---|
| A-raw | raw total | one per family | ≥3 waves, imputed | age − 9 |
| A-log | log(1+raw) | one per family | ≥3 waves, imputed | age − 9 |

These two capture the contrast that drove the main result: raw scores split off a
small severe class, while the log scale makes the "highest" class a large moderate
group. n = 8,331.

## Perturbations and replicate counts (tiered)

| Perturbation | What changes | Cheap methods (k-means, PAM, Ward, Gaussian mixture, kml) | GBTM / GMM |
|---|---|---|---|
| **Seeds / starts** | Same data and settings, new seed for the whole procedure (20 starts each) | 50 seeds | All 20 gridsearch starts recorded per fit (log-likelihood + solution), not extra seeds |
| **80% subsamples** | Random 80% of children, without replacement | 100 | 20 |
| **Site split-half** | 21 sites randomly split into halves; fit on each half and predict the other | 20 splits | 20 splits |
| **Measurement noise** | Jitter scores by their measurement error at two reliability levels | 20 per level | 5 per level |
| **Settings** | Method hyperparameters (below) | full-data fit + 20 subsamples per variant | full-data fit per variant |
| **Null benchmark** | Simulated data with no subgroups | 20 null datasets | 10 null datasets |

Notes:
- **Deterministic methods.** Ward and the Gaussian mixture don't depend on the
  seed, so they get no seed replicates; that is reported, not hidden. PAM at
  n > 2,000 uses CLARA, which is stochastic, so it does.
- **Seeds for GBTM/GMM.** Instead of whole refits, I'll rewrite the gridsearch to
  keep every start. For each fit we then report the share of starts reaching the
  best log-likelihood (within 0.01) and the ARI between the top distinct solutions.
  This is the standard "best log-likelihood replicated" check from the latent-class
  literature.
- **Site splits.** Uses prediction strength (Tibshirani & Walther, 2005). For each
  split, fit on half A, assign half B's children to A's classes, and compare with
  B's own fit (ARI); then repeat in the other direction. Children are assigned to
  another fit's classes by:
  - k-means, kml, Ward: nearest class mean trajectory
  - PAM: nearest medoid
  - Gaussian mixture: `predict.Mclust`
  - GBTM/GMM: posterior class from the half-A model (lcmm `predictClass`)
- **Noise.** The measurement error (SEM) is SD_wave × √(1 − r), on the raw scale.
  Jittered raw scores are rounded and floored at 0; the log anchor takes
  log1p of the jittered raw. The two levels are:
  1. the published test-retest reliability for CBCL Total Problems. The October run
     used .94; **.90 since 2026-11-07** (see "Open items");
  2. r = .80.

## Settings variations (stability only; the main spec curve is unchanged)

| Method | Base | Variants |
|---|---|---|
| Hierarchical | Ward.D2 | average linkage, complete linkage |
| PAM | Euclidean | Manhattan |
| Gaussian mixture | VVV | VVI (diagonal), EEE (shared covariance) |
| GBTM | quadratic | linear |
| GMM | quadratic, random intercept + slope | random intercept only (quadratic); linear with random intercept + slope |

Each variant is compared with the base solution on the full anchor data (ARI).
Cheap variants also get 20 subsamples, so their own stability can be compared
with the base method's.

## Null benchmark

- **Generating model.** A single-class growth model (quadratic in age, random
  intercept + slope) fitted to the anchor children on the **log(1+raw)** scale.
- **Simulation.** Each null dataset keeps the real children's ages, sex, site and
  income. CBCL trajectories are simulated from the 1-class model: random effects
  plus residual.
- **Log anchor.** Uses the simulated log scores directly.
- **Raw anchor.** Back-transforms with expm1, then rounds and floors at 0, so the
  null keeps the real floor and long right tail but has no subgroups.
- **Analysis.** The same pipeline as the real data: all 7 methods, k = 2-6, the
  same fit indices, relabelling and poverty model.
- **Comparison with real data:**
  - silhouette, relative entropy, and BIC gain from k to k+1;
  - which k each criterion selects;
  - class sizes and shapes (does the null also produce a small "severe" class?);
  - stability measured the same way (ARI across null replicates).
- **Poverty false-positive check.** Income is unrelated to the simulated CBCL, so
  the poverty RRR should sit at 1. Its spread across null datasets gives the
  false-positive rate.

## Headline metrics

1. **Solution-level ARI.** Mean (and IQR) ARI between replicate solutions, by
   method × k × perturbation × anchor:
   - subsamples: compared on the children each pair of replicates shares;
   - site splits: prediction-strength ARI;
   - seeds and noise: compared on all children.
2. **Child-level consistency.** After the usual relabelling (class K = highest mean
   CBCL), take the share of replicates in which each child lands in the highest
   class, out of the replicates that include that child. Reported as:
   - the distribution per method × k;
   - % "reliably high" (≥ .90);
   - % "never high" (≤ .10);
   - % on the boundary (in between).

Recorded but not headlined: per-class Jaccard, k-selection frequency, and fit
indices, so any of them can be added later without refitting.

**Poverty (secondary).** The RRR (highest vs lowest, per halving of income) is
stored for every replicate. Its spread is reported next to the main CI, alongside
the null false-positive distribution.

## Figures (main paper unless noted)

| Fig | Content |
|---|---|
| Stab-1 | Stability heatmap: method × k, mean ARI. Facets for perturbation (seeds, subsample, site split, noise r = .90 and .80) and anchor |
| Stab-2 | Child-level consistency: distribution of P(highest class) per child, method × k, raw vs log; % reliably high / boundary |
| Stab-3 | Null benchmark: real vs null silhouette, entropy and BIC gain by k and method; selected k in real vs null |
| Stab-4 | Poverty RRR: spread across subsamples and site splits next to the main CI, with the null distribution around 1 |
| Stab-S1 (supp) | Optimization: per-start log-likelihoods for GBTM/GMM; share of starts reaching the best |
| Stab-S2 (supp) | Settings: ARI of each variant with the base, and variant vs base subsample stability |

## Compute (local, 12 cores)

GBTM/GMM dominate. Per-task times from the full run (anchor, n = 8,331, k = 2-6):
GMM ≈ 73 min, GBTM ≈ 24 min. These scale roughly with n.

| Block | GMM CPU-h | GBTM CPU-h |
|---|---|---|
| Subsamples (20 × 2 anchors, n ≈ 6,700) | ~40 | ~13 |
| Site splits (20 × 2 halves × 2 anchors, n ≈ 4,200) | ~49 | ~16 |
| Null (10 × 2 anchors) | ~24 | ~8 |
| Noise (5 × 2 levels × 2 anchors) | ~24 | ~8 |
| Settings (3 GMM + 1 GBTM variants × 2 anchors) | ~7 | ~1 |
| Start recording (2 anchors) | ~2.5 | ~1 |
| **Total** | **~147** | **~47** |

That's about 195 CPU-hours, or **roughly 16-20 hours of wall time on 12 cores**.
The cheap methods add about 1-2 hours. As with the main run, the Mac needs to stay
plugged in with the lid open. I'll time one replicate of each type before
launching the whole thing.

Each GBTM/GMM replicate (one anchor × one perturbation × one replicate, all k) is
saved on its own, so an interruption loses at most one replicate (about 1 hour).

## Implementation

| File | Role |
|---|---|
| `07_stability_spec.R` | Replicate table: anchor × method × perturbation × replicate, with its seed |
| `R/perturb.R` | Subsample, site-split, noise and null-data generators; settings variants |
| `R/predict_class.R` | Assign new children to an existing fit's classes, per method |
| `R/fit_task.R` | Extended: settings arguments, and gridsearch that records every start |
| `08_stability_run.R` | Parallel, restartable, atomic saves; same pattern as 03 |
| `09_stability_summary.R` | ARI, child-level consistency, Jaccard, k selection, poverty RRR per replicate |
| `10_stability_figures.R` | Figures Stab-1 to Stab-4 and S1-S2 |

## As built (2026-10-24, second Mac: M5, 10 cores)

Run order:
```
Rscript 07_stability_spec.R full          # replicate table, null model, base solutions
Rscript 08_stability_run.R full 8         # restartable; filters e.g. anchor=raw,maxrep=1
Rscript 09_stability_summary.R full 8     # output/stability/summary.rds + output/summaries/summary_stability_full.txt
Rscript 10_stability_figures.R full       # output/figures/3_stability/
```
`smoke` in place of `full` runs the whole pipeline on 600 children with 2 replicates per type
(output/stability_smoke/), for checking code changes.

Changes from the plan above:
- **Sites: 22, not 21** (site22 has 21 anchor children). Splits are 11 vs 11 sites.
- **The Gaussian mixture is not deterministic at full n.** mclust initialises on a random
  2,000-child subset when n > 2,000, so it now gets 50 seed replicates like k-means, PAM
  and kml (2,860 replicates in all). `fit_llpa` now sets the seed; the main run's
  full-sample mixture fits were unseeded (see README flags).
- **Start recording converges the top 5 starts**, not only the best one. Every one of the
  20 starts' 10-iteration log-likelihood and classes is kept, and the 5 best are run to
  convergence, which is what the "best log-likelihood replicated" check needs. This adds
  about 4 full fits per anchor × method × k (~13 CPU-h). `N_CONVERGE` in R/perturb.R.
- **Null stability** is measured on 20 subsamples of null dataset 1 (cheap methods). ARI
  *between* null datasets would be ~0 by construction, since each null dataset has
  independent CBCL draws. GBTM/GMM null stability is not run (cost).
- **Null log anchor** uses log1p of the rounded raw draw (not the unrounded simulated log
  score), so both anchors cluster the same integer counts.
- **Paired design.** Subsamples, site halves, noise draws and null datasets are seeded by
  perturbation × replicate only, so all methods and both anchors see the same data.
  Settings-variant subsamples are subsamples 1-20.
- **SEM** from the anchor's wave SDs (17.3-17.9): 4.3 raw points at r = .94 (October run),
  5.5-5.7 at r = .90 (current), 7.9 at r = .80.

## Rerun (2026-11-07 to 11-08)

- Noise level 1 is r = .90 (`NOISE_R` in R/perturb.R; ptype `noise90`). Its replicates reuse
  the noise94 random draws (`P_CODE` 3), scaled to the new SEM.
- GMM replicates at k = 2-4 refitted with 10 converged starts (`R/gmm_stab_refit.R`, run after
  08 and before 09); k = 5-6 keep gridsearch. Originals in `output/stability/reps_gridsearch_gmm/`.
- Base solutions from `results_full_refit.rds` (07_stability_spec.R).
- October outputs: `output/archive/stability_archive_2026-10-25/`. Results: README "Stability rerun".
- **Null model** (log scale, age − 9): intercept 2.584, slope −0.058, quadratic 0.003;
  SD of random intercept 0.908 and slope 0.169; residual SD 0.524.

## Open items before building (status 2026-10-24)

- **Reliability: resolved 2026-11-07, r = .90.** Public summaries of the 2001 manual give
  test-retest r of .80-.94 across scales (mean .88), and none ties .94 to Total Problems.
  Jamie judged .94 too high. Achenbach & Rescorla (2000), Table 8-1, gives r = .90 for
  preschool Total Problems; Naar-King, Ellis & Frey (2004, p. 70) give .80-.94 across scales.
  `NOISE_R` in R/perturb.R is now c(noise90 = .90, noise80 = .80).
- **lcmm prediction: done.** latrend's `predictPostprob()` silently returns uniform
  probabilities for lcModelLcmmGMM; `predict_lcmm()` in R/predict_class.R calls
  `lcmm::predictClass` on the underlying hlme model instead (checked in
  archive/code/check_predict_lcmm.R on simulated data).
- **Ward:** nearest centroid, noted in R/predict_class.R.

Original list:

- Verify the CBCL Total Problems test-retest reliability value used for noise level 1.
- Confirm lcmm `predictClass` works on new data through the latrend model object;
  otherwise call `hlme` posterior predictions directly.
- Ward has no native "predict"; nearest centroid is a standard approximation and
  will be noted as one.
