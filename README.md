# ABCD clustering multiverse: CBCL total problems × baseline income

Started 2026-10-23. Asks how much subgroup solutions depend on the clustering
method and analytic choices, and whether the link between poverty and
membership in the highest-symptom class holds across that multiverse.

Data: ABCD 5.1, read in place from
`Dropbox/_Research_In_Progress/abcd-data-release-5.1/core`. Derived files are in
`data_derived/` (Dropbox is approved storage for this project). The project now
lives on the second Mac (`~/Desktop/RStudio/Pitt/`, DUC/MCW-approved; see docs/HANDOFF.md).
Code and docs, with no data, are in the private GitHub repository
`jlhanson5/abcd-clustering-multiverse`. `.gitignore` keeps out `data_derived/`,
`output/` and every data or figure file type. The clusterBMA port is a separate
repository (`jlhanson5/clusterBMA-port`, simulated data only), expected at
`../clusterBMA_port`.

## Run order

```
Rscript 00_install.R          # once
Rscript 01_data.R             # long.rds / wide.rds + output/logs/sample_flow.txt
Rscript 01b_race.R            # race.rds (race/ethnicity for the poverty models) + output/logs/race.txt
Rscript 02_spec.R             # specification table
Rscript 03_run.R pilot 12     # fit all tasks (restartable); "full" for the full sample
Rscript 04_poverty.R pilot    # relabel classes, multinomial poverty model per universe
Rscript 05_agreement.R pilot  # pairwise ARI / NMI
Rscript 06_figures.R pilot    # figures + output/summaries/summary_pilot.txt
Rscript 07_stability_spec.R full ... 10_stability_figures.R full   # stability; see docs/STABILITY_PLAN.md "As built"
Rscript R/gmm_refit.R ...; Rscript R/gmm_refit_merge.R full_refit  # GMM refit, mode full_refit (see "GMM refit")
Rscript 11_refit_figures.R    # extra full_refit figures (fig8-13)
Rscript R/gmm_stab_refit.R 8  # GMM stability replicates with converged starts, before 09
zsh R/run_cbma.sh; zsh R/run_cbma_avg.sh   # clusterBMA consensus (12_cbma.R x 16, then 12b_cbma_poverty.R and 13_cbma_summary.R)
Rscript 14_extra_figures.R    # single methods vs consensus (output/figures/5_single_vs_consensus)
Rscript 15_class_gradient.R   # every class vs the lowest (output/figures/6_class_gradient)
Rscript 16_sens_site_re.R 8   # sensitivity: site as a random effect (output/figures/7_sensitivity)
Rscript 17_permutation.R 500 8; Rscript 17_permutation.R summary   # spec-curve permutation test (~18 h)
Rscript 18_sens_weights.R 1   # sensitivity: ABCD population weights (needs abcd-general/abcd_p_demo.csv under ABCD_REL)
```

`12b_cbma_poverty.R` fits the poverty models (with and without race/ethnicity) on the stored consensus classes without reclustering.

`R/fit_task.R` holds the method fitters.

## Folder layout (reorganised 2026-11-09)

```
*.R, R/               pipeline scripts (numbered, run from the project root) and helpers
docs/                 HANDOFF.md, STABILITY_PLAN.md, DECISIONS_memory.md, FIGURES.md (figure index)
archive/              retired code, incl. the pre-race-adjustment scripts
data_derived/         wide.rds, long.rds, race.rds, weights.rds (participant-level; never leaves this Mac)
output/results/       results_<mode>.rds, agreement_<mode>.rds, spec_*.rds, cbma_poverty.rds, gmm_refit_summary
output/summaries/     summary_*.txt, one per figure script
output/figures/       1_spec_curve (primary run), 2_main_run, 3_stability, 4_clusterbma,
                      5_single_vs_consensus, 6_class_gradient, 7_sensitivity; FIGURES.md links to docs/FIGURES.md
output/permutation/   perm_####.rds, one per permutation (17_permutation.R; perm 0 = observed)
output/logs/          current logs; cbma/ (clusterBMA fits), archive/ (superseded)
output/universes/, stability/, cbma/, gmm_refit/   working files from the fitting steps
output/archive/       superseded outputs, incl. pre_race_2026-11-09/
```

Paths are set once in `00_setup.R` (`RES`, `SUMM`, `fig_dir()`). Result sections below that are dated before 2026-11-09 cite the age/sex/site model; their numbers match `output/archive/pre_race_2026-11-09/`, not the current summaries.

## Locked decisions

| Area | Decision |
|---|---|
| Release | 5.1 now; port to 6.x later |
| Outcome | CBCL total problems, parent report. Raw totals are the main scale; T-scores and log(1+raw) are forks |
| Waves | BL, Y1, Y2, Y3 (Y4 is only partly released in 5.1, n=4,677) |
| Inclusion | ≥3 of 4 waves (main) vs all 4 waves |
| Missing data | Main: all methods cluster the same single seeded pmm imputation. Fork: each method's own handling (lcmm uses all observed rows by likelihood; for the other methods that equals complete cases, i.e. the all-4 sample) |
| Families | Main: one child per family, seeded random pick. Fork: full sample, with family-clustered SEs in the poverty model |
| Time metric | Wave (0-3) vs age − 9 (years). This only changes GBTM and GMM; the other methods treat waves as positions |
| Methods (7) | k-means, PAM (CLARA; exact PAM when n ≤ 2,000), Ward, Gaussian mixture on the 4 wave scores (mclust VVV), kml, GBTM, GMM |
| GBTM / GMM | hlme via latrend, class-specific quadratic in time. GMM adds a random intercept and slope |
| k | 2-6 in the main analysis. The k-max 5 vs 6 fork is reported as a possible finding in itself. The supplement runs k 7-8 for the five cheap methods |
| Selection criteria | BIC (likelihood-based methods only) and silhouette. Every k still goes into the spec curve |
| Class assignment | Modal (hard) class for every method |
| Class order | Relabelled 1..K by mean raw CBCL across BL-Y3; class K is the "highest-symptom class" |
| Poverty exposure | log baseline household income at bin midpoints. The bottom bin (<$5k) is set to $2,500 and the top bin ($200k+) to $250,000 |
| Poverty model | Multinomial logit in every universe, with the lowest class as reference. Covariates: baseline age, sex, site and race/ethnicity (primary since 2026-11-09; ABCD's 5-level `race_ethnicity`, reference White). Age, sex and site only is the sensitivity analysis. The spec curve shows the highest-vs-lowest RRR per halving of income (−log2 income); every contrast is stored |
| Class-size floor | Solutions with any class under 5% are flagged and kept |
| Random starts | 20 per universe (k-means nstart, CLARA samples, kml redrawings, lcmm gridsearch reps) |
| Compute | Local (18 cores); pilot first |

## Exclusions (see output/logs/sample_flow.txt)

- 1,017 children with baseline income refused or answered "don't know" (777/999), or missing, are excluded, not imputed.
- 2 children with sex not coded 1 or 2 are excluded.
- 761 children with fewer than 3 of the 4 CBCL waves are excluded.
- Analysis sample: 10,088 children (9,003 with all 4 waves). 8,331 with one child per family.
- One child has no `interview_age` at any wave. Their ages are filled with the sample median at each wave (after imputation and flags, so nothing else changes). This was found when the full run's age-metric GBTM/GMM tasks failed; it does not affect the pilot.

## Flags: come back to these

- **Race/ethnicity** was added to the poverty models on 2026-11-09 and is now the primary model (see "Race/ethnicity adjustment"). Result sections dated before then report the age, sex and site model. One child with no race code at any wave is left out of the adjusted models only.
- **Random starts (20)** to revisit. Start-to-start stability is not yet recorded per universe; a stability figure needs that.
- **Single imputation.** Full multiple imputation (with class labels aligned across imputations) is deferred.
- **Modal assignment** ignores classification uncertainty in the likelihood-based methods (GBTM, GMM, Gaussian mixture). A three-step correction could be added later as a fork.
- **Gaussian likelihood for raw counts.** GBTM/GMM use a Gaussian hlme on every scale. The censored-normal / link-function version for raw counts is not implemented; the log scale fork partly covers skew.
- **Ward and PAM-at-small-n are deterministic.** Their 20 "starts" don't apply. The Gaussian mixture is deterministic only at n ≤ 2,000 (the pilot). Above that, mclust initialises on a random subset of 2,000 children (`mclust.options("subset")`), so full-sample fits depend on the RNG. `fit_llpa` set no seed in the main run, so those fits are not exactly reproducible: refits with a seed agree with them at ARI .28-.84 (checked 2026-10-24 on the anchor configs). `fit_llpa` now seeds; the stability analysis gives the Gaussian mixture 50 seed replicates. In `full_refit` the 84 Gaussian-mixture universes were refitted with seeds (`R/llpa_reseed.R`, 2026-11-07), so the primary run is reproducible; `results_full.rds` keeps the unseeded fits.
- **Silhouette** is computed on the imputed 4-wave score matrix (on the universe's scale) for every method, including native-missing lcmm fits, so all methods are scored in the same space. A seeded subsample of 3,000 is used when n > 3,000.
- **Site in small pilot classes.** With 21 site dummies and small classes, some multinomial fits separate. The pilot RRRs for tiny classes are unstable; class-gradient points with SE > 5 are dropped in the S1 figure.

## Race/ethnicity adjustment (2026-11-09; current results)

Race/ethnicity is now in every poverty model as the primary analysis: ABCD's 5-level `race_ethnicity` (White 55%, Black 13%, Hispanic 19%, Asian 2%, Other 10%; reference White; `01b_race.R`, `output/logs/race.txt`). The age, sex and site model is the sensitivity analysis (`*_unadj` columns, `fig13_race_sensitivity`, open points in `cbma5_poverty`). No reclustering: the classes are the same, only the poverty models were refit. One child with no race code at any wave is left out of the adjusted models only. The unadjusted reruns reproduce the archived results exactly.

Universes whose highest class has fewer than 10 children now have no RRR (`MIN_TOP_N`, 00_setup.R). This drops 2 GMM k = 6 universes per run (highest classes of 1-2 children), which separate once race is added (one had RRR 47,680). Fits that separate (SE zero, non-finite or above 5; `MAX_SE`) are also not estimable. This removes 6 stability solutions (highest classes of 12-22 children, RRRs of 1e37 and above, with and without race) and no spec-curve or clusterBMA estimate.

- **Primary run (`full_refit`, 652 universes, k 2-6):** median RRR 1.22 (IQR 1.12-1.28), CI above 1 in 92% (601), below 1 in 2. Without race: 1.16 (1.07-1.22), 85%, 19 below 1. Main run (`full`): 1.21, 93%.
- **Adjustment raises the estimate in almost every universe** (adjusted / unadjusted RRR median 1.05, IQR 1.03-1.07, the same on every scale). Race/ethnicity masks part of the income gradient rather than explaining it.
- **Scale still drives the result:** raw 1.37 (CI above 1 in 100%), T-score 1.22 (99%), log(1+raw) 1.10 (78%; 59% without race). The ordering and the conclusions of the sections below are unchanged; the numbers in them are the unadjusted ones.
- **Robustness subsets** (fig8): leaving out GMM k 5-6, unreplicated GMM, all mixture models, or small-class solutions gives medians of 1.20-1.22. By method 1.17 (Gaussian mixture) to 1.24 (GBTM); by k 1.15 (k = 2) to 1.25 (k = 6) (`summary_full_refit_extra.txt`).
- **clusterBMA consensus** (`output/results/cbma_poverty.rds`, from `12b_cbma_poverty.R`): every fit's RRR rose by 0.04-0.10. Raw 1.26-1.38, T-score 1.23-1.25, log 1.09-1.36, pooled all-scales 1.22 (1.17-1.28). The log equal and k-balanced fits (near-zero reference class) are now 1.09 (1.02-1.17) and 1.10 (1.04-1.16), with CIs above 1.
- **Input solutions** (32 per scale): raw median 1.41, T-score 1.21, log 1.11.
- **Full gradient** (`15_class_gradient.R`, `output/figures/6_class_gradient`, k 3-6): on raw scores the gradient rises steadily from the lowest class in 74% of universes, and a middle class falls below the lowest (CI below 1) in 6%. On T-scores and log(1+raw) a middle class is significantly *less* poor than the lowest class in 56% and 75% of universes (any dip: 77% and 87%). The lowest class there is a near-floor group (median 21% of children with mean raw CBCL 3.9 on T-scores, 12% with 2.2 on log; on raw it is 62% with 8.0) whose families are poorer than the mildly symptomatic children. The dip grows with k (CI below 1 in 22% of universes at k = 3, 61% at k = 6). The highest-vs-lowest contrast stays primary; this is the U-shape the clusterBMA log runs showed, and it is general on the transformed scales.
- **Stability replicates** (`summary_stability_full.txt`, `stab4_poverty`; 15,065 solutions, base variant): raw anchor median RRR 1.38, CI above 1 in 100% of 6,091 perturbed solutions (1.32, 100% without race); log anchor 1.12, 90% of 6,026 (1.07, 69%). Null datasets 0.99-1.00, CI above 1 in 1-1.5%. The weak cell is GMM on the log scale under site split-halves and subsamples (median 1.03-1.04, CI above 1 in 23-32%); every other method stays at 1.10 or above.

## Sensitivity analyses (2026-11-09; output/figures/7_sensitivity)

- **Site as a random effect** (`16_sens_site_re.R`, `summary_sens_site_re.txt`, `s1_site_random_effects`): highest vs lowest class as a binary logit (Begg & Gray) with a random intercept for site (glmer), same covariates. All 652 fits converged (12 singular). Median RRR 1.22 vs 1.22 primary, ratio to primary 1.00 on every scale; CI above 1 in 93% (raw 100%, T-score 99%, log 78%). Site SD 0.31 on the logit scale. Full-sample universes have no family term here, so their CIs are slightly too narrow. A first run that also gave them a random intercept for family did not work: most families have one child, 24 fits failed and raw-scale SEs were 16 times the fixed-effects SEs (kept as `*_famRE_attempt`).
- **Population weights** (`18_sens_weights.R`, `summary_sens_weights.txt`, `s3_population_weights`): primary model weighted by ABCD's ACS raked propensity weight (`acs_raked_propensity_score` in `abcd_p_demo`, baseline; `acspsw03` up to release 4.0), rescaled to mean 1, sandwich SEs (family-clustered in the full sample, by child otherwise). Every child in the sample has a weight. Unweighted refits with the same SEs reproduce the 327 full-sample primary estimates exactly. Weighting raises the RRR slightly (weighted / unweighted median 1.02, 10th-90th percentile 1.00-1.03) and widens SEs by about 13%. Median RRR 1.24 vs 1.22; CI above 1 in 93% vs 92% (raw 100%, T-score 99%, log 79%); no universe has a CI below 1.
- **Specification-curve permutation test** (`17_permutation.R`, Simonsohn, Simmons & Nelson 2020): 500 Freedman-Lane permutations of income residuals (on age, sex, site, race), one shuffle shared by all 652 universes, primary model refit in each. Statistics: median RRR, share of universes with CI above 1, Stouffer's Z; p = (b + 1) / (B + 1). Result (2026-11-10): observed median RRR 1.22 vs null 1.00 (95% range 0.96-1.04); CI above 1 in 92% of universes vs 0.2% (95% range 0-19%); Stouffer's Z 160 vs -1 (95% range -33 to 33). No permutation reached the observed value on any statistic, so p = 1/501 = 0.002 for all three (the smallest possible with 500 permutations). Figure: s2_permutation.

## traj: tested, not added (2026-11-10; supplement)

traj (CRAN `traj` 3.1.0; Leffondre et al. 2004 measures + spectral clustering) was tested as an eighth method before adding it to the spec curve and clusterBMA. Scripts: `R/traj_stability_test.R` (plain traj) and `R/traj_variants.R` (fixes); summaries `output/summaries/summary_traj_test.txt` and `summary_traj_variants.txt`. Same anchors (raw, log), the same 20 80% subsamples and 5 null datasets as the stability analysis, k 2-6, ARI on shared children.

- **Plain traj splits off tiny classes.** On raw scores a class of 58 children (0.7%) appears in 81-95% of fits at k 3-6: the flat trajectories (44 children, 38 at 0 every wave), extreme on measure 9 (rate of crossing the fitted line). On log(1+raw), k = 2 on the full data is 8,321 vs 10 children; children with identical score patterns get identical measures, become each other's nearest neighbours and form a cut-off group that spectral clustering splits off first. With the comparison against the lowest-mean class, these groups would become the reference class. Fit seeds do not matter (identical solutions).
- **Median subsample ARI, plain traj:** raw 0.96 / 0.98 / 0.96 / 0.78 / 0.85 at k 2-6 (k = 2 bimodal, lower quartile about 0); log 0.03 / 0.88 / 0.76 / 0.81 / 0.83. k-means and kml: 0.85-0.99.
- **Fixes tried:** dropping measures 8 and 9 (alone or with jitter) fails: tiny classes in 71-100% of fits, because the remaining ties are untouched. traj's own outlier capping does not separate ties either. Jitter (U(-0.05, 0.05) on every score, one draw per child) removes the tiny classes on raw scores (smallest class 9-32%) with subsample ARI 0.99 / 0.98 / 0.87 / 0.94 / 0.70; k = 4 and 6 switch between two solutions. Jitter + capping is the same on raw. On log no fix works: jitter leaves tiny classes in 29-86% of fits, jitter + capping is stable only because every fit splits off the same 0.4% of children.
- **Decision (Jamie, 2026-11-12):** traj stays out of the spec curve, stability analysis and clusterBMA, and is reported in the supplement as tested and excluded. A usable version exists only on raw scores with a non-standard modification (jitter), which would leave traj as the only method not run on every scale.

## Pilot results (2026-10-23, n ≈ 2,000; see output/archive/pilot/summary_pilot.txt)

- 780 universes (660 main with k = 2-6, plus 120 supplement). 13 fits failed (all Gaussian mixture, unconstrained covariance, mostly at k ≥ 7). 74 had an empty class (nearly all GMM). 654 main universes have a poverty estimate.
- Poverty RRR (highest vs lowest class, per halving of income): median 1.13 (IQR 1.06-1.21). The CI is above 1 in 54% of universes and below 1 in 4%. k-max 5 vs 6 makes no difference.
- **Scale drives the result** (39% unique variance in log RRR; method 2.5%, k 3%). Raw scores split off a small severe class (median 8%, mean CBCL 53): RRR 1.22, with the CI above 1 in 95% of universes. log(1+raw) makes the "highest" class a large moderate group (32%, mean 35): RRR 1.05, CI above 1 in 15%.
- Agreement between methods: kml-k-means ARI .92, GBTM-kml .72. The Gaussian mixture and GMM agree with everything else at only .19-.41.
- k selection: silhouette picks k = 2 almost everywhere. GBTM's BIC always picks the largest k allowed, so the k-max fork fully determines GBTM's selected k.
- The class gradient (S1) rises with class mean CBCL for every method except GMM, which is noisy.

## Full-sample results (2026-10-23; see output/summaries/summary_full.txt, output/figures/2_main_run)

The run took about 8.5 hours of wall time on 12 cores in two passes: the age-metric tasks were refit after the missing-age fix, and about 2 hours were lost to lid-closed sleep. Full-sample GMM took about 75-265 min per task.

- 780 universes. The only fit failures are 16 Gaussian-mixture fits that couldn't fit the unconstrained (VVV) model. 75 universes have an empty class (nearly all GMM). 653 main universes have a poverty estimate.
- Poverty RRR (highest vs lowest class, per halving of income): median 1.16 (IQR 1.07-1.21). The CI is above 1 in 85% of universes and below 1 in 4% (23). k-max 5 vs 6: CI above 1 in 83% vs 85%. The magnitudes match the pilot; the pilot was underpowered.
- **Scale is the dominant fork** (49% unique variance in log RRR; method 8%, k 6%, time/family/sample about 0):
  - raw: the highest class is median 8% of children with mean CBCL 53. RRR 1.25, CI above 1 in 100%.
  - T-score: 13%, mean 47. RRR 1.16, CI above 1 in 95%.
  - log(1+raw): 32%, mean 35. RRR 1.05, CI above 1 in 59%.
- **All 23 reversals (CI below 1) are mixture models on log or T scales**: 19 GMM-log, 2 Gaussian-mixture-log and 2 T-scale universes. In them the "highest" class is a majority group (median 56%, mean CBCL 25), so the comparison is no longer severe vs low.
- By method, the share of universes with the CI above 1: PAM 98%, k-means/kml 95%, GBTM 93%, Ward 87%, Gaussian mixture 81%, GMM 65%.
- Agreement: k-means-kml ARI .91, kml-GBTM .74. The Gaussian mixture and GMM agree with everything at .18-.34.
- k selection: silhouette picks k = 2 for every method except some GMM/Gaussian-mixture configs. GBTM's BIC picks the k-max in 69 of 72 configs, so the cap decides it.
- Class gradient (S1): log RRR rises roughly linearly with class mean CBCL for every method. The GMM panel is noisier, and its loess curve is pulled up by one extreme class (mean CBCL ~125).

## Stability results (2026-10-25; see output/archive/stability_archive_2026-10-25/summary_stability_full.txt and figures/)

The anchors are one child per family, ≥3 waves, imputed (n = 8,331), on the raw and log(1+raw) scales, with k = 2-6. There were 2,860 replicates (15,413 solutions); design and changes are in docs/STABILITY_PLAN.md. The run took about 36 hours on 8 cores. Fit failures were Gaussian-mixture VVV fits with no solution (272, mostly at k 5-6 and in the noisy data) and 3 raw-null GMM datasets (all k). 3 of 2,800 site-split prediction directions failed (degenerate raw GBTM/GMM training fits at k 5-6).

Values below are means over k = 2-6, given as log / raw anchor.
- **Seeds:** k-means, kml and GBTM come back identically (ARI ≈ 1.00). PAM (CLARA): .77-.98. The Gaussian mixture falls off with k (raw .97 → .62; log .95 → .39).
- **80% subsamples:** kml .95/.97, k-means .95/.96, GBTM .91/.96, PAM .88/.87, GMM .64/.88, Gaussian mixture .63/.76, Ward .43/.53.
- **Site split-half** (ARI between the test half's own fit and the training half's predicted classes): k-means .83/.88, kml .82/.88, GBTM .76/.88, PAM .78/.79, GMM .39/.80, Gaussian mixture .55/.68, Ward .45/.58.
- **Measurement noise** is the hardest test, especially on the log scale. At r = .80, ARI is .17-.41 on log and .27-.80 on raw. Even at r = .94 it is .25-.60 on log.
- **Settings** matter more than resampling for Ward and the Gaussian mixture. Average or complete linkage vs Ward: ARI .00-.23 (average linkage on raw splits off outliers). Mixture EEE vs VVV: .08-.19; VVI vs VVV: .35-.67. PAM Manhattan vs Euclidean: .74-.78. GBTM linear vs quadratic: .95-.99. GMM variants: .35-.67.
- **Child level** (pooled over seeds, subsamples, site halves and noise): on the raw scale at k = 2, 19% of children are reliably in the highest class with k-means (≥ .90 of replicates) and only 2% are on the boundary. 94% of the main run's highest-class children are reliably high. On the log scale the highest class is about half the sample (k-means 53% reliably high). Ward and GMM have much larger boundary groups (log k = 2: 48% and 32%).
- **Start recording (GBTM/GMM):** GBTM reaches the same best log-likelihood from all 5 converged starts at every k except raw k = 6 (1 of 5), and it reproduces the main-run solution exactly. **GMM does not.** The main run's GMM fits were local optima at log k = 2, 3 and 6 and raw k = 6 (the recorded fit's BIC is lower by 54, 197, 472 and 3,588). At log k = 5 the main run beat all 5 recorded starts by 399. At log k = 2 the better optimum gives RRR 0.93 (CI below 1) where the main run had 1.01. So the GMM-log reversals in the spec curve depend on which local optimum gridsearch reaches.
- **Null benchmark:** the methods do not reject the null.
  - Silhouette picks k = 2 in the null just as in the real data.
  - On the raw scale the null gives *higher* silhouettes than the real data for k-means, kml, Ward, GBTM and GMM. A skewed, floored count with no subgroups splits off a small "severe" class by itself (null highest class: median 3% of children, mean CBCL 116; real: 8%, mean 58).
  - On the log scale real silhouettes are only slightly above the null's.
  - The relative entropy of GBTM on raw is higher in the real data (~.9) than in the null (~.35).
- **Poverty:**
  - Raw anchor: the CI is above 1 in 100% of replicates under every perturbation (median RRR 1.30-1.33).
  - Log anchor: the CI is above 1 in 47% of site halves (half the sample), 65% of subsamples, 85% of seeds and 92-95% under noise (median 1.06-1.10).
  - Null: the RRR is centred on 1.00 and the CI excludes 1 in 3.6% (raw) and 4.8% (log) of null solutions, close to the nominal 5%. So the poverty gradient is not an artefact of the clustering.

## GMM optima experiment (2026-11-03; see output/archive/gmm_optima_experiment/summary_gmm_optima.txt, archive/code/gmm_optima.R)

50 random starts on each GMM anchor (raw and log, n = 8,331) × k = 2-6, each run to the end, with no 10-iteration screen deciding which to converge. The run was paused about 4.5 days by sleep, so raw k = 6 times are inflated.

| k | Log: starts at best, converged | Raw: starts at best, converged | Minutes per start |
|---|---|---|---|
| 2 | 33, 50 | 50, 50 | under 1 |
| 3 | 20, 50 | 3, 50 | about 1 |
| 4 | 17, 50 | 5, 49 | about 3 |
| 5 | 1, 35 | 5, 41 | about 15 |
| 6 | 1, 3 | 1, 8 | about 50 |

- **Gridsearch does not improve with more starts.** The 10-iteration log-likelihood barely predicts the final one (rank correlation -.27 to .23). At log k = 3, gridsearch's chance of reaching the best falls from .40 with 1 start to .17 with 20. The main run's 20 starts reach the best at log k = 2 and 4 only some of the time, and rarely or never at log k = 3 and 5-6 or raw k = 3-6.
- **Converging every start does work where the best can be identified.** For a 95% chance: about 3 starts at log k = 2, 6 at log k = 3, 8 at log k = 4, 49 at raw k = 3 (CI on the last is wide, up to about 240).
- **k = 5-6 and raw k = 4-5 are not reliably estimated.** Each best at k = 5-6 (except raw k = 5) was hit once. Most k = 6 starts hit the 500-iteration cap. Log k = 6's best fits worse than log k = 5's, which cannot be the true maximum. At raw k = 4 and 5 the best converged solution is the 3-class one, and one non-converged run at each reached a log-likelihood about 676 higher, so better solutions exist.
- **The better optima use small classes.** Raw k = 3: 3.1%, 14%, 83%. Raw k = 6 and log k = 6 include classes of 1.2-1.8%.
- **Poverty RRR moves.** Raw k = 3: 1.19 (main) to 1.32 (best). Log k = 2: 1.01 to 0.93 (CI below 1). Log k = 3: 0.92 (CI below 1) to 0.99. On log the "highest" class at these optima is 74-81% of children with mean CBCL 20-22, so the reading that log reversals compare a majority group with a low group holds, but which universes reverse depends on the optimum reached.
- **Main-run GMM (180 universes):** 72 have an empty class and 52 did not converge. The stability analysis used the same gridsearch, so GMM's replicate ARIs mix data sensitivity with optimiser luck. GBTM is not affected: it reproduced the main-run solutions.

## GMM refit (2026-11-05; see output/summaries/summary_full_refit.txt, output/figures/1_spec_curve, R/gmm_refit.R, R/gmm_refit_merge.R)

All 36 GMM tasks at k = 2-4 were refitted with 50 random starts each, every start run to the end (5,400 starts, no errors, 380 CPU-hours, 50 hours on 8 cores). The best converged solution replaces the main-run fit in the `full_refit` mode. k = 5-6 keep the main-run fit and are flagged as not reliably estimated. `results_full.rds` and `output/universes/full/` are unchanged, and all 600 non-GMM universes have identical estimates in the two modes. To rebuild: `R/gmm_refit_merge.R full_refit`, then `04_poverty.R`, `05_agreement.R` and `06_figures.R` with mode `full_refit`.

**Primary run (decided 2026-11-07): `full_refit`.** Report the main run as sensitivity. Write-ups should cite `output/summaries/summary_full_refit.txt` and `output/figures/1_spec_curve/`; the Full-sample results section above describes the main run.

**The multiverse-level result does not change.**

| Main universes (k 2-6) | Main run | Refit |
|---|---|---|
| Median RRR (IQR) | 1.16 (1.07-1.21) | 1.16 (1.07-1.22) |
| CI above 1, of 653 | 553 (85%) | 553 (85%) |
| CI below 1 | 23 | 22 |
| Scale, unique share of variance in log RRR | 49% | 51% |
| Method, unique share | 8% | 6% |
| Empty class | 75 | 58 |
| Smallest class under 5% | 203 | 219 |

Keeping only replicated GMM solutions with no empty class gives the same answer: median 1.16, CI above 1 in 87%, below 1 in 4% (554 universes). GMM as a whole: median RRR 1.16 to 1.17, CI above 1 in 65% both times.

**Individual GMM universes do change.** The refit was never worse by BIC, and it was better by more than 1 point in 6 of 36 universes at k = 2, 31 at k = 3 and 33 at k = 4. Median BIC gains at k = 3-4: raw 1,217 and 1,428, log 124 and 355, T-score 50 and 113. Median ARI with the main-run classes: 1.00 at k = 2, .48 at k = 3, .56 at k = 4. CI category changed in 18 of 108 universes.

| Median RRR, main to refit; universes with CI below 1, main to refit (of 12) | k = 2 | k = 3 | k = 4 |
|---|---|---|---|
| Raw | 1.19 to 1.19; 0 to 0 | 1.19 to 1.30; 0 to 0 | 1.18 to 1.24; 0 to 0 |
| T-score | 1.22 to 1.22; 0 to 0 | 1.17 to 1.21; 0 to 0 | 1.17 to 1.19; 0 to 0 |
| log(1+raw) | 0.97 to 0.93; 6 to 12 | 0.93 to 0.99; 9 to 4 | 0.96 to 0.96; 4 to 2 |

- **Log k = 2 now reverses in every configuration.** The 6 that were not below 1 in the main run were all local optima (BIC gains 54-136). One, `gmm_log_wave_full_all4_k2`, went from CI above 1 to CI below 1 (1.05 to 0.92), the only universe to cross from one side to the other.
- **The better optima use smaller classes.** Universes with a class under 5%: k = 3 from 16 to 27 of 36, k = 4 from 30 to 35. Empty classes at k = 3 fall from 11 to 0, and at k = 4 from 17 to 11. Non-converged fits at k = 4 fall from 6 to 0.
- **BIC among k = 2-4, excluding empty-class solutions:** log and T-score pick k = 4 in all 12 configurations; raw picks k = 3 in 11 and k = 4 in 1.

**How well determined the solutions are.** `gmm_flag` marks each GMM universe. "replicated" means at least 2 starts reached the best converged log-likelihood within 0.5, and no non-converged start in the refit beat it.

| k | Replicated | Best found once | Median starts at best |
|---|---|---|---|
| 2 | 36 | 0 | 50 |
| 3 | 30 | 6 (5 raw by wave) | 4 |
| 4 | 19 | 17 (raw 8, T-score 6, log 3) | 2 |

Experiment anchors: log k = 2-4 and raw k = 2-3 reach the same best log-likelihoods in fresh starts.

- **Raw k = 4 remains not reliably estimated.** 11 of 12 best solutions have an empty class. At the anchor the refit's best (-128764, found once) beats the experiment's best converged (-128958), but one non-converged experiment run reached -128282. So a better solution exists that the flag cannot see, because the flag only compares starts within the refit. Treat every raw k = 4 universe as flagged, whatever `gmm_flag` says.
- **T-score by wave at k = 4 mostly does not converge.** Converged starts of 50 at k = 4:

| | By age (of 300) | By wave (of 300) |
|---|---|---|
| Raw | 294 | 287 |
| T-score | 242 | 27 |
| Log | 300 | 300 |

  Failure takes the combination: neither T-scores nor wave alone does it. Most T-score-by-wave starts hit the 500-iteration cap; the cause was not investigated. Each such universe still has at least 2 converged starts, and 5 of the 6 are "best found once".

The stability analysis was rerun against `full_refit` with converged GMM starts on 2026-11-08 (next section). The Full-sample results section above describes the main run.

After the Gaussian-mixture reseed (2026-11-07, `R/llpa_reseed.R`) the `full_refit` multiverse has median RRR 1.16, CI above 1 in 555 of 654 universes (85%) and below 1 in 19. The GMM-only refit summary is kept in `output/archive/`.

## Stability rerun (2026-11-08; see output/summaries/summary_stability_full.txt, output/figures/3_stability)

Three changes from the October run, which is kept in `output/archive/stability_archive_2026-10-25/`:
- **Noise level 1 is r = .90, not .94** (decided 2026-11-07 from Achenbach & Rescorla 2000, Table 8-1, preschool Total Problems r = .90, and Naar-King, Ellis & Frey 2004, p. 70, .80-.94 across scales). SEM is 5.5-5.7 raw points. The noise90 replicates reuse the noise94 random draws, scaled to the new SEM.
- **GMM replicates at k = 2-4 use 10 converged starts** (`R/gmm_stab_refit.R`, `fit_gmm_conv` in R/fit_task.R), not gridsearch. k = 5-6 keep gridsearch. The originals are in `output/stability/reps_gridsearch_gmm/`.
- **Base solutions come from `full_refit`.**

2,860 replicates, 15,387 solutions, 313 fit errors (October 287). All 26 extra errors are Gaussian-mixture VVV fits on the noisier r = .90 data (101 vs 75). The 2 raw-null GMM datasets that hlme cannot fit fail at every k, as before.

- **Poverty: unchanged.** Raw anchor: CI above 1 in 100% of replicates under every perturbation (medians 1.19-1.46). Log anchor: medians 1.03-1.14. Null: medians 1.00, CI above 1 in at most 6%.
- **Noise at .90 vs .94:** replicates agree a little less (log k-means k = 2 ARI .51 vs .60; raw .82 vs .86). The poverty RRR barely moves (raw k-means 1.39 vs 1.40).
- **GMM with converged starts is more stable on log, unchanged on raw.** Log subsample ARI at k = 2/3/4: .88/.77/.76 (October .57/.70/.67). Log site split-half ARI: .23/.65/.72 (October .16/.43/.46). On log, GMM still shows no poverty effect across replicates (median RRR about 1.00, CI above 1 in 6-11%), as in the refit.
- Other methods are essentially as in October (their replicates were not refitted).

## clusterBMA consensus (2026-11-08; see output/summaries/summary_cbma.txt, output/figures/4_clusterbma)

`12_cbma.R` combines the anchor configuration's solutions (one child per family, ≥3 waves, imputed, n = 8,331) across methods and k with the native R port of clusterBMA (Forbes et al., 2023; `../clusterBMA_port`). Hard (modal) inputs, K_BMA = 6 (empty clusters drop out), 5,000 Adam steps, seed 1, vecLib BLAS. About 85 minutes per fit.

Figures (`13_cbma_summary.R`): cbma1 weights, cbma2 class trajectories, cbma3 uncertainty, cbma4 input-consensus ARI, cbma5 poverty RRR per fit, cbma6 ARI between weightings, cbma7 uncertainty by class. cbma7 shows that the small classes behind the log scale's outlying RRRs are the least certain: the S_Dbw 4% high class (mean uncertainty .50) and the equal-weight 7% near-zero class (.46), against about .20-.30 for the other log classes. On raw the large lowest class is the most certain (.11-.16) and the two higher classes sit at .26-.30.

- **Inputs:** k = 2-6 for all 7 methods, less errors, degenerate fits, GMM k = 5-6 and raw GMM k = 4: 32 per scale. The leave-outs are uneven by k (7 inputs at k = 2-3, 5-6 at k = 5-6) and all model-based (GMM; Gaussian mixture at k = 6 on T-score and log).
- **Weightings:** Calinski-Harabasz (the paper's default), S_Dbw, equal, averaged (mean of the three; the consensus is linear in the weights, so this equals averaging the three consensus matrices) and k-balanced (each k 1/5 of the weight, split equally within k; added 2026-11-08 because of the uneven inputs). Plus one pooled run: all 96 inputs from the three scales, equal weights.

| Scale | Clusters | High class (% of children) | RRR, high vs lowest class |
|---|---|---|---|
| Raw (5 weightings) | 3 | 8-15% | 1.22-1.33 |
| T-score (5) | 3 | 16-20% | 1.18-1.20 |
| log(1+raw) (5) | 3, or 4 with equal and k-balanced | 4-27% | 1.02-1.26 |
| All scales pooled, equal | 3 | 24% | 1.16 (1.12-1.21) |

- **Three flat trajectory classes** (low, moderate, high) under every weighting on raw and T-score; 6 clusters were never needed.
- **Weighting matters little on raw and T-score** (ARI between weightings .82-.95). S_Dbw is the outlier: a smaller, more extreme high class and a larger RRR. Calinski-Harabasz gives the T-score Gaussian mixture at k = 2 a weight of 0.00005 (its split is not variance-based).
- **Log is fragile** (ARI between weightings .46-.81). S_Dbw makes a 4% high class. Equal and k-balanced weights split off a near-zero class (7% and 12% of children, mean raw CBCL about 1-3 at every wave). With that class as the reference, the RRR is 1.02-1.04. Against the ordinary low class (equal weights), the near-zero class has RRR 1.12 (1.06-1.19) and the high class 1.14 (1.10-1.19): poverty is U-shaped at the bottom of the scale. For the write-up, compare with a fixed reference class rather than whichever class has the lowest mean.
- **GMM carries 4-14% of the total weight** in every run (fewer inputs, low index values). The consensus is driven by k-means, PAM, kml, GBTM and Ward, and four of those are distance-based on the same wave scores.
- **Allocation uncertainty** (1 − probability of the assigned class) averages .17-.20 on raw and .22-.27 on T-score and log. With hard inputs it measures disagreement between methods and k.
- **Convergence:** the objective still changed by about 1 part in 10,000 over the last 500 steps. A 10,000-step rerun of log/equal (2026-11-08, 162 min) lowered the objective by a further 0.03% and barely moved the labels: ARI .96 with the 5,000-step fit, the same 4 clusters (sizes 9/27/38/27 vs 7/28/38/27; 136 children moved from class 2 into the near-zero class 1), the same mean uncertainty (.263) and poverty RRR 1.03 (0.98-1.09) vs 1.02 (0.96-1.08). 5,000 steps is enough for these results. Files: `output/archive/smoke/cbma_smoke/log_equal.rds`, `output/logs/cbma/cbma_conv_check.log`.
- **Caveats:** only the port's Calinski-Harabasz path is compared with the original package. S_Dbw (and its minimise direction) is unit-tested only, and returns non-finite values on raw at small n (it is finite at n = 8,331). Equal, averaged and k-balanced weights are not in the original package; they skip the index and use the same consensus and factorisation code, which is parity-tested.

## Single methods vs consensus (2026-11-09; see output/summaries/summary_extra_figures.txt, output/figures/5_single_vs_consensus)

`Rscript 14_extra_figures.R` (after 13). It puts the 32 input solutions per scale (full_refit, anchor configuration, k 2-6) next to the 16 consensus fits. The poverty model is refit on every consensus class and reproduces the stored highest-vs-lowest RRRs exactly.
- **x1** poverty RRR of every input and consensus fit. Raw inputs spread 1.15-1.51 (median 1.34) and rise with k, because the highest class gets smaller and more severe. The consensus sits at 1.22-1.33, below most k 4-6 inputs. T-score inputs: 1.05-1.21 (median 1.14), consensus 1.18-1.20. On log the inputs are 0.93-1.19 (median 1.06) and the consensus spans almost the same range (1.02-1.26).
- **x2** size and mean CBCL of every class. The highest class runs from 2% (mean CBCL ~80) to 45% (~30) on raw, and up to 80% (~20) on log.
- **x3** crosswalk of the equal-weight consensus against each method's k = 3 solution.
- **x4** child-level uncertainty against mean CBCL: uncertainty peaks at the class borders (around CBCL 15-20 and 30-35 on raw).
- **x5** RRR of every class vs the lowest, against class mean CBCL. **U-shape:** in the log equal-weight and k-balanced fits the reference class is a near-zero group (7-12%), and the next two classes have RRR *below* 1 against it (0.88-0.90). The near-zero children come from poorer families than the mildly symptomatic ones, so the highest-vs-lowest contrast on these fits understates the gradient. See "Still open" in HANDOFF (fixed reference class).
- **x6** share of input solutions that put each child in their highest class. With equal weights, 6% (raw), 14% (T), 27% (log) and 14% (pooled) of children are in the highest class in at least half of the inputs; only 1.4% (raw) are in it in every input.
- **x7** ARI between all 16 consensus fits. Agreement is high within a scale (except log S_Dbw) and low across scales (raw vs log .16-.38). The pooled all-scales fit is closest to log CH (.83) and T-score (.59-.74).

## Package notes

- kml 2.5.0, longitudinalData 2.4.7 and clv 0.3-2.5 were archived from CRAN on 2026-01-16. `00_install.R` installs them from the archive with two minimal patches:
  - clv: deletes an unused `#include <R_ext/PrtUtil.h>`, which R 4.6 removed.
  - longitudinalData: drops the misc3d dependency (3-D plots only; misc3d needs XQuartz).
- rgl is built from source with `--disable-opengl` because this Mac has no XQuartz.
- latrend's GBTM wrapper ignores gridsearch initialisation. GBTM is therefore fitted through `lcMethodLcmmGMM(random = ~ -1)`, the same model, so it gets the same 20-start gridsearch.
- Exact versions: `output/logs/package_versions.txt`; session info: `output/logs/sessionInfo_<mode>.txt`.
- `renv.lock` (2026-11-12, R 4.6.1, 160 packages) records every package the scripts use at the installed versions, matching `package_versions.txt`. The project is not renv-activated; the scripts use the normal library. To rebuild elsewhere: `renv::restore(lockfile = "renv.lock")`, then `00_install.R` for kml, longitudinalData, clv (archived from CRAN, patched) and rgl (built without OpenGL).
