---
name: abcd-clustering-multiverse
description: ABCD 5.1 CBCL clustering multiverse x poverty (latrend engine); design decisions locked 2026-10-23 before any code
metadata:
  node_type: memory
  type: project
  originSessionId: 600bafe0-be4b-4019-8374-7488a95e72e3
  modified: 2026-10-23T16:25:26.804Z
---

Started 2026-10-23 from Jamie's pasted spec. The planned folder is Dropbox/RStudio/Pitt/ABCD_clustering_multiverse_20261023 (not created yet). Same 5.1 release and dictionary paths as [[abcd-asthma-cognition]].

Decisions Jamie locked on 2026-10-23:
- Release: 5.1 now, port to 6.x later.
- Wave window: BL-Y3 only (4 annual CBCL waves), because Y4 is partial in 5.1 (n=4,677).
- Families: main analysis keeps one random child per family (seeded); a fork uses the full sample with family-clustered SEs.
- Poverty target class: highest mean across waves; most-worsening (largest slope) class as a secondary fork.
- Engine: latrend 1.6.3. kml was archived from CRAN on 2026-01-16 because 'clv' was archived; Jamie chose to install kml and clv from the CRAN archive and pin them with renv.
- Missing data: main analysis gives every method the same multiply imputed data; a fork uses each method's native handling; the spec table marks invalid combinations.
- Scale: raw CBCL total problems is the main scale (skew 1.9); T-scores and log1p are forks.
- Poverty exposure: log of the baseline income bin midpoint (demo_comb_income_v2). Keep it simple, with no income-to-needs, ADI or hardship for now. The edge bins are $2,500 and $250,000. Exclude income refuse/don't know (777/999) and report the exclusion.
- Poverty-model covariates: baseline age, sex and site only. Flag in the write-up that there is no race/ethnicity adjustment; it may be added later.
- Spec curve: every k goes in (a variation on option A), not only the k a criterion selects. The k range is 2-6, with a k-max fork of 5 vs 6, which Jamie flagged as possibly its own finding. Supplement: k=7-8 for the cheap methods.
- Methods (7): k-means, PAM, Ward, mclust LLPA, kml, GBTM, GMM. HDBSCAN and two-step LMKM come later.
- Class-size floor: flag any solution with a class under 5%.
- Starts: 20 per universe. Jamie flagged this to revisit later.
- Storage: Dropbox is approved. No git.
- Jamie wants figures for the diagnostics too: convergence failures, class sizes, k range.
- Poverty model: a multinomial in every universe, with the lowest-mean CBCL class as reference. The main spec curve shows the highest-vs-lowest RRR per halving of income (from -log2 income), and every contrast is stored. Supplementary figure: each class's RRR plotted against its class mean CBCL, across all universes. Jamie found a binary "target vs all" OR confusing; income stays continuous.
- Time: age is the main metric, with wave as a fork.
- Compute: everything runs locally, accepting that it is slower.

5.1 facts checked on 2026-10-23: 9,807 families, with 4,050 children in multi-child families. BL income is binned (demo_comb_income_v2, 1,015 answered refuse or don't know). demo_roster_v2 contains 0 and 77 codes. ADI (led_l_adi) is baseline only. BPM in 5.1 is youth and teacher report only, with no parent BPM, and youth BPM starts at 6 months. Ages at a given wave span about 2 years.

**How to apply:** keep these defaults unless Jamie reopens them. The defaults proposed but not yet confirmed (age vs wave, INR method, compute) are in the session's closing message and the README once written.

2026-10-23 build: the project folder exists (00_install, 01_data, 02_spec, 03_run [restartable, atomic saves], 04_poverty, 05_agreement, 06_figures, R/fit_task.R, README with all flags). The pilot (n≈2,000; 132 tasks / 780 universes) is done. A Mac crash mid-run was recovered by restarting with caffeinate. Pilot result: median RRR 1.13. Scale is the dominant fork: raw isolates a small severe class (RRR 1.22) while log(1+raw) makes the "highest" class a large moderate one (RRR 1.05). Method and k matter little for the poverty effect but a lot for agreement (Gaussian mixture and GMM disagree with the others). GBTM's BIC always picks the k-max. The full-sample run has not been started.
2026-10-23 full run done (results_full.rds, figures/full, README "Full-sample results"). Median RRR 1.16, CI above 1 in 85% of universes. Scale explains 49% of the variance in log RRR (raw 1.25 vs log 1.05). All 23 reversals are mixture models on log/T, where the "highest" class is a majority group. Gotchas for reruns: one child has no interview_age (filled with wave medians in 01_data.R); full-sample GMM takes ~1.5-4.5 h per task; lid-closed sleep pauses runs even under caffeinate.
2026-10-24 stability interview (STABILITY_PLAN.md in the project). Goals: are the classes real, optimization reliability, and a null benchmark (Jamie did not pick "honest CIs", but the poverty RRR is tracked as a secondary outcome). Design:
- Anchors: raw + log, one per family, ge3_imp, age time metric, 7 methods, k 2-6.
- Perturbations: seeds, 80% subsamples, repeated site split-half (prediction strength), settings (linkage, PAM distance, mixture covariance, trajectory shape; stability only, not new spec-curve forks), noise at the published CBCL reliability and at .80.
- Tiered reps: cheap methods 100 subsamples / 50 seeds / 20 splits / 20 null; GBTM/GMM 20 / recorded starts / 20 / 10.
- Null: a 1-class growth model fit on log, back-transformed for raw.
- Headline metrics: solution-level ARI and child-level consistency.
- Placement: main paper. Next step: Jamie reviews the plan, then build and time it.
2026-10-24 (second Mac, M5 10 cores, project on ~/Desktop/RStudio; Jamie confirmed DUC/MCW-approved): built 07-10 + R/perturb.R, R/predict_class.R, fit_task.R extensions (settings args, lcmm_record). Smoke-tested. Findings: latrend predictPostprob silently returns uniform probs for lcmm models (use predict_lcmm); mclust is RNG-dependent at n > 2,000 (subset init), so main-run full-sample Gaussian-mixture fits were not reproducible; fit_llpa now seeds and llpa gets 50 seed reps (2,860 reps). Start recording converges top 5 of 20 starts. 22 sites, not 21. Reliability .94 still unverified.
2026-10-25 full stability run done (2,860 reps, ~36 h on 8 cores). Headlines: k-means/kml/GBTM stable (subsample ARI .91-.97, site .76-.88); Ward, Gaussian mixture and GMM much less so; noise at r=.80 is the hardest test (log scale ARI .17-.41). Settings change Ward/mixture solutions a lot. GMM main-run fits were local optima at several k (log k=2 better optimum gives RRR 0.93 vs 1.01). Null: methods find classes in no-subgroup data (raw null gives higher silhouettes than real and its own small severe class). Poverty: raw CI>1 in 100% of replicates; null false-positive 3.6-4.8% (nominal). Results in README "Stability results".

## 2026-11-05: GMM refit
Jamie chose to refit the 36 GMM tasks at k = 2-4, 50 converged starts each, keeping all 50 for T-score k = 4 despite the cost. Merged as mode full_refit; main-run files kept for comparison. Results in README "GMM refit".

2026-11-07: Jamie chose full_refit as the primary run for the write-up; the main run (results_full.rds) is reported as sensitivity.
2026-11-07: Jamie judged the .94 CBCL reliability too high and supplied sources; noise level 1 is now r = .90. Jamie asked for: a seeded rerun of the Gaussian-mixture fits (done, R/llpa_reseed.R), GMM stability replicates refitted with converged starts (done), more full_refit figures (11_refit_figures.R), and clusterBMA consensus across methods with three weightings (Calinski-Harabasz, S_Dbw, equal), plus averaged weights per scale and a pooled all-scales equal-weight run. Code pushed to private GitHub repos (no data).
2026-11-08: Jamie added k-balanced weighting after seeing the inputs are uneven by k (model-based methods missing at k 5-6). Jamie went ahead without adding clustering methods (Claude's recommendation for this paper); a skew-aware trajectory model (lcmm with a link function) is the first candidate if a reviewer asks. Stability rerun and clusterBMA results are in the README.
2026-11-09: Jamie made race/ethnicity adjustment the primary poverty model (ABCD race_ethnicity, 5 levels, reference White); age/sex/site only is now the sensitivity analysis. Poverty models refit for the spec curve, clusterBMA (on stored consensus classes, 12b_cbma_poverty.R) and stability, with no reclustering. Jamie chose to treat the RRR as not estimable when the highest class has fewer than 10 children (2 GMM k = 6 universes per run separated once race was added). Adjusted median RRR 1.22 (CI above 1 in 92%) vs 1.16 (85%) unadjusted. Pre-race outputs and code archived (output/archive/pre_race_2026-11-09/, archive/code/pre_race_2026-11-09/). Jamie asked for a full folder restructure: output/results, output/summaries, figure folders named by purpose (1_spec_curve ... 5_single_vs_consensus), docs/ for planning files, docs/FIGURES.md as the figure index.

2026-11-09 (later): Separation rule extended: a fit with SE zero, non-finite or above 5 (MAX_SE) is not estimable, after 6 stability solutions gave RRRs of 1e37+; no spec-curve or clusterBMA estimate affected. Jamie asked for three sensitivity/alternative analyses, none replacing the primary model: (1) Simonsohn-style specification-curve permutation test, 500 permutations chosen by Jamie (17_permutation.R, Freedman-Lane on income residuals, one shuffle across all universes); (2) site as a random effect (16_sens_site_re.R; site-only random intercept after a family random intercept proved unidentifiable in the full sample); (3) ABCD population weights (18_sens_weights.R; acs_raked_propensity_score from abcd_p_demo, baseline, the 5.x name of acspsw03). Children with missing income stay excluded (Jamie). The clusterBMA port goes on a public GitHub repo (Jamie). Race-framing wording, renv and main/supplement figure labels are on hold.

2026-11-12: Jamie asked to test traj (CRAN) before adding it. Plain traj splits off tiny classes (flat trajectories on raw, tied patterns on log); of four fixes only jitter on raw scores works, and k = 4 and 6 are bimodal there. Jamie chose supplement only: traj is reported as tested and excluded, not added to the spec curve, stability analysis or clusterBMA (README "traj: tested, not added"). The permutation test finished (p = 0.002 on all three statistics; README).

2026-11-12 (later): Jamie's answers to the open questions. Race/ethnicity methods wording and the U-shape text drafted in docs/draft_text.md (U-shape is a caveat, not a main finding). GMM k 5-6 stays in the spec curve. Figure main/supplement split: Jamie still thinking. Primary model (race-adjusted, fixed site) with site RE, population weights and age/sex/site-only as sensitivity: agreed. renv.lock added (not activated). Zenodo DOI: later. clusterBMA port history rewritten to the GitHub no-reply email (force-pushed). Main repo is already public (code and docs only, checked). No income imputation. Skew-aware trajectory model (lcmm with an estimated spline link): Jamie asked to do it; feasibility pilot first (R/lcmm_link_pilot.R), then a decision on adding it.
