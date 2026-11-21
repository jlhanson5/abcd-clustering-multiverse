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

## 2026-11-13: skew-aware trajectory model (lcmm, estimated link)

- Pilot (`R/lcmm_link_pilot.R`, full anchor data, 10 starts): the GBTM-style model (no random effects, `5-quant-splines` link) converged to the same solution from every converged start at k = 2-5 (classes >= 6%); k = 6 converged from 1 of 10 starts and added a 1.5% class. The GMM version (random intercept and slope) did not work: k = 3 and k = 4 converged from no start, most starts fell back to the 1-class likelihood, and BIC barely improved on one class (237,245 at k = 1, 237,078 at k = 4). GMM with a link goes in the supplement as tested and not feasible.
- Stability test (`R/lcmm_link_stability.R`, 20 subsamples + 5 nulls, 3 starts): median pairwise subsample ARI 0.975 / 0.958 / 0.945 / 0.961 / 0.916 at k = 2-6 (vs full-data fit 0.98 / 0.97 / 0.97 / 0.97 / 0.95); no class under 1% in any subsample fit; at k = 6 no start converged in 6 of 20 subsamples. Comparable to plain GBTM and k-means on the raw anchor (more stable than plain GBTM at k = 5: 0.96 vs 0.89).
- Null data: smallest class shrinks quickly (median 3.6% / 1.3% / 0.8% at k = 4-6, vs 9.4% / 6.2% / 1.7% in the real subsamples), but BIC falls with k in null data too. Pairwise ARI between null replicates is 0, as expected for independent datasets, so it is not informative.
- Decision (Jamie, 2026-11-13): add it as an eighth method, `gbtm_link`, on raw scores only, to the spec curve, the stability analysis (raw anchor) and clusterBMA (raw and pooled). GMM with a link goes in the supplement as tested and not feasible.
- Integration (2026-11-13/14, `R/run_link_integration.sh`; pre-link results in `output/archive/pre_link_2026-11-13/`): 12 tasks, 60 universes, no errors. Primary run now 712 universes: median RRR 1.22, CI above 1 in 93% (661); gbtm_link alone 1.26. Stability: subsample ARI .98/.96/.95/.96/.87, site prediction strength .89/.86/.76/.71/.42, but noise-90 ARI only .74/.61/.52/.48/.42 (plain GBTM .85-.73). Raw clusterBMA: 37 inputs, gbtm_link 13% of the weight, consensus ARI .93-.98 with the pre-link fits, RRRs within 0.01. Permutation p = 0.002 on all three statistics. Population weights and site random effect unchanged.
- The rerun clusterBMA fits used R's reference BLAS (libRblas.0.dylib; 200-230 min per fit) where the earlier ones used vecLib (about 85 min). T-score and log fits were not rerun.
- 11_refit_figures.R had 7 methods hard-coded in Fig 8; fixed to `length(METHODS)`.

## 2026-11-14: growth-curve k-means (GCKM) pilot

Jamie asked which other methods to consider (latrend). With 4 waves, only GCKM (`lcMethodGCKM`: lmer random intercept and slope, then k-means on each child's standardised coefficients) adds a different family; DTW, functional, ZIP/count GBTM, anchored k-medoids, NPRM and MixAK were judged redundant, unsuited to 4 waves or too slow. Pilot `R/gckm_pilot.R` (latrend defaults written out, 20 k-means starts; linear and quadratic random effects; raw and log anchors; full data, 20 subsamples, 5 nulls), summary `output/summaries/summary_gckm_pilot.txt`. Runs in about a minute on 4 cores.
- Linear (latrend default): median pairwise subsample ARI raw .99/.98/.97/.97/.94, log .87/.97/.96/.77/.88 at k 2-6; passes the proposed clusterBMA rule (>= .75 at every k on both anchors), narrowly at log k = 5. Quadratic is less stable (raw k = 6 .74).
- Not redundant: ARI with k-means, kml, PAM and GBTM base solutions .26-.57. Standardising gives the slope as much weight as the level, so GCKM groups children by change as well as severity. Unstandardised coefficients would be dominated by the intercept and close to k-means.
- Null data: on raw the smallest null class shrinks to 0.1-1% at k 3-6 (real 3.6-13%); on log the null classes are as balanced as the real ones (GCKM partitions a continuous cloud either way). Null pairwise ARI is 0 and uninformative (independent null datasets).
- Site split-half and poverty (added the same day): held-out children get BLUP coefficients from the training half's fixed effects and variance components (checked: reproduces coef(m) to 1e-13 and the fitted classes exactly). Linear, age time: site prediction strength raw .90/.83/.76/.66/.58, log .72/.80/.69/.55/.55 (k-means raw .91/.87/.77/.70/.51, log .92/.87/.85/.60/.57). Poverty RRR, full data: raw 1.23/1.28/1.41/1.24/1.24, log 1.14-1.19; below k-means at k 5-6 on raw (1.55) because GCKM's highest-mean class is less severe (7%, mean CBCL 47-49).
- **Visit timing problem with age as time.** In 20 raw null datasets the poverty CI was above 1 in 10-12% of solutions at k 3-6 (k-means 1%, kml 2.5%, GBTM 0% in the stability nulls). Follow-up span is slightly longer in poorer families (partial r = .05 after baseline age and site), and the standardised slope picks it up. With wave number (0-3) as time ("linwave") the rate is 0-6% (pooled about 1%), and subsample stability is better on log (.98/.98/.95/.92/.86). Wave-time site prediction strength raw .90/.81/.76/.51/.66, log .89/.81/.60/.69/.53; RRR raw 1.24/1.28/1.41/1.32/1.29, log 1.02-1.21. If GCKM is added, use wave time.
- Open: whether to add it (Jamie to decide). Discussed: it groups children by change as well as level, so it answers a somewhat different question from the other methods.

## 2026-11-15: clusterBMA inputs finalised

Decided (Jamie, 2026-11-15). Replaces the all-methods input set, which becomes a sensitivity analysis.
- **Rule** (stability and estimability only, set without reference to the consensus results): a method goes in if it fits on all three scales, is estimable at every k, and has median pairwise subsample ARI >= .75 at every k on both stability anchors (raw, log).
- **Methods:** k-means, kml, PAM, GBTM (age time, quadratic). Out: Ward (subsample ARI .40-.58; linkage-dependent), Gaussian mixture (falls to .36-.55 above k = 3; VVV failures), GMM (not estimable at k = 5-6), GBTM (link) (raw only; noise-sensitive; stays in the spec curve), traj (excluded 2026-11-12). GCKM not added.
- **k = 2-8.** k-means, kml and PAM already have k = 7-8 (supplement universes); GBTM needs new fits at k = 7-8. The rule must also hold at k = 7-8, so all four need stability replicates at k = 7-8 before the consensus runs.
- **Variations:** each method's one-per-family fit and its all-children fit (siblings included), both on the anchor configuration (>= 3 waves, imputed), scored on the 8,331 anchor children. Not used: PAM Manhattan (fails the rule: raw k = 6 .69, log k = 5 .74), GBTM linear (near-duplicate of quadratic, ARI .95-.99), GBTM wave time and observed waves (would give GBTM four times the inputs of the other methods).
- **Inputs:** 4 methods x 7 k x 2 = 56 per scale, 168 pooled; every method and k equally represented, so k-balanced weighting is no longer needed. Weightings: Calinski-Harabasz, S_Dbw, equal, averaged; pooled all-scales run with equal weights.
- **To run:** GBTM k = 7-8 fits (3 scales x 2 family settings), stability at k = 7-8, then all clusterBMA fits, 12b-14 and the class-gradient figures. Not started; the machine is busy with the cognition / positive affect project.
- GCKM (decided by Jamie, 2026-11-15): left out of the spec curve as well as clusterBMA; reported in the supplement as piloted and not added (like traj). The clusterBMA reruns start after the cognition / positive affect project's current jobs finish.

## 2026-11-18: running the 2026-11-15 input set
- **GBTM k = 7-8 fitted** (R/gbtm_k78.R: 6 tasks x 2 k, task seeds as 03_run.R gives them, about 23-69 min per fit). Merged into the full_refit task files (originals moved to output/archive/pre_k78_2026-11-18/). 02_spec.R gives these six tasks ks = 2:8, so they are supplement universes like the k-means, kml and PAM k = 7-8 fits; then 04_poverty.R full_refit was rerun. Results at k = 2-6 are unchanged; 05, 06 and 11 were not rerun (follow-up if the spec curve should show these).
- **Degenerate fits:** GBTM raw, all children, ends up with empty classes at k = 7 (5 classes found) and k = 8 (7 found). Every other fit (log, t, one per family) is non-degenerate, smallest class 1.2-5.5%. 12_cbma.R leaves degenerate fits out (they are listed in `left_out`), which would give raw 54 inputs and the pooled run 166. Going ahead on that basis (left out, as both projects already handle degenerate fits; raw 54 inputs, pooled 166): asked Jamie 2026-11-18, no answer yet, revisit if they want otherwise. The rule's "estimable at every k" is met by the one-per-family fits on all scales.
- **Stability at k = 7-8:** 08_stability_run.R takes a `k=7|8` filter key (methods kmeans, kml, pam and gbtm only; same REPS design and replicate seeds as k = 2-6) and writes to output/stability/reps_k78/ (log stability_full_k78.log), so the k = 2-6 replicates are untouched. Subsample replicates run first. R/cbma_rule_k28.R checks the .75 rule at k = 2-8 (k = 2-6 from summary.rds, k = 7-8 from reps_k78), and `09_stability_summary.R full <cores> k78` summarises the full k = 7-8 set (summary_k78.rds).
- **clusterBMA:** 12_cbma.R takes `CBMA_SET` = rule (default; output/cbma, K_BMA = 8) or allmethods (output/cbma_allmethods, the old set, K_BMA = 6). R/run_cbma_k28.sh runs the 13 fits two at a time, then 12b, 13, 14 and 15. 04_poverty.R takes POV_CORES (default 8).

## 2026-11-18: main and supplement figures
Decided by Jamie on 2026-11-18 ("Yes, log both", relayed by the main session). There are six main figures:
1. fig1_spec_curve
2. fig7_decision_importance
3. fig10_top_class
4. stab1_heatmap and stab3_null as one figure
5. cbma2_trajectories and cbma5_poverty as one figure
6. g1_gradient_lines

Every other figure goes in the supplement. That includes fig13_race_sensitivity, the 7_sensitivity figures s1-s3 (the permutation statistic is also reported in the text), single-vs-consensus x1-x7 and the diagnostics. The two combined figures are drawn at the redraw after the k = 2-8 clusterBMA fits. docs/FIGURES.md has a Role column to match.
- Zenodo stays "at submission": Jamie links Zenodo to GitHub, then a release is tagged. Nothing to do before then.
- Jamie also confirmed (in the cognition project) that the rule set is primary and the all-methods set is the sensitivity analysis. The .75 rule is applied as written, with no readmitting of near misses.

## 2026-11-18: .75 rule at k = 2-8, stopped before the clusterBMA fits
R/cbma_rule_k28.R uses the k = 7-8 subsample replicates in output/stability/reps_k78/: 100 per cell for k-means, kml and PAM and 20 for GBTM, with the same design as k = 2-6. Output is in output/summaries/summary_cbma_rule_k28.txt.

Median pairwise subsample ARI at k = 7 and k = 8 (log anchor, then raw):

| Method | log k7 | log k8 | raw k7 | raw k8 | Result |
|---|---|---|---|---|---|
| k-means | .897 | .772 | .902 | .950 | Pass; lowest is log k8 .772 |
| kml | .875 | .777 | .908 | .952 | Pass; lowest is log k8 .777 |
| PAM | .571 | .653 | .600 | .659 | **Fail** at k = 7 and 8 on both anchors |
| GBTM | **.632** | .898 | .938 | .901 | **Fail** at log k7 only |

At k = 2-6, all four methods passed.

The task said to stop and report to Jamie if any method failed, so the clusterBMA fits have not started. Jamie has already said the rule is applied as written, with no readmitting of near misses. Read that way, PAM and GBTM drop out and the set becomes k-means and kml, which is 28 inputs per scale. Another reading would cap k at the highest k where all four methods pass, which is k = 6. That choice is Jamie's.

A first run of the check also counted the null-subsample replicates as subsamples. The file filter was fixed before the table above was produced.

The other k = 7-8 replicates for k-means, kml and PAM are done (942). The GBTM replicates for the other perturbations (site, noise, null, starts) are on hold until Jamie decides.

## 2026-11-18: clusterBMA input sets after the k = 2-8 rule check
Decided by Jamie on 2026-11-18 (relayed by the main session; matches the cognition project). Jamie answered "yes" to applying the rule as written and asked not to cap k at 6.
- **Primary (`CBMA_SET=rule`, output/cbma/):** k-means and kml at k = 2-8, using the one-per-family and all-children fits. That is 28 inputs per scale and 84 pooled, with K_BMA = 8.
- **Sensitivity (`CBMA_SET=four`, output/cbma_four/):** all four methods that passed at k = 2-6 (k-means, kml, PAM, GBTM), at k = 2-8. That is 56 inputs per scale; raw has 54 because the two degenerate GBTM all-children fits are left out. Pooled has 166 inputs, with K_BMA = 8.
- The 2026-11-08 all-methods fits (16 fits, K_BMA = 6) were moved to output/cbma_allmethods/. Their summaries, figures (4_clusterbma, 5_single_vs_consensus, 6_class_gradient) and cbma_poverty.rds were moved to output/archive/pre_k28_2026-11-18/. They are archived, not reported.
- R/run_cbma_k28.sh runs the 13 fits for a set. It runs at most two at a time, and the second starts only with memory headroom (memory_pressure free >= 45% and compressor < 3 GB), because another project's targets job holds about 5.5 GB. For the primary set it then runs 12b, 13, 14 and 15. The primary fits started on 2026-11-18 at 22:30; the sensitivity fits follow.

## 2026-11-19: S_Dbw undefined at k = 8: no S_Dbw or averaged fit on raw and log
Decided by Jamie on 2026-11-19 ("OK, let's do those", relayed by the main session). This follows the cognition project's precedent for poa at full n.

clusterCrit's S_Dbw is undefined for some k = 8 inputs:

| Set | Scale | Undefined inputs |
|---|---|---|
| Primary | raw | NaN for all four k = 8 inputs (k-means and kml, one per family and all children) |
| Primary | log | Inf for k-means and kml one-per-family k = 8 |
| Sensitivity (four methods) | raw | Same as the primary set, plus Inf for GBTM one-per-family k = 7 and k = 8 |
| Sensitivity (four methods) | log | Same as the primary set, plus NaN for GBTM all-children k = 8 |

All t inputs are finite in both sets. cbma_weights needs finite, positive values, so the raw S_Dbw fit failed at weighting.

- **Design, both sets, 9 fits each:**
  - Calinski-Harabasz and equal weights on raw, t and log.
  - S_Dbw and averaged on t only.
  - The pooled all-scales run with equal weights.
- Raw and log have no S_Dbw fit and no averaged fit. The averaged fit is not redefined as the mean of CH and equal, since that would be a new rule.
- R/cbma_set.R holds this design, read by 12b, 13, 14 and 15. The four fits that are not run are named in the summaries and figure subtitles. Any other missing fit stops the scripts. R/run_cbma_k28.sh runs the 9 fits.
- The sensitivity set writes to output/cbma_four, results/cbma_poverty_four.rds, figures/4_clusterbma_four and 5_single_vs_consensus_four, and the summaries with a _four suffix.
- Run notes: raw ch took 423 min and t ch 304 min (K_BMA = 8). On 2026-11-19 a second fit started twice on brief memory dips (12:28, then 12:47 with the Mac thrashing). Both were stopped within 20 min. Fits now run one at a time while another project's targets/crew job runs or swap is growing.

## 2026-11-20: reporting the consensus poverty gradient (primary set)
Decided by Jamie on 2026-11-20, relayed by the main session. This applies to the primary clusterBMA set, whose 9 fits finished on 2026-11-19 at 22:30 (12b, 13, 14 and 15 rerun).
- **Main text:** every consensus class against the lowest, on every scale, as 12b_cbma_poverty.R gives it. All six log and T-score fits split off a near-zero lowest class (mean raw CBCL 2.2-2.4, 13-15% of children). Class 2 is less poor than that class, about 0.89 per halving of income (1 / 1.10-1.13 below; CIs exclude 1), so the gradient is U-shaped. Raw and pooled fits rise steadily: top-class RRR 1.41 (CH) and 1.57 (equal) on raw, and 1.37 pooled.
- The U-shape is reported as a property of those scales, as in the 2026-11-12 caveat framing. The reference class is not switched only on the scales that show it, since that would read as a choice made after seeing the results.
- **Supplement:** a class 2 reference table for the log and T-score fits only, from R/cbma_poverty_ref2.R. Output is in results/cbma_poverty_ref2.rds and summaries/summary_cbma_poverty_ref2.txt; the four-method set gets _four versions once its fits finish. Race-adjusted results vs class 2:

  | Fit | Class 1 (near-zero) | Class 3 | Top class |
  |---|---|---|---|
  | log ch | 1.12 (1.06-1.18) | 1.07 | 1.22 (1.16-1.27) |
  | log equal | 1.12 (1.06-1.18) | 1.07 | 1.22 (1.16-1.27) |
  | t ch | 1.12 (1.07-1.19) | 1.05 | 1.25 (1.20-1.31) |
  | t sdbw | 1.10 (1.04-1.15) | 1.08 | 1.27 (1.21-1.32) |
  | t equal (5 classes) | 1.13 (1.07-1.19) | 1.05 | 1.43 (1.34-1.52) |
  | t avg | 1.12 (1.07-1.18) | 1.08 | 1.26 (1.20-1.31) |

  t equal also has class 4 at 1.18.
- **The runner's "Poverty RRR (highest vs lowest)" log line** is 12_cbma.R's stored model, which is age, sex and site with no race (`rrr_unadj` in summary_cbma.txt). For example raw equal shows 1.46 there against 1.57 adjusted, and t equal 1.18 against 1.26. 13, 14 (x1, x5) and 15 (g3) report the race-adjusted model from 12b, as their text says. 12b checks that its unadjusted fit reproduces the stored value.
- The x5 and g3 subtitles now describe the near-zero class on log and T-score. The primary 14 and 15 figures need a rerun to pick that up. That waits until the four-method fits are not using the memory; the four-method run will draw its own 14.

### 2026-11-20 four-method set: class 2 table and 15
- Four-method (sensitivity) set finished 08:45: 9 fits, 12b, 13, 14. `CBMA_SET=four Rscript R/cbma_poverty_ref2.R` gives its class 2 table (summary_cbma_poverty_ref2_four.txt).
- The near-zero lowest class holds in 5 of 6 log/T fits (mean raw CBCL 2.1-2.2, class 1 vs class 2 RRR 1.12-1.15, race-adjusted). The exception is T-score S_Dbw: 3 classes, lowest class mean 5.2 (38% of children), class 1 vs 2 RRR 0.97 (0.93-1.01), so no U-shape. 14's x5 subtitle and the ref2 header now say which fits have it, from the data.
- 15_class_gradient.R is primary-only by design: g1 and g2 use the single-solution universes, which do not depend on the clusterBMA set, and only g3 reads a consensus set; its outputs have no set suffix, so a four-set run would overwrite the primary. The four-set per-class gradient is in 14's x5 (5_single_vs_consensus_four) and summary_extra_figures_four.txt.

### 2026-11-21 GBTM k = 7-8 stability finished
- All 84 remaining GBTM k = 7-8 replicates finished 02:08 (124 in reps_k78/, 0 errors). `09_stability_summary.R full 8 k78` and `10_stability_figures.R full k78` reran; stab1 and stab3 now show k 2-8 for k-means, kml, PAM and GBTM, and R/main_figures.R redrew Main 4.
- .75 rule at k = 2-8 is unchanged: k-means passes (lowest median subsample ARI .772, log k8) and kml passes (.777, log k8). PAM fails (.571, log k7) and GBTM fails (.632, log k7). The primary clusterBMA set stays k-means and kml.

### 2026-11-21 Stability weighting checked as a third scheme: not run
- Idea: weight each input by its median pairwise subsample ARI (output/stability/rule_k28.rds), proportional and normalised like the index weights. It is defined for every input, unlike S_Dbw, and matches the .75 inclusion rule. Not available for T-scores (stability anchors are raw and log only).
- Check (R/cbma_stab_weights_check.R, no fit needed because the consensus matrix is linear in the weights): stability weights are nearly flat (max/min 1.09 rule raw, 1.28 rule log, 1.64 four raw, 1.73 four log; CH 1.98-2.38) and correlate .57-.85 with the CH weights. They move the consensus matrix 9-37% as far from equal weighting as CH does (largest cell change .012-.051 vs CH .104-.118).
- CH and equal fits already agree closely (four log ARI .89), so a stability-weighted fit would sit between them, close to the equal fit. Decision: no stability-weighted fits; the equal-weight fits stand in for it. Methods text: one sentence with these numbers.
