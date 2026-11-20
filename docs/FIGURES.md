# Figure index

Every figure in `output/figures/`, by folder. Since 2026-11-09 every poverty RRR is
adjusted for age, sex, site and race/ethnicity (primary). The age/sex/site-only model
is the sensitivity analysis: it is shown in fig13, as open points in cbma5, and in the
summaries (`*_unadj`). Pre-race versions of every figure are in
`output/archive/pre_race_2026-11-09/figures/` (old folder names).

RRR = relative risk ratio of being in the highest- vs lowest-symptom class per halving of
income. Universes whose highest class has fewer than 10 children have no RRR (not estimable).

Roles decided by Jamie 2026-11-18 (docs/DECISIONS_memory.md): six main figures, everything else
supplement, including the whole four-method clusterBMA sensitivity set (`*_four` folders).

## Main figures

| Main | Figure | File |
|---|---|---|
| 1 | Specification curve | `1_spec_curve/fig1_spec_curve.png` |
| 2 | Decision importance | `1_spec_curve/fig7_decision_importance.png` |
| 3 | Highest-class size vs RRR | `1_spec_curve/fig10_top_class.png` |
| 4 | Stability: A stab1_heatmap, B stab3_null | `0_main/main4_stability.png` |
| 5 | clusterBMA, primary set: A cbma2_trajectories, B g3_consensus_gradient (every class vs the lowest) | `0_main/main5_clusterbma.png` |
| 6 | Poverty gradient across classes | `6_class_gradient/g1_gradient_lines.png` |

Main 4 and 5 are drawn by `R/main_figures.R` from the plot objects (`*.ggrds`) that
10_stability_figures.R, 13_cbma_summary.R and 15_class_gradient.R save next to each PNG; rerun
it after any of them. Main 5B shows the whole gradient, not highest vs lowest (2026-11-20);
its subtitle names the lowest class as the reference.
Their single panels stay in 3_stability and 4_clusterbma.

## 1_spec_curve: the primary run (`full_refit`; 06_figures.R + 11_refit_figures.R)

| Figure | Shows | Role |
|---|---|---|
| fig1_spec_curve | RRR in every universe (k 2-6), sorted, with the analytic choices below | **Main 1** |
| fig2_agreement_methods | Mean ARI between methods, same data and k | Supplement |
| fig3_agreement_full | ARI between every pair of universes | Supplement |
| fig4_alluvial | The same children across 8 methods, one configuration, k = 4 | Supplement |
| fig5_k_selection | k chosen by BIC and silhouette; k-max 5 vs 6 | Supplement |
| fig6_diagnostics | Fit failures, empty classes, smallest class | Supplement |
| fig7_decision_importance | Share of variance in log RRR explained by each choice | **Main 2** |
| fig8_robustness | RRR median, IQR and % CI above 1 across subsets of the multiverse | Supplement |
| fig9_rrr_by_method | Every universe's RRR, by method and scale | Supplement |
| fig10_top_class | Size of the highest class vs RRR (why scale matters) | **Main 3** |
| fig11_trajectories | Class mean CBCL by wave, anchor configuration, k = 3 | Supplement |
| fig12_gmm_refit | GMM: main-run vs refit RRR, and how often starts converge | Supplement |
| fig13_race_sensitivity | RRR with vs without race/ethnicity adjustment, every universe | Supplement |
| figS1_class_gradient | Every class's RRR vs its mean CBCL, by method | Supplement |
| figS2_k_extension | RRR at k 2-8, five fast methods | Supplement |

## 2_main_run: before the GMM and Gaussian-mixture refits (`full`; 06_figures.R)

fig1-7, S1, S2 as above, for the original main run. **Sensitivity** only.

## 3_stability: replicates (09/10_stability_*.R)

| Figure | Shows | Role |
|---|---|---|
| stab1_heatmap | Mean ARI between replicates by method and k, per perturbation | Main 4A |
| stab2_children | Share of children reliably in the highest class vs on the boundary | Supplement |
| stab3_null | Real vs null data: silhouette, entropy, k selected | Main 4B |
| stab4_poverty | RRR across subsamples, site halves and noise, against the null | Supplement |
| stabS1_starts | GBTM/GMM: log-likelihood of each start vs the best | Supplement |
| stabS2_settings | Settings variants vs base | Supplement |

## 4_clusterbma: consensus across methods and k (13_cbma_summary.R)

Primary set: k-means and kml, k 2-8 (28 inputs per scale). Nine fits: Calinski-Harabasz and
equal on every scale, S_Dbw and averaged on T-score only, pooled all-scales with equal weights
(raw and log S_Dbw / averaged not run by design, 2026-11-19). `4_clusterbma_four/` has the same
figures for the four-method sensitivity set (k-means, kml, PAM, GBTM): **Supplement**.

| Figure | Shows | Role |
|---|---|---|
| cbma1_weights | Weight of every input, per scale and weighting | Supplement |
| cbma2_trajectories | Consensus classes: mean CBCL by wave, class sizes | Main 5A |
| cbma3_uncertainty | Per-child allocation uncertainty | Supplement |
| cbma4_agreement | ARI of each input with the consensus | Supplement |
| cbma5_poverty | RRR, highest vs lowest class, for the 9 consensus fits (open points: without race) | Supplement |
| cbma6_weighting_agreement | ARI between weightings on the same scale | Supplement |
| cbma7_uncertainty_by_class | Mean uncertainty per consensus class | Supplement |

## 5_single_vs_consensus: single methods next to clusterBMA (14_extra_figures.R)

| Figure | Shows | Role |
|---|---|---|
| x1_rrr_inputs_vs_consensus | RRR of every input solution and every consensus fit | Supplement |
| x2_class_map | Size and mean CBCL of every class | Supplement |
| x3_crosswalk | Consensus classes vs each method's k = 3 classes | Supplement |
| x4_child_uncertainty | Child uncertainty vs the child's mean CBCL | Supplement |
| x5_class_gradient | RRR of every class vs the lowest, against class mean CBCL | Supplement |
| x6_highest_share | Share of inputs placing each child in the highest class | Supplement |
| x7_consensus_agreement | ARI between the 9 consensus fits | Supplement |

`5_single_vs_consensus_four/`: the same for the four-method set, **Supplement**.

## 6_class_gradient: every class vs the lowest (15_class_gradient.R)

The primary contrast stays highest vs lowest-mean class; these show the whole gradient.
Supplementary table: log and T-score consensus classes vs class 2
(`summaries/summary_cbma_poverty_ref2.txt`, `_four` for the sensitivity set; R/cbma_poverty_ref2.R).

| Figure | Shows | Role |
|---|---|---|
| g1_gradient_lines | Each universe's RRR across its classes, against class mean CBCL; scale x k | **Main 6** |
| g2_gradient_shape | Share of universes whose gradient rises steadily vs has a middle class below the lowest; scale x method | Supplement |
| g3_consensus_gradient | Every class of the 9 primary clusterBMA fits vs the lowest, 95% CI (primary set only) | Main 5B |

## 7_sensitivity: alternative models and the permutation test (16-18)

| Figure | Shows | Script | Role |
|---|---|---|---|
| s1_site_random_effects | Each universe's RRR, primary model vs binary logit with site as a random effect; by scale | 16_sens_site_re.R | Supplement |
| s1_site_random_effects_famRE_attempt | Same, from the first run with a family random effect in full-sample universes (not used; see README) | 16_sens_site_re.R (earlier version) | Supplement |
| s2_permutation | Observed median RRR, share with CI above 1 and Stouffer's Z against 500 permuted curves | 17_permutation.R summary | Supplement (statistic reported in the text) |
| s3_population_weights | Each universe's RRR, unweighted vs weighted to the ACS population; by scale | 18_sens_weights.R | Supplement |
