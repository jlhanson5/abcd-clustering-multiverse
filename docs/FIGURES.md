# Figure index

Every figure in `output/figures/`, by folder. Since 2026-11-09 every poverty RRR is
adjusted for age, sex, site and race/ethnicity (primary). The age/sex/site-only model
is the sensitivity analysis: it is shown in fig13, as open points in cbma5, and in the
summaries (`*_unadj`). Pre-race versions of every figure are in
`output/archive/pre_race_2026-11-09/figures/` (old folder names).

RRR = relative risk ratio of being in the highest- vs lowest-symptom class per halving of
income. Universes whose highest class has fewer than 10 children have no RRR (not estimable).

## 1_spec_curve: the primary run (`full_refit`; 06_figures.R + 11_refit_figures.R)

| Figure | Shows | Role |
|---|---|---|
| fig1_spec_curve | RRR in every universe (k 2-6), sorted, with the analytic choices below | **Main result** |
| fig2_agreement_methods | Mean ARI between methods, same data and k | Main |
| fig3_agreement_full | ARI between every pair of universes | Supplement |
| fig4_alluvial | The same children across 8 methods, one configuration, k = 4 | Main |
| fig5_k_selection | k chosen by BIC and silhouette; k-max 5 vs 6 | Supplement |
| fig6_diagnostics | Fit failures, empty classes, smallest class | Supplement |
| fig7_decision_importance | Share of variance in log RRR explained by each choice | Main |
| fig8_robustness | RRR median, IQR and % CI above 1 across subsets of the multiverse | Main |
| fig9_rrr_by_method | Every universe's RRR, by method and scale | Main |
| fig10_top_class | Size of the highest class vs RRR (why scale matters) | Main |
| fig11_trajectories | Class mean CBCL by wave, anchor configuration, k = 3 | Main |
| fig12_gmm_refit | GMM: main-run vs refit RRR, and how often starts converge | Supplement |
| fig13_race_sensitivity | RRR with vs without race/ethnicity adjustment, every universe | **Sensitivity** |
| figS1_class_gradient | Every class's RRR vs its mean CBCL, by method | Supplement |
| figS2_k_extension | RRR at k 2-8, five fast methods | Supplement |

## 2_main_run: before the GMM and Gaussian-mixture refits (`full`; 06_figures.R)

fig1-7, S1, S2 as above, for the original main run. **Sensitivity** only.

## 3_stability: replicates (09/10_stability_*.R)

| Figure | Shows | Role |
|---|---|---|
| stab1_heatmap | Mean ARI between replicates by method and k, per perturbation | Main |
| stab2_children | Share of children reliably in the highest class vs on the boundary | Main |
| stab3_null | Real vs null data: silhouette, entropy, k selected | Main |
| stab4_poverty | RRR across subsamples, site halves and noise, against the null | Main |
| stabS1_starts | GBTM/GMM: log-likelihood of each start vs the best | Supplement |
| stabS2_settings | Settings variants vs base | Supplement |

## 4_clusterbma: consensus across methods and k (13_cbma_summary.R)

| Figure | Shows | Role |
|---|---|---|
| cbma1_weights | Weight of every input, per scale and weighting | Main |
| cbma2_trajectories | Consensus classes: mean CBCL by wave, class sizes | Main |
| cbma3_uncertainty | Per-child allocation uncertainty | Supplement |
| cbma4_agreement | ARI of each input with the consensus | Supplement |
| cbma5_poverty | RRR for all 16 consensus fits (open points: without race) | Main |
| cbma6_weighting_agreement | ARI between weightings on the same scale | Supplement |
| cbma7_uncertainty_by_class | Mean uncertainty per consensus class | Supplement |

## 5_single_vs_consensus: single methods next to clusterBMA (14_extra_figures.R)

| Figure | Shows | Role |
|---|---|---|
| x1_rrr_inputs_vs_consensus | RRR of every input solution and every consensus fit | Main |
| x2_class_map | Size and mean CBCL of every class | Supplement |
| x3_crosswalk | Consensus classes vs each method's k = 3 classes | Supplement |
| x4_child_uncertainty | Child uncertainty vs the child's mean CBCL | Supplement |
| x5_class_gradient | RRR of every class vs the lowest, against class mean CBCL | Main |
| x6_highest_share | Share of inputs placing each child in the highest class | Supplement |
| x7_consensus_agreement | ARI between all 16 consensus fits | Supplement |

"Main" and "Supplement" are suggestions for the write-up, not decisions.

## 6_class_gradient: every class vs the lowest (15_class_gradient.R)

The primary contrast stays highest vs lowest-mean class; these show the whole gradient.

| Figure | Shows | Role |
|---|---|---|
| g1_gradient_lines | Each universe's RRR across its classes, against class mean CBCL; scale x k | Main |
| g2_gradient_shape | Share of universes whose gradient rises steadily vs has a middle class below the lowest; scale x method | Main |
| g3_consensus_gradient | Every class of the 16 clusterBMA fits, 95% CI | Supplement |

## 7_sensitivity: alternative models and the permutation test (16-18)

| Figure | Shows | Script |
|---|---|---|
| s1_site_random_effects | Each universe's RRR, primary model vs binary logit with site as a random effect; by scale | 16_sens_site_re.R |
| s1_site_random_effects_famRE_attempt | Same, from the first run with a family random effect in full-sample universes (not used; see README) | 16_sens_site_re.R (earlier version) |
| s2_permutation | Observed median RRR, share with CI above 1 and Stouffer's Z against 500 permuted curves | 17_permutation.R summary |
| s3_population_weights | Each universe's RRR, unweighted vs weighted to the ACS population; by scale | 18_sens_weights.R |
