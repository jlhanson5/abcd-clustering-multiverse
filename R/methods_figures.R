# =============================================================================
# R/methods_figures.R (2026-11-22): methods figures. Text and simulated data only.
#   m1_pipeline.png    the analysis pipeline (counts from docs/DECISIONS_memory.md)
#   m2_clusterbma.png  how clusterBMA combines solutions (toy data; ../clusterBMA_port)
#   m3_stability.png   the stability rule on real-like and null data (toy data)
#   Rscript R/methods_figures.R
# Output: output/figures/0_methods/
# =============================================================================

PORT <- normalizePath(file.path("..", "clusterBMA_port"))
source(file.path(PORT, "tools", "pipeline_diagram.R"))
figd <- file.path("output", "figures", "0_methods"); dir.create(figd, recursive = TRUE, showWarnings = FALSE)

stages <- list(
  list(title = "Data", lines = c(
    "ABCD Study, CBCL total problems at baseline and years 1, 2 and 3",
    "Anchor sample: 8,331 children, one per family, at least 3 waves, missing waves imputed",
    "Family income, age, sex, site, race/ethnicity")),
  list(title = "Clustering multiverse", lines = c(
    "852 universes",
    "8 methods: k-means, PAM, Ward, longitudinal LPA, kml, GBTM, GBTM (link), growth mixture",
    "k = 2-8",
    "Scale: raw, T-score, log(1+raw)",
    "Waves: all 4, at least 3 imputed, at least 3 observed",
    "One child per family or all children; wave or age as time")),
  list(title = "Stability", lines = c(
    "80% subsamples, site split-halves, added noise, null data, seeds and settings",
    "Inclusion rule: median pairwise subsample ARI >= .75 at every k = 2-8 (raw and log)",
    "Passes: k-means, kml",
    "Null data (no subgroups) are about as stable: stability shows reproducibility, not structure")),
  list(title = "clusterBMA consensus", lines = c(
    "Primary inputs: k-means and kml, k = 2-8, one-per-family and all-children fits: 28 per scale",
    "Weights: Calinski-Harabasz, S_Dbw, equal, averaged; all scales pooled",
    "Up to 8 classes; each child gets class probabilities",
    "9 fits; sensitivity set adds PAM and GBTM (9 fits)")),
  list(title = "Poverty gradient", lines = c(
    "Multinomial model of consensus class on income (per halving), age, sex, site, race/ethnicity",
    "Every class vs the lowest class",
    "Supplement: class 2 as the reference on log and T-score")))
checks <- c(
  "Consensus method: evidence accumulation and median partition on the same inputs and weights",
  "Input set: four-method sensitivity set (adds PAM and GBTM)",
  "Weighting: stability weights (nearly equal; equal-weight fits stand in)",
  "Site random effects, permutation test, population weights",
  "Single solutions vs the consensus; race-unadjusted models")
draw_pipeline(stages, checks, title = "Analysis pipeline",
              subtitle = "Many clustering solutions, filtered by stability, combined by clusterBMA, then related to family income.",
              file = file.path(figd, "m1_pipeline.png"), width = 16, wrap = 44)

# Toy-data figures: drawn by the port, copied here
owd <- setwd(PORT); system2("Rscript", "tools/methods_figures.R", stdout = FALSE); setwd(owd)
file.copy(file.path(PORT, "figures", "m_clusterbma.png"), file.path(figd, "m2_clusterbma.png"), overwrite = TRUE)
file.copy(file.path(PORT, "figures", "m_stability.png"), file.path(figd, "m3_stability.png"), overwrite = TRUE)
cat("Figures written to", figd, "\n")
