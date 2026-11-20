# =============================================================================
# R/main_figures.R: the two combined main-text figures (split decided 2026-11-18,
# docs/DECISIONS_memory.md; list in docs/FIGURES.md).
#   Main 4  stab1_heatmap (A) over stab3_null (B)          (10_stability_figures.R)
#   Main 5  cbma2_trajectories (A) over cbma5_poverty (B)  (13_cbma_summary.R, primary set)
# Reads the plot objects 10 and 13 save next to each PNG (<name>.ggrds), so run those
# first. Aggregates only.
#   Rscript R/main_figures.R
# Output: output/figures/0_main/main4_stability.png, main5_clusterbma.png
# =============================================================================

source("00_setup.R")
suppressPackageStartupMessages(library(patchwork))

figd <- fig_dir("0_main"); dir.create(figd, recursive = TRUE, showWarnings = FALSE)
obj <- function(dir, name) {
  f <- file.path(FIG, dir, paste0(name, ".ggrds"))
  if (!file.exists(f)) stop(sprintf("%s missing: rerun the script that draws %s", f, name))
  readRDS(f)
}
tag_theme <- theme(plot.tag = element_text(face = "bold", size = 16))

combine <- function(a, b, out) {
  theme_set(a$theme)   # 10 and 13 use the same theme
  w <- max(a$w, b$w); h <- a$h + b$h
  # wrap_elements keeps an inner patchwork (stab3 has two rows) as one tagged panel
  p <- wrap_elements(full = a$plot) / wrap_elements(full = b$plot) +
    plot_layout(heights = c(a$h, b$h)) + plot_annotation(tag_levels = "A") & tag_theme
  ggsave(file.path(figd, paste0(out, ".png")), p, width = w, height = h, dpi = 300)
  cat(sprintf("%s: %.1f x %.1f in\n", file.path(figd, paste0(out, ".png")), w, h))
}

combine(obj(FIG_DIR[["stability"]], "stab1_heatmap"), obj(FIG_DIR[["stability"]], "stab3_null"), "main4_stability")
combine(obj(FIG_DIR[["cbma"]], "cbma2_trajectories"), obj(FIG_DIR[["cbma"]], "cbma5_poverty"), "main5_clusterbma")
