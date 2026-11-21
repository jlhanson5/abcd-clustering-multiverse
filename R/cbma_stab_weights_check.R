# R/cbma_stab_weights_check.R (2026-11-21): stability weighting as a third clusterBMA scheme, checked
#   Rscript R/cbma_stab_weights_check.R
# How far do stability weights move the consensus matrix from equal weights, relative to
# how far CH weights move it? Consensus is linear in weights: C_a - C_b = sum (w_a - w_b) S_m.
# Run from the CBCL root. Aggregates only.
source("00_setup.R")
suppressPackageStartupMessages(library(mclust))
source(file.path("..", "clusterBMA_port", "R", "cbma.R"))
stab <- readRDS("output/stability/rule_k28.rds")$table
res <- readRDS(file.path(RES, "results_full_refit.rds"))
for (SET in c("rule", "four")) for (sc in c("raw", "log")) {
  ch <- readRDS(file.path(OUT, if (SET == "rule") "cbma" else "cbma_four", paste0(sc, "_ch.rds")))
  inp <- ch$inputs
  ari <- stab$median_ari[match(paste(sc, inp$method, inp$k), paste(stab$anchor, stab$method, stab$k))]
  stopifnot(!anyNA(ari))
  w_s <- ari / sum(ari); w_e <- rep(1 / nrow(inp), nrow(inp)); w_c <- inp$weight
  ids <- ch$ids
  D_se <- D_ce <- 0; n2 <- 0
  # accumulate the difference matrices input by input (co-assignment S_m = A A')
  Ds <- matrix(0, length(ids), length(ids)); Dc <- Ds
  for (m in seq_len(nrow(inp))) {
    a <- res$assign[[inp$universe_id[m]]]
    A <- cbma_hard_to_prob(a$class[match(ids, a$src_subject_id)])
    S <- tcrossprod(A)
    Ds <- Ds + (w_s[m] - w_e[m]) * S; Dc <- Dc + (w_c[m] - w_e[m]) * S
  }
  cat(sprintf("%s %s: %d inputs; stability weights ratio %.2f, CH ratio %.2f; ||C_stab - C_equal|| / ||C_ch - C_equal|| = %.2f; max |C_stab - C_equal| = %.4f, max |C_ch - C_equal| = %.4f; cor(w_stab, w_ch) = %.2f\n",
              SET, sc, nrow(inp), max(w_s) / min(w_s), max(w_c) / min(w_c),
              sqrt(sum(Ds^2)) / sqrt(sum(Dc^2)), max(abs(Ds)), max(abs(Dc)), cor(w_s, w_c)))
  rm(Ds, Dc, S); invisible(gc())
}
