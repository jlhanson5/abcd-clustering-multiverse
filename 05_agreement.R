# =============================================================================
# 05_agreement.R: pairwise agreement between every pair of universes.
#   Rscript 05_agreement.R pilot|full
# Adjusted Rand index (ARI) and normalised mutual information (NMI, sqrt
# normalisation), computed on the children the two universes share. Agreement
# does not depend on class labels, so the relabelling in 04 does not matter here.
# Output: output/results/agreement_<mode>.rds (ARI, NMI, n_common matrices; universe order).
# =============================================================================

source("00_setup.R")
suppressPackageStartupMessages(library(parallel))

mode <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(mode)) mode <- "pilot"
res  <- readRDS(file.path(RES, paste0("results_", mode, ".rds")))

keep <- res$fits |>
  filter(is.na(error), !degenerate | is.na(degenerate)) |>
  arrange(match(method, c("kmeans", "pam", "ward", "llpa", "kml", "gbtm", "gmm")),
          k, scale, time, family, config)
uids <- keep$universe_id
ids  <- sort(unique(unlist(lapply(res$assign[uids], `[[`, "src_subject_id"))))
A <- sapply(uids, \(u) { a <- res$assign[[u]]; a$class[match(ids, a$src_subject_id)] })

pair_stats <- function(x, y) {
  kx <- max(x); ky <- max(y)
  tab <- matrix(tabulate((x - 1L) * ky + y, kx * ky), kx, ky, byrow = TRUE)
  n <- length(x)
  # ARI
  c2 <- function(v) sum(v * (v - 1) / 2)
  sij <- c2(tab); sa <- c2(rowSums(tab)); sb <- c2(colSums(tab)); tot <- n * (n - 1) / 2
  expct <- sa * sb / tot
  ari <- (sij - expct) / ((sa + sb) / 2 - expct)
  # NMI
  p <- tab / n; px <- rowSums(p); py <- colSums(p)
  nz <- p > 0
  mi <- sum(p[nz] * log(p[nz] / outer(px, py)[nz]))
  hx <- -sum(px[px > 0] * log(px[px > 0])); hy <- -sum(py[py > 0] * log(py[py > 0]))
  c(ari = ari, nmi = mi / sqrt(hx * hy), n = n)
}

U <- length(uids)
t0 <- Sys.time()
rows <- mclapply(seq_len(U), function(i) {
  out <- matrix(NA_real_, U, 3)
  for (j in seq_len(U)) {
    if (j < i) next
    ok <- !is.na(A[, i]) & !is.na(A[, j])
    out[j, ] <- pair_stats(A[ok, i], A[ok, j])
  }
  out
}, mc.cores = 10)

ARI <- NMI <- NC <- matrix(NA_real_, U, U, dimnames = list(uids, uids))
for (i in seq_len(U)) {
  ARI[i, ] <- rows[[i]][, 1]; NMI[i, ] <- rows[[i]][, 2]; NC[i, ] <- rows[[i]][, 3]
}
sym <- function(M) { M[lower.tri(M)] <- t(M)[lower.tri(M)]; M }
ARI <- sym(ARI); NMI <- sym(NMI); NC <- sym(NC)

saveRDS(list(ARI = ARI, NMI = NMI, n_common = NC, universes = keep),
        file.path(RES, paste0("agreement_", mode, ".rds")))
cat(sprintf("Agreement: %d universes, %d pairs, %.1f min\n", U, U * (U - 1) / 2,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
