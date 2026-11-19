# =============================================================================
# R/cbma_rule_k28.R: the clusterBMA inclusion rule at k = 2-8 (decided 2026-11-15,
# docs/DECISIONS_memory.md): median pairwise subsample ARI >= .75 at every k on both
# stability anchors (raw, log). k = 2-6 from output/stability/summary.rds (09); k = 7-8
# computed here from the subsample replicates in output/stability/reps_k78/ (08 with
# k=7|8), with 09's ARI on shared children, so the rule can be checked before the rest of
# the k = 7-8 replicates finish. Aggregates only.
#   Rscript R/cbma_rule_k28.R
# Output: output/summaries/summary_cbma_rule_k28.txt, output/stability/rule_k28.rds
# =============================================================================

source("00_setup.R")
source("R/perturb.R")

RULE <- 0.75
METH <- c("kmeans", "kml", "pam", "gbtm")
SOUT <- file.path(OUT, "stability")
base <- readRDS(file.path(SOUT, "anchor.rds")); N <- length(base$ids)
fullv <- function(idx, cls) { v <- rep(NA_integer_, N); v[idx] <- cls; v }
ari <- function(x, y) {            # 09_stability_summary.R's formula
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]
  n <- length(x); if (n < 2) return(NA_real_)
  x <- match(x, sort(unique(x))); y <- match(y, sort(unique(y)))
  kx <- max(x); ky <- max(y)
  tab <- matrix(tabulate((x - 1L) * ky + y, kx * ky), kx, ky, byrow = TRUE)
  c2 <- function(v) sum(v * (v - 1) / 2)
  sa <- c2(rowSums(tab)); sb <- c2(colSums(tab)); ex <- sa * sb / (n * (n - 1) / 2)
  den <- (sa + sb) / 2 - ex
  if (den == 0) NA_real_ else (c2(tab) - ex) / den
}

lo <- readRDS(file.path(SOUT, "summary.rds"))$pairs |>
  filter(ptype == "subsample", variant == "base", j > 0, method %in% METH) |>
  group_by(anchor, method, k) |>
  summarise(reps = n_distinct(c(i, j)), pairs = n(), median_ari = median(ari, na.rm = TRUE), .groups = "drop")

files <- list.files(file.path(SOUT, "reps_k78"), "^(raw|log)_[a-z]+_subsample_base_r\\d+\\.rds$", full.names = TRUE)
sol <- bind_rows(lapply(files, function(f) {
  r <- readRDS(f); if (r$rep$ptype != "subsample") return(NULL)
  bind_rows(lapply(r$fits, \(o) if (length(o$parts$all))
    tibble(anchor = r$rep$anchor, method = r$rep$method, rep = r$rep$rep, k = o$fit$k[1],
           v = list(fullv(o$parts$all$idx, o$parts$all$cls)))))
}))
hi <- if (!nrow(sol)) NULL else sol |> group_by(anchor, method, k) |> group_modify(function(g, key) {
  R <- nrow(g); ij <- t(combn(R, 2))
  a <- vapply(seq_len(nrow(ij)), \(p) ari(g$v[[ij[p, 1]]], g$v[[ij[p, 2]]]), 0)
  tibble(reps = R, pairs = length(a), median_ari = median(a, na.rm = TRUE))
}) |> ungroup()

expect <- c(kmeans = 100, kml = 100, pam = 100, gbtm = 20)
tab <- bind_rows(lo, hi) |>
  mutate(complete = reps == expect[method], pass = median_ari >= RULE) |>
  arrange(match(method, METH), anchor, k)
verdict <- tab |> group_by(method) |>
  summarise(ks = paste(sort(unique(k)), collapse = ","), all_complete = all(complete) && n() == 14,
            min_median = min(median_ari), at = paste0(anchor[which.min(median_ari)], " k", k[which.min(median_ari)]),
            pass = all(pass) && n() == 14, .groups = "drop")
saveRDS(list(table = tab, verdict = verdict, rule = RULE), file.path(SOUT, "rule_k28.rds"))

sumf <- file.path(SUMM, "summary_cbma_rule_k28.txt")
out <- c(sprintf("clusterBMA inclusion rule, k = 2-8 (%s): median pairwise subsample ARI >= %.2f at every k, raw and log anchors.",
                 format(Sys.time(), "%Y-%m-%d %H:%M"), RULE),
         "k = 2-6 from output/stability/summary.rds; k = 7-8 from output/stability/reps_k78/.", "",
         capture.output(print(tab |> select(-pass, -complete) |> mutate(median_ari = round(median_ari, 3)) |>
                                pivot_wider(id_cols = c(method, anchor), names_from = k, values_from = median_ari,
                                            names_prefix = "k"), n = 50, width = 200)), "",
         "Replicates per cell (k = 7-8):",
         capture.output(print(filter(tab, k >= 7) |> select(method, anchor, k, reps, pairs), n = 50)), "",
         capture.output(print(verdict, width = 200)))
writeLines(out, sumf); cat(out, sep = "\n")
