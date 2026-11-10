# =============================================================================
# R/traj_stability_test.R: is traj (CRAN 'traj' 3.1.0) stable enough to add as a
# method? Exploratory, added 2026-11-10; not part of the pipeline yet.
#   Rscript R/traj_stability_test.R run [cores]     # restartable, one file per run
#   Rscript R/traj_stability_test.R summary
# traj: 19 measures per trajectory (Leffondre et al. 2004; trajMeasures default, i.e.
# all but the midpoint measure), standardised, then spectral clustering on a binary
# k-nearest-neighbour graph (trajClusters, nstart = 50). Waves equidistant (0-3), as
# for the other non-GBTM/GMM methods. k fixed at 2-6 (no traj k selection).
# Same anchor data and perturbation draws as the stability analysis (R/perturb.R):
#   seed       3 fit seeds on the full anchor data (cut from 10: the first two were identical)
#   subsample  20 80% subsamples (draws 1-20, the ones every other method got)
#   null       5 null datasets (draws 1-5)
#   cap        1 full-data fit with cap.outliers = TRUE (settings check)
# Metric, as in 09_stability_summary.R: ARI between two solutions on the children
# they share, all pairs of replicates within (anchor, perturbation, k); plus each
# replicate vs the seed-1 full-data fit. Compared with the other methods' pairs from
# output/stability/summary.rds.
# Output: output/traj_test/*.rds; output/summaries/summary_traj_test.txt
# =============================================================================

source("00_setup.R")
source("R/perturb.R")
suppressPackageStartupMessages({ library(traj); library(parallel) })

args <- commandArgs(trailingOnly = TRUE)
TDIR <- file.path(OUT, "traj_test"); dir.create(TDIR, showWarnings = FALSE)
SOUT <- file.path(OUT, "stability")
base <- readRDS(file.path(SOUT, "anchor.rds"))
nm   <- readRDS(file.path(SOUT, "null_model.rds"))
N    <- length(base$ids)

RUNS <- bind_rows(
  tibble(ptype = "seed", rep = 1:3),
  tibble(ptype = "subsample", rep = 1:20),
  tibble(ptype = "null", rep = 1:5),
  tibble(ptype = "cap", rep = 1L)) |>
  tidyr::crossing(anchor = STAB_ANCHORS) |>
  arrange(match(ptype, c("seed", "subsample", "cap", "null")), match(anchor, STAB_ANCHORS), rep) |>
  mutate(file = file.path(TDIR, sprintf("%s_%s_r%02d.rds", anchor, ptype, rep)))

run_one <- function(i) {
  r <- RUNS[i, ]
  if (file.exists(r$file)) return(invisible(NULL))
  idx <- if (r$ptype == "subsample") draw_subsample(N, r$rep) else seq_len(N)
  R <- if (r$ptype == "null") simulate_null(base, nm, r$rep) else base$R
  X <- R[idx, , drop = FALSE]; if (r$anchor == "log") X <- log1p(X)
  t0 <- Sys.time()
  m <- trajMeasures(unname(X), cap.outliers = r$ptype == "cap")
  cls <- sapply(K_STAB, \(k) {
    set.seed(if (r$ptype == "seed") SEED + r$rep else SEED + 1)
    tryCatch(as.integer(trajClusters(m, nclusters = k)$partition[, 2]), error = \(e) rep(NA_integer_, length(idx)))
  })
  colnames(cls) <- paste0("k", K_STAB)
  out <- list(run = r, idx = idx, cls = cls, mins = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  saveRDS(out, paste0(r$file, ".tmp")); file.rename(paste0(r$file, ".tmp"), r$file)
  cat(sprintf("%s %s %s r%02d done (%.1f min)\n", format(Sys.time(), "%m-%d %H:%M"), r$anchor, r$ptype, r$rep, out$mins))
}

if (length(args) >= 1 && args[1] == "run") {
  cores <- if (length(args) >= 2) as.integer(args[2]) else 1
  invisible(mclapply(seq_len(nrow(RUNS)), run_one, mc.cores = cores, mc.preschedule = FALSE))
  quit(save = "no")
}

# ---- Summary -------------------------------------------------------------------------------
ari <- function(x, y) {   # ARI on shared, non-missing children; x, y full-length (N) vectors
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) < 2) return(NA_real_)
  mclust::adjustedRandIndex(x[ok], y[ok])
}
fullv <- function(idx, v) { out <- rep(NA_integer_, N); out[idx] <- v; out }
recs <- lapply(RUNS$file[file.exists(RUNS$file)], readRDS)
sumf <- file.path(SUMM, "summary_traj_test.txt"); unlink(sumf)
say  <- function(...) { cat(..., "\n"); cat(..., "\n", file = sumf, append = TRUE) }
pt   <- function(x) paste(capture.output(print(x, n = 200, width = 200)), collapse = "\n")
say(sprintf("traj stability test (%s): %d of %d runs; median %.1f min per run (k 2-6).",
            format(Sys.time(), "%Y-%m-%d %H:%M"), length(recs), nrow(RUNS), median(sapply(recs, `[[`, "mins"))))

V <- bind_rows(lapply(recs, \(x) tibble(anchor = x$run$anchor, ptype = x$run$ptype, rep = x$run$rep,
                                        k = K_STAB, v = lapply(seq_along(K_STAB), \(j) fullv(x$idx, x$cls[, j])))))
ref <- V |> filter(ptype == "seed", rep == 1) |> select(anchor, k, ref = v)
pairs <- V |> group_by(anchor, ptype, k) |> group_modify(\(g, key) {
  if (nrow(g) < 2) return(tibble(ari = numeric()))
  ij <- t(combn(nrow(g), 2)); tibble(ari = apply(ij, 1, \(p) ari(g$v[[p[1]]], g$v[[p[2]]])))
}) |> ungroup()
vsref <- V |> inner_join(ref, by = c("anchor", "k")) |> filter(!(ptype == "seed" & rep == 1)) |>
  mutate(ari = map2_dbl(v, ref, ari)) |> select(anchor, ptype, rep, k, ari)
sizes <- V |> mutate(min_share = map_dbl(v, \(x) { t <- table(x); min(t) / sum(t) }),
                     n_cls = map_int(v, \(x) length(unique(na.omit(x))))) |>
  group_by(anchor, ptype, k) |> summarise(min_class_share = round(median(min_share), 3), k_found_lt_k = sum(n_cls < k), .groups = "drop")

say("\nPairwise ARI between replicates (shared children), median [IQR], traj:")
tab <- pairs |> group_by(anchor, ptype, k) |>
  summarise(n_pairs = n(), median = round(median(ari, na.rm = TRUE), 3), q25 = round(quantile(ari, .25, na.rm = TRUE), 3),
            q75 = round(quantile(ari, .75, na.rm = TRUE), 3), .groups = "drop") |> left_join(sizes, by = c("anchor", "ptype", "k"))
say(pt(tab))
say("\nEach run vs the seed-1 full-data fit (median ARI):")
say(pt(vsref |> group_by(anchor, ptype, k) |> summarise(median = round(median(ari, na.rm = TRUE), 3), .groups = "drop") |>
         tidyr::pivot_wider(names_from = k, values_from = median, names_prefix = "k")))

# The other methods, same metric and the same draws (seed reps 1-50 / subsample draws 1-100 there).
st <- readRDS(file.path(SOUT, "summary.rds"))
other <- st$pairs |> filter(variant == "base", j > 0, ptype %in% c("seed", "subsample", "null")) |>
  group_by(anchor, method, ptype, k) |> summarise(median = round(median(ari, na.rm = TRUE), 3), .groups = "drop")
cmp <- bind_rows(other, tab |> transmute(anchor, method = "traj", ptype, k, median)) |>
  filter(ptype %in% c("seed", "subsample", "null")) |>
  tidyr::pivot_wider(names_from = method, values_from = median) |> arrange(anchor, ptype, k)
say("\nMedian pairwise ARI, traj next to the other methods (same anchors and perturbations):")
say(pt(cmp))
