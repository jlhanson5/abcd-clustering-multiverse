# =============================================================================
# R/traj_variants.R: follow-up to R/traj_stability_test.R (added 2026-11-10).
#   Rscript R/traj_variants.R run [cores]     # restartable, one file per run
#   Rscript R/traj_variants.R summary
# The base traj test split off tiny classes (raw: the 44 flat trajectories, extreme on
# measure 9; log: ~10 children with identical patterns, cut off in the kNN graph).
# Two fixes, and both together:
#   nolin   drop measures 8 and 9 (variance explained by, and rate of crossing, the
#           affine approximation; both degenerate for a flat trajectory)
#   jitter  add U(-0.05, 0.05) noise to every score on the anchor scale (after log1p for
#           log) to break ties; one draw for the full anchor data (fixed seed), so a child
#           has the same noise in every replicate and the noise is not a source of
#           instability between replicates
#   both    nolin + jitter
#   jitcap  jitter + cap.outliers = TRUE (all measures); added 2026-11-10 after nolin
#           failed (its tiny classes are children with tied measures, which capping
#           cannot separate)
# Runs per variant and anchor (raw, log): the seed-1 full-data fit, the same 20
# subsamples and 5 null datasets as the base test. 208 runs, ~7 min each on 1 core.
# Output: output/traj_test/variants/*.rds; output/summaries/summary_traj_variants.txt
# =============================================================================

source("00_setup.R")
source("R/perturb.R")
suppressPackageStartupMessages({ library(traj); library(parallel) })

args <- commandArgs(trailingOnly = TRUE)
BDIR <- file.path(OUT, "traj_test")
VDIR <- file.path(BDIR, "variants"); dir.create(VDIR, showWarnings = FALSE, recursive = TRUE)
SOUT <- file.path(OUT, "stability")
base <- readRDS(file.path(SOUT, "anchor.rds"))
nm   <- readRDS(file.path(SOUT, "null_model.rds"))
N    <- length(base$ids)
ALL_M  <- c(1:10, 12:20)                 # trajMeasures default
VARIANTS <- list(nolin = list(meas = setdiff(ALL_M, 8:9), jit = FALSE),
                 jitter = list(meas = ALL_M, jit = TRUE),
                 both = list(meas = setdiff(ALL_M, 8:9), jit = TRUE),
                 jitcap = list(meas = ALL_M, jit = TRUE, cap = TRUE))
JIT <- local({ set.seed(SEED + 9e5); matrix(runif(N * ncol(base$R), -0.05, 0.05), N) })

RUNS <- bind_rows(tibble(ptype = "seed", rep = 1L), tibble(ptype = "subsample", rep = 1:20),
                  tibble(ptype = "null", rep = 1:5)) |>
  tidyr::crossing(variant = names(VARIANTS), anchor = STAB_ANCHORS) |>
  arrange(match(ptype, c("seed", "subsample", "null")), rep, match(variant, names(VARIANTS)), match(anchor, STAB_ANCHORS)) |>
  mutate(file = file.path(VDIR, sprintf("%s_%s_%s_r%02d.rds", variant, anchor, ptype, rep)))

run_one <- function(i) {
  r <- RUNS[i, ]; v <- VARIANTS[[r$variant]]
  if (file.exists(r$file)) return(invisible(NULL))
  idx <- if (r$ptype == "subsample") draw_subsample(N, r$rep) else seq_len(N)
  R <- if (r$ptype == "null") simulate_null(base, nm, r$rep) else base$R
  X <- if (r$anchor == "log") log1p(R) else R
  if (v$jit) X <- X + JIT
  X <- X[idx, , drop = FALSE]
  t0 <- Sys.time()
  m <- trajMeasures(unname(X), measures = v$meas, cap.outliers = isTRUE(v$cap))
  cls <- sapply(K_STAB, \(k) {
    set.seed(SEED + 1)
    tryCatch(as.integer(trajClusters(m, nclusters = k)$partition[, 2]), error = \(e) rep(NA_integer_, length(idx)))
  })
  colnames(cls) <- paste0("k", K_STAB)
  out <- list(run = r, idx = idx, cls = cls, mins = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  saveRDS(out, paste0(r$file, ".tmp")); file.rename(paste0(r$file, ".tmp"), r$file)
  cat(sprintf("%s %s %s %s r%02d done (%.1f min)\n", format(Sys.time(), "%m-%d %H:%M"), r$variant, r$anchor, r$ptype, r$rep, out$mins))
}

if (length(args) >= 1 && args[1] == "run") {
  cores <- if (length(args) >= 2) as.integer(args[2]) else 1
  n_test <- if (length(args) >= 3) as.integer(args[3]) else nrow(RUNS)
  invisible(mclapply(seq_len(n_test), run_one, mc.cores = cores, mc.preschedule = FALSE))
  quit(save = "no")
}

# ---- Summary -------------------------------------------------------------------------------
ari <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) < 2) return(NA_real_)
  mclust::adjustedRandIndex(x[ok], y[ok])
}
fullv <- function(idx, v) { out <- rep(NA_integer_, N); out[idx] <- v; out }
bfiles <- list.files(BDIR, "^(raw|log)_(seed_r01|subsample|null).*\\.rds$", full.names = TRUE)
recs <- c(lapply(bfiles, \(f) { x <- readRDS(f); x$run$variant <- "base"; x }),
          lapply(RUNS$file[file.exists(RUNS$file)], readRDS))
sumf <- file.path(SUMM, "summary_traj_variants.txt"); unlink(sumf)
say  <- function(...) { cat(..., "\n"); cat(..., "\n", file = sumf, append = TRUE) }
pt   <- function(x) paste(capture.output(print(x, n = 300, width = 200)), collapse = "\n")
say(sprintf("traj variants (%s): %d of %d variant runs (+ %d base runs).", format(Sys.time(), "%Y-%m-%d %H:%M"),
            sum(file.exists(RUNS$file)), nrow(RUNS), length(bfiles)))

V <- bind_rows(lapply(recs, \(x) tibble(variant = x$run$variant, anchor = x$run$anchor, ptype = x$run$ptype, rep = x$run$rep,
                                        k = K_STAB, v = lapply(seq_along(K_STAB), \(j) fullv(x$idx, x$cls[, j])))))
pairs <- V |> filter(ptype != "seed") |> group_by(variant, anchor, ptype, k) |> group_modify(\(g, key) {
  if (nrow(g) < 2) return(tibble(ari = numeric()))
  ij <- t(combn(nrow(g), 2)); tibble(ari = apply(ij, 1, \(p) ari(g$v[[p[1]]], g$v[[p[2]]])))
}) |> ungroup()
ref <- V |> filter(ptype == "seed") |> select(variant, anchor, k, ref = v)
vsref <- V |> filter(ptype == "subsample") |> inner_join(ref, by = c("variant", "anchor", "k")) |>
  mutate(ari = map2_dbl(v, ref, ari))
sizes <- V |> filter(ptype != "null") |>
  mutate(min_share = map_dbl(v, \(x) { t <- table(x); min(t) / sum(t) }),
         n_tiny = map_int(v, \(x) { t <- table(x); sum(t / sum(t) < 0.01) })) |>
  group_by(variant, anchor, k) |>
  summarise(min_class_share = round(median(min_share), 3), share_runs_tiny_class = round(mean(n_tiny > 0), 2), .groups = "drop")
lvl <- c("base", names(VARIANTS))

say("\nSubsample stability: median pairwise ARI [IQR], median ARI vs the full-data fit, smallest class (median share),",
    "share of fits with a class under 1%:")
tab <- pairs |> filter(ptype == "subsample") |> group_by(variant, anchor, k) |>
  summarise(median = round(median(ari, na.rm = TRUE), 3), q25 = round(quantile(ari, .25, na.rm = TRUE), 3),
            q75 = round(quantile(ari, .75, na.rm = TRUE), 3), .groups = "drop") |>
  left_join(vsref |> group_by(variant, anchor, k) |> summarise(vs_full = round(median(ari, na.rm = TRUE), 3), .groups = "drop"),
            by = c("variant", "anchor", "k")) |>
  left_join(sizes, by = c("variant", "anchor", "k")) |>
  arrange(anchor, k, match(variant, lvl))
say(pt(tab))
say("\nNull data: median pairwise ARI between null replicates (high = traj finds the same structure in data without classes):")
say(pt(pairs |> filter(ptype == "null") |> group_by(variant, anchor, k) |> summarise(median = round(median(ari, na.rm = TRUE), 3), .groups = "drop") |>
         tidyr::pivot_wider(names_from = k, values_from = median, names_prefix = "k") |> arrange(anchor, match(variant, lvl))))
say("\nFull-data fits: ARI of each variant with base traj, same anchor and k:")
say(pt(ref |> inner_join(ref |> filter(variant == "base") |> select(anchor, k, b = ref), by = c("anchor", "k")) |>
         filter(variant != "base") |> mutate(ari = round(map2_dbl(ref, b, ari), 3)) |> select(variant, anchor, k, ari) |>
         tidyr::pivot_wider(names_from = k, values_from = ari, names_prefix = "k")))

st <- readRDS(file.path(SOUT, "summary.rds"))
other <- st$pairs |> filter(variant == "base", j > 0, ptype == "subsample") |>
  group_by(anchor, method, k) |> summarise(median = round(median(ari, na.rm = TRUE), 3), .groups = "drop")
cmp <- bind_rows(other, tab |> transmute(anchor, method = paste0("traj_", variant), k, median)) |>
  tidyr::pivot_wider(names_from = method, values_from = median) |> arrange(anchor, k)
say("\nMedian pairwise subsample ARI, traj variants next to the other methods:")
say(pt(cmp))
