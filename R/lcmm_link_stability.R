# =============================================================================
# R/lcmm_link_stability.R: stability test of the skew-aware GBTM (lcmm with an
# estimated link function, no random effects), added 2026-11-13. Exploratory, not
# part of the pipeline yet. Follows R/traj_stability_test.R.
#   Rscript R/lcmm_link_stability.R run [cores] [n_test]   # restartable, one file per run
#   Rscript R/lcmm_link_stability.R summary
# Model as in R/lcmm_link_pilot.R (model "lcgm"): class-specific quadratic in age
# (centred at AGE_CENTER), link = "5-quant-splines", raw CBCL scores only (the link
# estimates the transformation, so there is no separate log anchor). k = 2-6.
# Starts: 1-class model, N_STARTS random starts of 10 iterations from it, each run to
# convergence; best converged start kept (the pilot found the same solution from
# every converged start at k = 2-5, so 3 starts are enough).
# Runs (same draws as the stability analysis, R/perturb.R):
#   subsample  20 80% subsamples (draws 1-20)
#   null       5 null datasets (draws 1-5; same as every other method, though their upper
#              tail is heavier than the real data: 99th percentile 131 vs 82)
# Reference: the pilot's full-data fits (output/lcmm_link/lcgm_k*.rds, 10 starts).
# Metric: ARI on shared children, all pairs of replicates within (perturbation, k),
# and each subsample vs the full-data fit; compared with GBTM on the raw anchor from
# output/stability/summary.rds.
# Output: output/lcmm_link/stability/*.rds; output/summaries/summary_lcmm_link_stability.txt
# =============================================================================

source("00_setup.R")
source("R/perturb.R")
suppressPackageStartupMessages({ library(lcmm); library(parallel) })

args <- commandArgs(trailingOnly = TRUE)
LDIR <- file.path(OUT, "lcmm_link")
SDIR <- file.path(LDIR, "stability"); dir.create(SDIR, showWarnings = FALSE, recursive = TRUE)
SOUT <- file.path(OUT, "stability")
base <- readRDS(file.path(SOUT, "anchor.rds"))
nm   <- readRDS(file.path(SOUT, "null_model.rds"))
N    <- length(base$ids)
N_STARTS <- 3

RUNS <- bind_rows(tibble(ptype = "subsample", rep = 1:20), tibble(ptype = "null", rep = 1:5)) |>
  mutate(file = file.path(SDIR, sprintf("%s_r%02d.rds", ptype, rep)))

fit_k <- function(common, m1, k, n) {
  runs <- lapply(seq_len(N_STARTS), function(s) tryCatch({
    sh <- do.call(lcmm, c(common, list(mixture = ~ time + I(time^2), ng = k, maxiter = 10, B = quote(random(.m1)))))
    m <- do.call(lcmm, c(common, list(mixture = ~ time + I(time^2), ng = k, B = sh$best)))
    list(m = m, ll = m$loglik, conv = m$conv)
  }, error = function(e) list(m = NULL, ll = NA_real_, conv = NA_real_)))
  ll <- vapply(runs, `[[`, numeric(1), "ll"); conv <- vapply(runs, `[[`, numeric(1), "conv")
  okc <- which(conv %in% 1)
  if (all(is.na(ll))) return(list(cls = rep(NA_integer_, n), ll = NA_real_, bic = NA_real_, n_conv = 0L, at_best = 0L))
  b <- if (length(okc)) okc[which.max(ll[okc])] else which.max(ll)
  m <- runs[[b]]$m; p <- m$pprob
  list(cls = as.integer(p$class[match(seq_len(n), p$id)]), ll = m$loglik, bic = m$BIC, conv = m$conv,
       n_conv = length(okc), at_best = if (length(okc)) sum(ll[okc] > ll[b] - 0.5) else 0L)
}

run_one <- function(i) {
  r <- RUNS[i, ]
  if (file.exists(r$file)) return(invisible(NULL))
  idx <- if (r$ptype == "subsample") draw_subsample(N, r$rep) else seq_len(N)
  R <- if (r$ptype == "null") simulate_null(base, nm, r$rep) else base$R
  X <- R[idx, , drop = FALSE]; A <- base$ages[idx, , drop = FALSE]; n <- length(idx)
  long <- data.frame(id = rep(seq_len(n), times = ncol(X)), time = as.vector(A) - AGE_CENTER, y = as.vector(X))
  long <- long[order(long$id, long$time), ]
  common <- list(fixed = y ~ time + I(time^2), random = ~ -1, subject = "id", data = long,
                 link = "5-quant-splines", verbose = FALSE)
  t0 <- Sys.time()
  set.seed(SEED + 700 + i)
  m1 <- do.call(lcmm, common); assign(".m1", m1, envir = globalenv())
  fits <- lapply(K_STAB, \(k) fit_k(common, m1, k, n))
  cls <- sapply(fits, `[[`, "cls"); colnames(cls) <- paste0("k", K_STAB)
  info <- tibble(k = K_STAB, ll = sapply(fits, `[[`, "ll"), bic = sapply(fits, `[[`, "bic"),
                 n_conv = sapply(fits, `[[`, "n_conv"), at_best = sapply(fits, `[[`, "at_best"))
  out <- list(run = r, idx = idx, cls = cls, info = info, mins = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  saveRDS(out, paste0(r$file, ".tmp")); file.rename(paste0(r$file, ".tmp"), r$file)
  cat(sprintf("%s %s r%02d done (%.1f min; converged starts k2-6: %s)\n", format(Sys.time(), "%m-%d %H:%M"),
              r$ptype, r$rep, out$mins, paste(info$n_conv, collapse = "/")))
}

if (length(args) >= 1 && args[1] == "run") {
  cores <- if (length(args) >= 2) as.integer(args[2]) else 1
  n_test <- if (length(args) >= 3) as.integer(args[3]) else nrow(RUNS)
  invisible(mclapply(seq_len(n_test), \(i) tryCatch(run_one(i), error = \(e) cat("run", i, "error:", conditionMessage(e), "\n")),
                     mc.cores = cores, mc.preschedule = FALSE))
  quit(save = "no")
}

# ---- Summary (aggregates only) ---------------------------------------------------------------
ari <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) < 2) return(NA_real_)
  mclust::adjustedRandIndex(x[ok], y[ok])
}
fullv <- function(idx, v) { out <- rep(NA_integer_, N); out[idx] <- v; out }
recs <- lapply(RUNS$file[file.exists(RUNS$file)], readRDS)
sumf <- file.path(SUMM, "summary_lcmm_link_stability.txt"); unlink(sumf)
say  <- function(...) { cat(..., "\n"); cat(..., "\n", file = sumf, append = TRUE) }
pt   <- function(x) paste(capture.output(print(x, n = 200, width = 200)), collapse = "\n")
say(sprintf("Skew-aware GBTM (lcmm, 5-quant-splines link) stability test (%s): %d of %d runs; median %.0f min per run (k 2-6, %d starts).",
            format(Sys.time(), "%Y-%m-%d %H:%M"), length(recs), nrow(RUNS), median(sapply(recs, `[[`, "mins")), N_STARTS))

V <- bind_rows(lapply(recs, \(x) tibble(ptype = x$run$ptype, rep = x$run$rep, k = K_STAB,
                                        v = lapply(seq_along(K_STAB), \(j) fullv(x$idx, x$cls[, j])))))
ref <- tibble(k = K_STAB, ref = lapply(K_STAB, \(k) readRDS(file.path(LDIR, sprintf("lcgm_k%d.rds", k)))$cls))
pairs <- V |> group_by(ptype, k) |> group_modify(\(g, key) {
  if (nrow(g) < 2) return(tibble(ari = numeric()))
  ij <- t(combn(nrow(g), 2)); tibble(ari = apply(ij, 1, \(p) ari(g$v[[p[1]]], g$v[[p[2]]])))
}) |> ungroup()
vsref <- V |> filter(ptype == "subsample") |> inner_join(ref, by = "k") |> mutate(ari = map2_dbl(v, ref, ari))
sizes <- V |> filter(ptype == "subsample") |>
  mutate(min_share = map_dbl(v, \(x) { t <- table(x); min(t) / sum(t) }),
         n_tiny = map_int(v, \(x) { t <- table(x); sum(t / sum(t) < 0.01) })) |>
  group_by(k) |> summarise(min_class_share = round(median(min_share), 3), share_runs_tiny_class = round(mean(n_tiny > 0), 2), .groups = "drop")
conv <- bind_rows(lapply(recs, \(x) mutate(x$info, ptype = x$run$ptype))) |> group_by(ptype, k) |>
  summarise(runs_no_converged_start = sum(n_conv == 0), median_at_best = median(at_best), .groups = "drop")

say("\nSubsample stability: median pairwise ARI [IQR], median ARI vs the full-data fit, smallest class, share of fits with a class under 1%:")
tab <- pairs |> filter(ptype == "subsample") |> group_by(k) |>
  summarise(median = round(median(ari, na.rm = TRUE), 3), q25 = round(quantile(ari, .25, na.rm = TRUE), 3),
            q75 = round(quantile(ari, .75, na.rm = TRUE), 3), .groups = "drop") |>
  left_join(vsref |> group_by(k) |> summarise(vs_full = round(median(ari, na.rm = TRUE), 3)), by = "k") |>
  left_join(sizes, by = "k")
say(pt(tab))
say("\nNull data: median pairwise ARI between null replicates (high = the model finds the same structure without classes):")
say(pt(pairs |> filter(ptype == "null") |> group_by(k) |> summarise(median = round(median(ari, na.rm = TRUE), 3))))
say("\nConvergence: runs where no start converged, median starts at the best solution (of 3):")
say(pt(conv))

st <- readRDS(file.path(SOUT, "summary.rds"))
other <- st$pairs |> filter(variant == "base", j > 0, ptype %in% c("subsample", "null"), anchor == "raw") |>
  group_by(method, ptype, k) |> summarise(median = round(median(ari, na.rm = TRUE), 3), .groups = "drop")
cmp <- bind_rows(other, pairs |> group_by(ptype, k) |> summarise(median = round(median(ari, na.rm = TRUE), 3), .groups = "drop") |>
                   mutate(method = "gbtm_link")) |>
  tidyr::pivot_wider(names_from = method, values_from = median) |> arrange(ptype, k)
say("\nMedian pairwise ARI, skew-aware GBTM next to the other methods on the raw anchor:")
say(pt(cmp))
