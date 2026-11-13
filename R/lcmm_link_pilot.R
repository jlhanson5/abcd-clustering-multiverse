# =============================================================================
# R/lcmm_link_pilot.R: skew-aware trajectory model, feasibility pilot (2026-11-12).
#   Rscript R/lcmm_link_pilot.R <model> <k> [n_starts]     # one fit; writes one file
#   Rscript R/lcmm_link_pilot.R summary
# lcmm::lcmm with an estimated link function: the observed CBCL score is a monotone
# transformation (I-splines, 5 knots at quantiles) of a latent Gaussian outcome, so
# the model estimates the transformation instead of assuming raw or log scores.
# Same class structure as GBTM / GMM in R/fit_task.R: class-specific quadratic in age
# (centred at AGE_CENTER); model "lcgm" has no random effects (GBTM analogue), "gmm"
# a random intercept and slope.
# Data: the stability raw anchor (one child per family, 4 waves, N = 8,331). Starts as
# in fit_gmm_conv: 1-class model, n_starts random starts of 10 iterations, best run to
# convergence (here every start is run to convergence, as for the GMM refit).
# Output: output/lcmm_link/<model>_k<k>.rds (classes, loglik, BIC, entropy, timings);
#         summary to output/summaries/summary_lcmm_link_pilot.txt (aggregates only)
# =============================================================================

source("00_setup.R")
suppressPackageStartupMessages(library(lcmm))
args <- commandArgs(trailingOnly = TRUE)
LDIR <- file.path(OUT, "lcmm_link"); dir.create(LDIR, showWarnings = FALSE)
base <- readRDS(file.path(OUT, "stability", "anchor.rds"))
N <- length(base$ids)
long <- data.frame(id = rep(seq_len(N), times = 4),
                   time = as.vector(base$ages) - AGE_CENTER,
                   y = as.vector(base$R))
long <- long[order(long$id, long$time), ]

if (length(args) >= 2 && args[1] != "summary") {
  model <- args[1]; k <- as.integer(args[2]); n_starts <- if (length(args) >= 3) as.integer(args[3]) else 10
  out_f <- file.path(LDIR, sprintf("%s_k%d.rds", model, k))
  rnd <- if (model == "gmm") ~ time else ~ -1
  common <- list(fixed = y ~ time + I(time^2), random = rnd, subject = "id", data = long,
                 link = "5-quant-splines", verbose = FALSE)
  set.seed(SEED + 500 + k)
  t0 <- Sys.time()
  m1 <- do.call(lcmm, common)
  t1 <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  if (k == 1) {
    res <- list(model = model, k = 1, m = m1, ll = m1$loglik, bic = m1$BIC, conv = m1$conv, mins = t1)
  } else {
    assign(".m1", m1, envir = globalenv())
    runs <- lapply(seq_len(n_starts), function(s) tryCatch({
      ts <- Sys.time()
      sh <- do.call(lcmm, c(common, list(mixture = ~ time + I(time^2), ng = k, maxiter = 10, B = quote(random(.m1)))))
      m <- do.call(lcmm, c(common, list(mixture = ~ time + I(time^2), ng = k, B = sh$best)))
      cat(sprintf("%s %s k%d start %d: ll %.1f conv %d (%.1f min)\n", format(Sys.time(), "%H:%M"), model, k, s,
                  m$loglik, m$conv, as.numeric(difftime(Sys.time(), ts, units = "mins"))))
      list(m = m, ll = m$loglik, conv = m$conv)
    }, error = function(e) { cat("start", s, "error:", conditionMessage(e), "\n"); list(m = NULL, ll = NA_real_, conv = NA_real_) }))
    ll <- vapply(runs, `[[`, numeric(1), "ll"); conv <- vapply(runs, `[[`, numeric(1), "conv")
    okc <- which(conv %in% 1)
    b <- if (length(okc)) okc[which.max(ll[okc])] else which.max(ll)
    m <- runs[[b]]$m; p <- m$pprob
    res <- list(model = model, k = k, cls = as.integer(p$class[match(seq_len(N), p$id)]),
                post = as.matrix(p[, -(1:2)]), ll = m$loglik, bic = m$BIC, conv = m$conv,
                ll_starts = ll, conv_starts = conv, n_conv = length(okc),
                at_best = if (length(okc)) sum(ll[okc] > ll[b] - 0.5) else 0L,
                link = list(est = m$estimlink), mins_1class = t1,
                mins = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  }
  saveRDS(res, paste0(out_f, ".tmp")); file.rename(paste0(out_f, ".tmp"), out_f)
  cat(sprintf("%s %s k%d done (%.1f min)\n", format(Sys.time(), "%m-%d %H:%M"), model, k, res$mins))
  quit(save = "no")
}

# ---- Summary (aggregates only) ---------------------------------------------------------------
sumf <- file.path(SUMM, "summary_lcmm_link_pilot.txt"); unlink(sumf)
say <- function(...) { cat(..., "\n"); cat(..., "\n", file = sumf, append = TRUE) }
fs <- list.files(LDIR, "\\.rds$", full.names = TRUE)
for (f in fs) {
  r <- readRDS(f)
  if (r$k == 1) { say(sprintf("%s k1: ll %.1f BIC %.1f conv %d (%.1f min)", r$model, r$ll, r$bic, r$conv, r$mins)); next }
  sh <- prop.table(table(r$cls)); mr <- tapply(rowMeans(base$R), r$cls, mean)
  ent <- 1 + sum(r$post * log(pmax(r$post, 1e-12))) / (N * log(r$k))
  say(sprintf("%s k%d: ll %.1f BIC %.1f; %d/%d starts converged, %d at best; %.0f min (1-class %.1f)",
              r$model, r$k, r$ll, r$bic, r$n_conv, length(r$ll_starts), r$at_best, r$mins, r$mins_1class))
  say(sprintf("   class shares %s; mean raw CBCL %s; entropy %.2f",
              paste(round(sh, 3), collapse = " / "), paste(round(mr, 1), collapse = " / "), ent))
}
