# =============================================================================
# R/gmm_optima.R: how many random starts does GMM need to reach its best optimum?
#   Rscript R/gmm_optima.R <n_starts> <cores> [cells]
#   cells: comma-separated anchor:k, e.g. log:2,raw:6 (default: both anchors, k 2-6)
# Started 2026-10-26 after the stability run showed the main-run GMM fits were
# local optima at several k (README "Stability results").
#
# Each start is the GMM procedure with one start, run to the end: a 10-iteration
# hlme run from random(minit) (as lcmm::gridsearch does), then convergence from
# where that run ended. Recording both log-likelihoods for every start gives
#   - the share of starts reaching each optimum,
#   - P(gridsearch with R starts reaches the best): it converges only the start
#     with the best 10-iteration log-likelihood,
#   - P(converging all R starts reaches the best).
# Anchor data as in the stability analysis (one per family, ge3_imp, age - 9).
# One .rds per start in output/gmm_optima/ (restartable, atomic saves).
# =============================================================================

source("00_setup.R")
source("R/fit_task.R")
source("R/perturb.R")
suppressPackageStartupMessages(library(parallel))

args  <- commandArgs(trailingOnly = TRUE)
n_starts <- as.integer(args[1])
cores <- as.integer(args[2])
cells <- if (length(args) >= 3) strsplit(args[3], ",")[[1]] else
  paste0(rep(c("raw", "log"), each = 5), ":", 2:6)
dest <- file.path(OUT, "gmm_optima")
dir.create(dest, showWarnings = FALSE)

base <- readRDS(file.path(OUT, "stability", "anchor.rds"))
f <- lcmm_formulas(gmm = TRUE)

for (cell in cells) {
  anchor <- sub(":.*", "", cell); k <- as.integer(sub(".*:", "", cell))
  d <- make_d(list(idx = seq_along(base$ids), R = base$R, ages = base$ages), anchor, "gmm", base$ids)
  dl <- d$long; dl$id <- match(dl$id, d$ids)
  b <- list(fixed = f$fixed, random = f$random, subject = "id", data = dl, verbose = FALSE)
  # The 1-class model is fitted once per cell and bound globally before forking,
  # so every worker's random(.lcmm_minit) sees it (see lcmm_record in R/fit_task.R).
  assign(".lcmm_minit", do.call(lcmm::hlme, c(b, list(ng = 1))), envir = globalenv())

  out <- file.path(dest, sprintf("%s_k%d_s%03d.rds", anchor, k, seq_len(n_starts)))
  todo <- which(!file.exists(out))
  cat(sprintf("%s  %s k=%d: %d of %d starts to run\n", format(Sys.time(), "%H:%M"), anchor, k,
              length(todo), n_starts))
  res <- mclapply(todo, function(s) {
    set.seed(as.integer(SEED + 9e5 + 1e3 * k + s + if (anchor == "log") 5e4 else 0))
    t0 <- Sys.time()
    short <- do.call(lcmm::hlme, c(b, list(mixture = f$mixture, ng = k, maxiter = 10,
                                          B = quote(random(.lcmm_minit)))))
    full <- do.call(lcmm::hlme, c(b, list(mixture = f$mixture, ng = k, B = short$best)))
    p <- full$pprob
    cls <- as.integer(p$class[match(seq_along(d$ids), p$id)])
    rec <- list(anchor = anchor, k = k, start = s, ll_short = short$loglik,
                ll_full = full$loglik, conv = full$conv, bic = full$BIC, niter = full$niter,
                props = tabulate(cls, k) / length(cls), cls = cls,
                mins = as.numeric(difftime(Sys.time(), t0, units = "mins")))
    saveRDS(rec, paste0(out[s], ".tmp")); file.rename(paste0(out[s], ".tmp"), out[s])
    sprintf("%s s%03d %.1f min", cell, s, rec$mins)
  }, mc.cores = cores, mc.preschedule = FALSE)
  errs <- vapply(res, inherits, logical(1), "try-error")
  if (any(errs)) { cat("worker errors:\n"); print(res[errs]) }
}
