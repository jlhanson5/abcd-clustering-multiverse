# =============================================================================
# R/gmm_stab_refit.R: refit the GMM stability replicates at k = 2-4 with every start
# converged (fit_gmm_conv in R/fit_task.R), replacing the gridsearch fits.
#   Rscript --vanilla R/gmm_stab_refit.R [cores] [n_starts] [rep_regex] [reps_dir]
# Defaults: 8 cores, GMM_CONV_STARTS (R/perturb.R), every GMM replicate, the full
# run's output/stability/reps. Run after 08_stability_run.R has written every
# replicate. Decided 2026-11-07 (README "GMM refit"; the main-run GMM fits were
# often local optima, so replicate ARIs mixed data sensitivity with optimiser luck).
#
# - The "starts" replicates are left alone: they record gridsearch's own behaviour.
# - k = 5-6 keep the gridsearch fits: they are not reliably estimated (README).
# - Same data, fit seed and record layout as 08_stability_run.R, so 09 reads them
#   unchanged; fit rows gain fitter = "converged", n_conv and at_best.
# - Before a file is first replaced, its original is copied to reps_gridsearch_gmm/.
# - Restartable: a record already refitted (rec$gmm_conv) is skipped.
# =============================================================================

source("00_setup.R")
source("R/fit_task.R")
source("R/perturb.R")
source("R/predict_class.R")
suppressPackageStartupMessages(library(parallel))

args     <- commandArgs(trailingOnly = TRUE)
cores    <- if (length(args) >= 1) as.integer(args[1]) else 8
n_starts <- if (length(args) >= 2) as.integer(args[2]) else GMM_CONV_STARTS
rx       <- if (length(args) >= 3) args[3] else "."
SOUT     <- file.path(OUT, "stability")
dest     <- if (length(args) >= 4) args[4] else file.path(SOUT, "reps")
KS       <- 2:4
arch     <- file.path(dirname(dest), "reps_gridsearch_gmm")
dir.create(arch, showWarnings = FALSE)

reps <- readRDS(file.path(SOUT, "spec_reps.rds")) |>
  filter(method == "gmm", ptype != "starts", grepl(rx, rep_id))
base <- readRDS(file.path(SOUT, "anchor.rds"))
nm   <- readRDS(file.path(SOUT, "null_model.rds"))
logf <- file.path(LOGS, "gmm_stab_refit.log")

refit_rep <- function(rep) {
  f <- file.path(dest, paste0(rep$rep_id, ".rds"))
  if (!file.exists(f)) return(sprintf("%s: missing, run 08 first", rep$rep_id))
  rec <- readRDS(f)
  if (!is.null(rec$gmm_conv)) return(sprintf("%s: already refitted", rep$rep_id))
  a <- file.path(arch, basename(f)); if (!file.exists(a)) file.copy(f, a)
  t0 <- Sys.time()
  parts <- rep_data(rep, base, nm)
  fargs <- rep$args[[1]]
  ds <- lapply(parts, make_d, scale = rep$anchor, method = "gmm", ids_all = base$ids)
  for (k in KS) {
    res <- list(); fit <- list(); out <- list()
    for (h in names(ds)) {
      t1 <- Sys.time()
      r <- tryCatch(do.call(fit_gmm_conv, c(list(ds[[h]], k, rep$fit_seed, n_starts), fargs)),
                    error = function(e) list(error = conditionMessage(e)))
      row <- tibble(rep_id = rep$rep_id, k = k, part = h, n = length(ds[[h]]$ids), seed = rep$fit_seed,
                    secs = as.numeric(difftime(Sys.time(), t1, units = "secs")), fitter = "converged")
      if (!is.null(r$error)) { fit[[h]] <- mutate(row, error = r$error); next }
      props <- tabulate(r$cls, k) / length(r$cls)
      fit[[h]] <- mutate(row, error = NA_character_, converged = r$converged, bic = r$bic, aic = r$aic,
                         entropy = r$entropy, asw = asw(ds[[h]]$X, r$cls, rep$fit_seed),
                         k_found = sum(props > 0), min_prop = min(props),
                         small_class = min(props) < CLASS_FLOOR, n_conv = r$n_conv, at_best = r$at_best)
      r$map <- relabel_map(r$cls, parts[[h]]$R)
      res[[h]] <- r
      out[[h]] <- list(idx = parts[[h]]$idx, cls = unname(as.integer(r$map[as.character(r$cls)])))
    }
    pred <- NULL
    if (rep$ptype == "site" && all(c("A", "B") %in% names(res))) {
      p <- function(from, to) {
        cl <- tryCatch(predict_new("gmm", res[[from]], ds[[from]], ds[[to]], fargs),
                       error = function(e) rep(NA_integer_, length(ds[[to]]$ids)))
        unname(as.integer(res[[from]]$map[as.character(cl)]))
      }
      pred <- list(AonB = p("A", "B"), BonA = p("B", "A"))
    }
    rec$fits[[match(k, K_STAB)]] <- list(fit = bind_rows(fit), parts = out, pred = pred, starts = NULL)
  }
  rec$gmm_conv <- list(ks = KS, n_starts = n_starts, date = Sys.time())
  saveRDS(rec, paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f)
  msg <- sprintf("%s  %-40s %6.1f min", format(Sys.time(), "%m-%d %H:%M"), rep$rep_id,
                 as.numeric(difftime(Sys.time(), t0, units = "mins")))
  cat(msg, "\n", file = logf, append = TRUE)
  msg
}

# Full-sample replicates first (slowest), site halves last.
reps <- reps[order(match(reps$ptype, c("noise90", "noise80", "null", "setting", "subsample", "site"))), ]
cat(sprintf("%d GMM replicates, k = %s, %d starts each, %d cores\n", nrow(reps),
            paste(KS, collapse = "-"), n_starts, cores))
done <- mclapply(seq_len(nrow(reps)), \(i) refit_rep(reps[i, ]), mc.cores = cores, mc.preschedule = FALSE)
errs <- vapply(done, inherits, logical(1), "try-error")
cat(sprintf("Finished: %d ok, %d worker errors\n", sum(!errs), sum(errs)))
if (any(errs)) print(done[errs])
