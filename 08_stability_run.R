# =============================================================================
# 08_stability_run.R: fit every stability replicate, in parallel, one .rds each.
#   Rscript 08_stability_run.R full|smoke [cores] [filters]
# filters: comma-separated key=value, keys anchor, method, ptype, maxrep. E.g.
#   Rscript 08_stability_run.R full 8 anchor=raw,maxrep=1    # timing pass
#   Rscript 08_stability_run.R full 8 method=gmm,ptype=site
# Replicates whose .rds already exists are skipped, so the run can be restarted;
# each one is saved under a temp name and renamed (same pattern as 03_run.R).
#
# A replicate record holds, per k:
#   fit    one row per part (all, or A / B for site splits): fit indices
#   parts  per part: idx (rows of the anchor) and classes, relabelled 1..K by
#          mean raw CBCL of the data the replicate clustered (the 04 rule)
#   pred   site splits: A's classes predicted for B's children (AonB) and the
#          reverse, in the predicting fit's relabelled classes
#   starts GBTM/GMM start recording: per-start log-likelihoods and classes
# =============================================================================

source("00_setup.R")
source("R/fit_task.R")
source("R/perturb.R")
source("R/predict_class.R")
suppressPackageStartupMessages(library(parallel))

args  <- commandArgs(trailingOnly = TRUE)
mode  <- if (length(args) >= 1) args[1] else "full"
cores <- if (length(args) >= 2) as.integer(args[2]) else 8
filt  <- if (length(args) >= 3) strsplit(strsplit(args[3], ",")[[1]], "=") else list()
stopifnot(mode %in% c("full", "smoke"))
if (mode == "smoke") N_STARTS <- 5
SOUT <- file.path(OUT, if (mode == "smoke") "stability_smoke" else "stability")

reps <- readRDS(file.path(SOUT, "spec_reps.rds"))
base <- readRDS(file.path(SOUT, "anchor.rds"))
nm   <- readRDS(file.path(SOUT, "null_model.rds"))
for (f in filt) {
  reps <- if (f[1] == "maxrep") filter(reps, rep <= as.integer(f[2]))
          else filter(reps, .data[[f[1]]] %in% strsplit(f[2], "\\|")[[1]])
}

dest <- file.path(SOUT, "reps")
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(LOGS, paste0("stability_", mode, ".log"))

# ---- One replicate ----------------------------------------------------------------
fit_rep <- function(rep) {
  parts <- rep_data(rep, base, nm)
  a <- rep$args[[1]]
  linkage <- if (is.null(a$linkage)) "ward.D2" else a$linkage
  fargs <- a[setdiff(names(a), "linkage")]
  ds <- lapply(parts, make_d, scale = rep$anchor, method = rep$method, ids_all = base$ids)
  trees <- if (rep$method == "ward") lapply(ds, \(d) hclust(dist(d$X), method = linkage))
  fitter <- if (rep$ptype == "starts") {
    function(d, k, seed) lcmm_record(d, k, seed, gmm = rep$method == "gmm", n_converge = N_CONVERGE)
  } else get(paste0("fit_", rep$method))

  lapply(K_STAB, function(k) {
    res <- list(); fit <- list(); out <- list()
    for (h in names(ds)) {
      t0 <- Sys.time()
      r <- tryCatch(
        if (rep$method == "ward") fitter(ds[[h]], k, rep$fit_seed, trees[[h]])
        else do.call(fitter, c(list(ds[[h]], k, rep$fit_seed), fargs)),
        error = function(e) list(error = conditionMessage(e)))
      row <- tibble(rep_id = rep$rep_id, k = k, part = h, n = length(ds[[h]]$ids),
                    seed = rep$fit_seed,
                    secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
      if (!is.null(r$error)) { fit[[h]] <- mutate(row, error = r$error); next }
      props <- tabulate(r$cls, k) / length(r$cls)
      fit[[h]] <- mutate(row, error = NA_character_, converged = r$converged,
                         bic = r$bic, aic = r$aic, entropy = r$entropy,
                         asw = asw(ds[[h]]$X, r$cls, rep$fit_seed),
                         k_found = sum(props > 0), min_prop = min(props),
                         small_class = min(props) < CLASS_FLOOR)
      r$map <- relabel_map(r$cls, parts[[h]]$R)
      res[[h]] <- r
      out[[h]] <- list(idx = parts[[h]]$idx, cls = unname(as.integer(r$map[as.character(r$cls)])))
    }
    pred <- NULL
    if (rep$ptype == "site" && all(c("A", "B") %in% names(res))) {
      p <- function(from, to) {
        cl <- tryCatch(predict_new(rep$method, res[[from]], ds[[from]], ds[[to]], fargs),
                       error = function(e) rep(NA_integer_, length(ds[[to]]$ids)))
        unname(as.integer(res[[from]]$map[as.character(cl)]))
      }
      pred <- list(AonB = p("A", "B"), BonA = p("B", "A"))
    }
    starts <- if (rep$ptype == "starts" && !is.null(res$all))
      res$all[c("starts", "start_cls", "full_cls")]
    list(fit = bind_rows(fit), parts = out, pred = pred, starts = starts)
  })
}

# ---- Run --------------------------------------------------------------------------
todo <- reps[!file.exists(file.path(dest, paste0(reps$rep_id, ".rds"))), ]
cat(sprintf("[%s] %d of %d replicates to fit on %d cores\n", mode, nrow(todo), nrow(reps), cores))

done <- mclapply(seq_len(nrow(todo)), function(i) {
  rep <- todo[i, ]
  t0 <- Sys.time()
  out <- fit_rep(rep)
  f <- file.path(dest, paste0(rep$rep_id, ".rds"))
  saveRDS(list(rep = rep, fits = out), paste0(f, ".tmp"))
  file.rename(paste0(f, ".tmp"), f)
  n_err <- sum(vapply(out, \(o) sum(!is.na(o$fit$error)), integer(1)))
  msg <- sprintf("%s  %-48s %7.1f min  %s", format(Sys.time(), "%H:%M:%S"), rep$rep_id,
                 as.numeric(difftime(Sys.time(), t0, units = "mins")),
                 if (n_err) paste(n_err, "fit errors") else "")
  cat(msg, "\n", file = logf, append = TRUE)
  msg
}, mc.cores = cores, mc.preschedule = FALSE)

errs <- vapply(done, inherits, logical(1), "try-error")
cat(sprintf("Finished: %d ok, %d worker errors\n", sum(!errs), sum(errs)))
if (any(errs)) print(done[errs])
writeLines(capture.output(sessionInfo()), file.path(LOGS, paste0("sessionInfo_stability_", mode, ".txt")))
